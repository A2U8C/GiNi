import sys
import os
import subprocess
import argparse
import pandas as pd
import nipype
import nipype.pipeline.engine as pe
import nipype.interfaces.utility as niu
from nipype.interfaces.utility import Function

#


################################Review
'''

##Make the temporary files and directories outside the package
###Job Scheduler https://nipype.readthedocs.io/en/0.11.0/users/plugins.html
### Put a flag for the type of job scheduler

If there is any error while executing, make sure thepackage mentions the error

hqw (In qsub) (How to hold jobs for the future ones)

How to jump directly to a step and contoinue the process from there. (Make the path ending with unique string to identify the process till where it was executed.)

Preparation of qced files: Use python dataframe (replace the column names )

Do check "call" for the status of a specific job
"check_call" 
'''



'''
nipypecli crash .pklz

'''

######################################



sys.path.append('/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/')

from gini_main import * ##file_checker
from LDSC_Module.ldsc_manager import * ##General_Munge

metal_cmd_path="/ifs/loni/faculty/njahansh/nerds/ravi/genetics/random-metal-0.1.0/executables/metal"


   
main_path=""

def metal_checker(studies_list)->bool:
    set_all_traits=set() # reference traits list; To cross check that trait_list are same for all the studies
    
    trait_study_dict={}#Create a dict for all traits, trait_i->[list of traits_paths for all studies ]

    #change it to the "-"; Mention  in the documentation
    sep_btw_study_trait="___" #Seperator between study_name___trait_name.regenie 


    #Get list of traits files from 1st study for reference
    # print(studies_list)
    with open(studies_list[0], 'r') as file:
        trait_files_study_1 = file.read().splitlines()
    
    #Prepare the reference dict from study 1
    for trait_i_study_1 in trait_files_study_1:
        file_checker(trait_i_study_1)# To check if the file path exists
        trait_i=trait_i_study_1.split("/")[-1].split(".")[0].split(sep_btw_study_trait)[-1]
        set_all_traits.add(trait_i) #Preparing reference set of traits
        trait_study_dict[trait_i]=[trait_i_study_1]
    
    #Check if for each study, the traits are same
    for study_i in studies_list[1:]:
        set_all_traits_i=set() #set for studi_i; later checked against the reference set_all_traits
        file_checker(study_i)# To check if the file path exists
        with open(study_i, 'r') as file:
            trait_files_study_i = file.read().splitlines()

        for trait_j_study_i in trait_files_study_i:
            trait_j=trait_j_study_i.split("/")[-1].split(".")[0].split(sep_btw_study_trait)[-1]
            set_all_traits_i.add(trait_j)

            if trait_j not in trait_study_dict:
                print("Different traits in different studies, Make sure the traits across the studies are same.")
                ######Mention the file which gives the error
                sys.exit()
            trait_study_dict[trait_j].append(trait_j_study_i)

      
        if set_all_traits!=set_all_traits_i:
            print("Different traits in different studies, Make sure the traits across the studies are same.")
            sys.exit()

    return trait_study_dict
        

         
    pass



def metal_execution_function(trait_name:str,trait_study_loc_paths_str:str,meta_type:str):
    import os
    import subprocess
    
    import sys  
    
    metal_cmd_path="/ifs/loni/faculty/njahansh/nerds/ravi/genetics/random-metal-0.1.0/executables/metal"
    metal_file_content="""#Metal script for running a random effects meta-analysis from REGENIE summary stats

# Options for describing input files ...
   MARKERLABEL      ID
   ALLELELABELS     ALLELE1 ALLELE0
   EFFECTLABEL      BETA
# Options for inverse variance weighted meta-analysis ...
   STDERRLABEL      SE
   SCHEME           STDERR
# Options to enable tracking of allele frequencies ...
   AVERAGEFREQ      ON
   MINMAXFREQ       ON
   FREQLABEL        A1FREQ
# Options to enable tracking of user defined variables ...
   CUSTOMVARIABLE   TotalN
   LABEL            TotalN AS N
# Automatic genomic control correction of input statistics ...
   GENOMICCONTROL   ON
# Specify the names of the QCed data files
{PROCESSFILE_contents}
# Options for general analysis control ...
   EFFECT_PRINT_PRECISION 6
   STDERR_PRINT_PRECISION 6
   OUTFILE                {abs_meta_path}{out_name} RE_{out_name}.tbl
   ANALYZE                {meta_type}

"""


    process_file_i=""
    qced_files=[]
    
    trait_study_loc_paths=trait_study_loc_paths_str.split("abcd123__")

    print("*-*-*--**-*-*---**-*-*-**-*--*-*-*--*-*-**-*========================================================================")

    extras_metal_qced_files="/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/Exta_temp_files/METAL_Out/METAL_QCED_file/"#os.path.join(os.getcwd(), "Exta_temp_files/METAL_Out/METAL_QCED_file/")
    extras_metal_input_files="/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/Exta_temp_files/METAL_Out/METAL_input_files/"#os.path.join(os.getcwd(), "Exta_temp_files/METAL_Out/METAL_input_files/")
    extras_metal_output_files="/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/Exta_temp_files/METAL_Out/METAL_output_files/"#os.path.join(os.getcwd(), "Exta_temp_files/METAL_Out/METAL_output_files/")
    for study_i in trait_study_loc_paths: 
        study_i_f_name=study_i.split("/")[-1].split(".")[0] #study_i__trait name

        ###Big files are created in the intermediate step
        ## Change the hard coded $6 column ($6 column is MAF value)
        ## Change the column names according to the input format
        study_i_prep_qced_cmd=f'''echo "CHROM GENPOS ID ALLELE0 ALLELE1 A1FREQ N TEST BETA SE CHISQ LOG10P EXTRA P" > {extras_metal_qced_files+study_i_f_name}.qced'''
        study_i_final_qced_cmd=f"""awk '{{ if ($6>=0.01 && $6<=.99) print $0}}' {study_i} >> {extras_metal_qced_files+study_i_f_name}.qced"""
        study_i_gzip_qced_cmd=f"""gzip -9 {extras_metal_qced_files+study_i_f_name}.qced"""

        subprocess.call(study_i_prep_qced_cmd, shell=True)
        subprocess.call(study_i_final_qced_cmd, shell=True)
        subprocess.call(study_i_gzip_qced_cmd, shell=True)
        qced_files.append(extras_metal_qced_files+study_i_f_name+".qced.gz")
        #######################################

        process_file_i+=f"PROCESSFILE\t{extras_metal_qced_files+study_i_f_name}.qced.gz \n"
    


    # f=base_f.format(PROCESSFILE_contents=process_file_i,abs_meta_path=extras_metal_output_files,out_name=trait_i)
    f=metal_file_content.format(PROCESSFILE_contents=process_file_i,abs_meta_path=extras_metal_output_files,meta_type=meta_type,out_name=trait_name)
    out_file_metal=extras_metal_input_files+"/RE_"+trait_name+".metal"


    with open(out_file_metal, 'w') as n_f:
        n_f.write(f)


    metal_cmd=f"""{metal_cmd_path} < """+out_file_metal+f""" > {extras_metal_output_files}/RE.metal.{trait_name}.out"""

    
    subprocess.call(metal_cmd, shell=True)
    


    # for qced_file_i in qced_files:
    #     os.remove(qced_file_i) # To remove the temporary files created during the METAL input preparation

    new_tbl_file=extras_metal_output_files+trait_name+"1RE_"+trait_name+".tbl"
    gzip_tbl_file=f"""gzip -9 {new_tbl_file}"""
    subprocess.call(gzip_tbl_file, shell=True)
    out_zipped_metal=new_tbl_file+".gz"

    return out_zipped_metal
    






def metal_improved_execution_function(trait_name:str,trait_study_loc_paths_str:str,meta_type:str):
    # import time

    # time.sleep(5)
    # print(trait_name,meta_type)
    # return trait_study_loc_paths_str

    import os
    import subprocess
    import json
    import sys  
    import pandas as pd
    with open('/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/association_format.json', 'r') as json_file:
        json_association = json.load(json_file)

    
    metal_cmd_path="/ifs/loni/faculty/njahansh/nerds/ravi/genetics/random-metal-0.1.0/executables/metal"
    metal_file_content="""#Metal script for running a random effects meta-analysis from REGENIE summary stats

# Options for describing input files ...
   MARKERLABEL      ID
   ALLELELABELS     ALLELE1 ALLELE0
   EFFECTLABEL      BETA
# Options for inverse variance weighted meta-analysis ...
   STDERRLABEL      SE
   SCHEME           STDERR
# Options to enable tracking of allele frequencies ...
   AVERAGEFREQ      ON
   MINMAXFREQ       ON
   FREQLABEL        A1FREQ
# Options to enable tracking of user defined variables ...
   CUSTOMVARIABLE   TotalN
   LABEL            TotalN AS N
# Automatic genomic control correction of input statistics ...
   GENOMICCONTROL   ON
# Specify the names of the QCed data files
{PROCESSFILE_contents}
# Options for general analysis control ...
   EFFECT_PRINT_PRECISION 6
   STDERR_PRINT_PRECISION 6
   OUTFILE                {abs_meta_path}{out_name} RE_{out_name}.tbl
   ANALYZE                {meta_type}

"""


    process_file_i=""
    qced_files=[]
    
    trait_study_loc_paths=trait_study_loc_paths_str.split("abcd123__")

    print("*-*-*--**-*-*---**-*-*-**-*--*-*-*--*-*-**-*========================================================================")

    extras_metal_qced_files="/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/Exta_temp_files/METAL_Out/METAL_QCED_file/"#os.path.join(os.getcwd(), "Exta_temp_files/METAL_Out/METAL_QCED_file/")
    extras_metal_input_files="/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/Exta_temp_files/METAL_Out/METAL_input_files/"#os.path.join(os.getcwd(), "Exta_temp_files/METAL_Out/METAL_input_files/")
    extras_metal_output_files="/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/Exta_temp_files/METAL_Out/METAL_output_files/"#os.path.join(os.getcwd(), "Exta_temp_files/METAL_Out/METAL_output_files/")
    for study_i in trait_study_loc_paths: 

        compl_split_study_i=study_i.split("/")[-1].split(".")
        study_i_f_name=compl_split_study_i[0] #study_i__trait name



        gwas_format_name=compl_split_study_i[-1] if compl_split_study_i[-1].lower()!="gz" else compl_split_study_i[-2]
        dict_col_header_mapper=json_association[gwas_format_name]['format_dict']


        df_study_i=pd.read_csv(study_i,sep=" ")
        df_study_i.rename(columns=dict_col_header_mapper,inplace=True)

        if gwas_format_name=="regenie":    
            df_study_i["P"]= (1/10)**df_study_i['LOG10P']

        if "EXTRA" in df_study_i.columns:
            df_study_i["EXTRA"]=df_study_i["EXTRA"].fillna("NA")

        # print(df_study_i.columns)
        df_study_i=df_study_i[(df_study_i["A1FREQ"]>=0.01) & (df_study_i["A1FREQ"]<=.99)]
        df_study_i[["CHROM","GENPOS","ID","ALLELE0","ALLELE1","A1FREQ","N","TEST","BETA","SE","CHISQ","LOG10P","EXTRA","P"]].to_csv(extras_metal_qced_files+study_i_f_name+".qced",sep="\t",index=False)

        study_i_gzip_qced_cmd=f"""gzip -9 {extras_metal_qced_files+study_i_f_name}.qced"""

        subprocess.call(study_i_gzip_qced_cmd, shell=True)
        qced_files.append(extras_metal_qced_files+study_i_f_name+".qced.gz")
        #######################################

        process_file_i+=f"PROCESSFILE\t{extras_metal_qced_files+study_i_f_name}.qced.gz \n"
    


    # f=base_f.format(PROCESSFILE_contents=process_file_i,abs_meta_path=extras_metal_output_files,out_name=trait_i)
    f=metal_file_content.format(PROCESSFILE_contents=process_file_i,abs_meta_path=extras_metal_output_files,meta_type=meta_type,out_name=trait_name)
    out_file_metal=extras_metal_input_files+"/RE_"+trait_name+".metal"


    with open(out_file_metal, 'w') as n_f:
        n_f.write(f)


    metal_cmd=f"""{metal_cmd_path} < """+out_file_metal+f""" > {extras_metal_output_files}/RE.metal.{trait_name}.out"""

    
    subprocess.call(metal_cmd, shell=True)
    


    for qced_file_i in qced_files:
        os.remove(qced_file_i) # To remove the temporary files created during the METAL input preparation

    new_tbl_file=extras_metal_output_files+trait_name+"1RE_"+trait_name+".tbl"
    gzip_tbl_file=f"""gzip -9 {new_tbl_file}"""
    subprocess.call(gzip_tbl_file, shell=True)
    out_zipped_metal=new_tbl_file+".gz"

    return out_zipped_metal
    





def metal_main_program(trait_studies_input,hpc_bool:bool,meta_type:str,machineName:str):


    extras_folder_path = os.path.join(os.getcwd(), "Exta_temp_files")
    if not (os.path.exists(extras_folder_path) and os.path.isdir(extras_folder_path)): # Check if folder exists; to add temporary folder in it
        os.makedirs(extras_folder_path)
    
    extras_metal_files = os.path.join(os.getcwd(), "Exta_temp_files/METAL_Out")
    if not (os.path.exists(extras_metal_files) and os.path.isdir(extras_metal_files)): # Check if Metal extras folder exists; to add temporary files in it
        os.makedirs(extras_metal_files)

    extras_metal_qced_files = os.path.join(os.getcwd(), "Exta_temp_files/METAL_Out/METAL_QCED_file")
    if not (os.path.exists(extras_metal_qced_files) and os.path.isdir(extras_metal_qced_files)): # Check if Metal extras folder exists; to add temporary files in it
        os.makedirs(extras_metal_qced_files)
        
    extras_metal_input_files = os.path.join(os.getcwd(), "Exta_temp_files/METAL_Out/METAL_input_files")
    if not (os.path.exists(extras_metal_input_files) and os.path.isdir(extras_metal_input_files)): # Check if Metal extras folder exists; to add temporary files in it
        os.makedirs(extras_metal_input_files)

    extras_metal_output_files = os.path.join(os.getcwd(), "Exta_temp_files/METAL_Out/METAL_output_files")
    if not (os.path.exists(extras_metal_output_files) and os.path.isdir(extras_metal_output_files)): # Check if Metal extras folder exists; to add temporary files in it
        os.makedirs(extras_metal_output_files)

    extras_LDSC_Munge_files = os.path.join(os.getcwd(), "Exta_temp_files/LDSC_Out/Munged_results")
    if not (os.path.exists(extras_LDSC_Munge_files) and os.path.isdir(extras_LDSC_Munge_files)): # Check if LDSC extras folder exists; to add temporary files in it
        os.makedirs(extras_LDSC_Munge_files)

    extras_LDSC_Heri_files = os.path.join(os.getcwd(), "Exta_temp_files/LDSC_Out/Heritability_Results")
    if not (os.path.exists(extras_LDSC_Heri_files) and os.path.isdir(extras_LDSC_Heri_files)): # Check if LDSC extras folder exists; to add temporary files in it
        os.makedirs(extras_LDSC_Heri_files)

    metal_output_file_paths=[]

    if hpc_bool:
        with open(extras_metal_input_files+"/metal_hpc_trait_study_locs.txt", 'w') as n_f:
            for trait_i_name in trait_studies_input:
                n_f.write(f'{trait_i_name}:{",".join(trait_studies_input[trait_i_name])}\n')
                metal_output_file_paths.append(extras_metal_output_files+trait_i_name+"1RE_"+trait_i_name+".tbl.gz")
            

        traits_cnt=len(trait_studies_input.keys())
        # {'trait_1':[path_1_study1,path_1_study2]}
        ##runnow.q
        # template_el=f'{main_path}/HPC_Module/sge_file_prep.sh'
        # qsub_args_el=f'-q {machineName} -t 1:{traits_cnt} {main_path}/HPC_Module/sge_file_prep.sh Meta_Analysis {extras_metal_input_files}/metal_hpc_trait_study_locs.txt {meta_type}'
        # nipype.pipeline.engine.workflows.Workflow.run(plugin='SGE',plugin_args=dict(template=template_el,qsub_args=qsub_args_el))
         
        # qsub_cmd=f"qsub -V -q runnow.q -t 1:{traits_cnt} {main_path}/HPC_Module/sge_file_prep.sh Meta_Analysis {extras_metal_input_files}/metal_hpc_trait_study_locs.txt {meta_type}"
        # subprocess.call(qsub_cmd, shell=True)

        all_trait_studies_input_items=trait_studies_input.items()

        all_trait_studies_input_keys=[k for k,v in all_trait_studies_input_items]#[:1] #[trait_1,trait_2.....trait_n][::-1]
        all_trait_studies_input_values=["abcd123__".join(v) for k,v in all_trait_studies_input_items]#[:1]  #[ukb_trait_1abcd123__abcd_trait_1,ukb_trait_2abcd123__abcd_trait_2.....ukb_trait_nabcd123__abcd_trait_n]
        meta_full=[meta_type for i in all_trait_studies_input_values]#[:1] #["RANDOM","RANDOM",....."RANDOM"]
        

        workflow_metal_mungng = pe.Workflow(name='METAL_Munging_wf2')
        workflow_metal_mungng.base_dir = "/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/Exta_temp_files/"



        iterable_new_list=[[all_trait_studies_input_keys[i],all_trait_studies_input_values[i],meta_full[i]] for i in range(len(meta_full))]#len(meta_full)
        
        
        inputnode = pe.Node(
            niu.IdentityInterface(
                fields=['trait_name', 'trait_study_loc_paths_str', 'meta_type']
                ),
            synchronize=True,
            iterables=[('trait_name','trait_study_loc_paths_str','meta_type'),iterable_new_list],
            name='inputnode'
        )
        

        node1_METAL_Results = pe.Node(Function(input_names=['trait_name', 'trait_study_loc_paths_str', 'meta_type'],
                      output_names=['out_zipped_metal'],
                      function=metal_improved_execution_function),
             name='node1_METAL_Results')


        node2_LDSC_Munge = pe.Node(Function(input_names=['filePathInp','kwargs'],
                      output_names=['LDSC_Munge_out'],
                      function=General_Munge),
             name='node2_LDSC_Munge')
        node2_LDSC_Munge.inputs.kwargs=dict(N_col="TotalN",snp="MarkerName",frq="Freq1",signed_sumstats="EffectARE",out=extras_LDSC_Munge_files+"/",a1="Allele1",a2="Allele2",p="PvalueARE")
        
        
        node3_LDSC_Heritability = pe.Node(Function(input_names=['filePathInp_Munged'],
                      output_names=['LDSC_Heritability_out'],
                      function=HeritabilityLDSC),
             name='node3_LDSC_Heritability')
        # node3_LDSC_Heritability.inputs.kwargs_Munged=dict(out=extras_LDSC_Heri_files+"/")
        

        connections=[(inputnode, node1_METAL_Results, [('trait_name', 'trait_name')]),
                     (inputnode, node1_METAL_Results, [('trait_study_loc_paths_str', 'trait_study_loc_paths_str')]),
                     (inputnode, node1_METAL_Results, [('meta_type', 'meta_type')]),
                     (node1_METAL_Results, node2_LDSC_Munge, [('out_zipped_metal', 'filePathInp')]), #out_zipped_metal is the output from node1_METAL_Results
                     (node2_LDSC_Munge, node3_LDSC_Heritability, [('LDSC_Munge_out', 'filePathInp_Munged')])] 
        
        workflow_metal_mungng.connect(connections) 
        
        workflow_metal_mungng.write_graph(graph2use='orig', dotfilename='/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/graph_orig.dot')
       
        workflow_metal_mungng.run('SGE',plugin_args=dict(dont_resubmit_completed_jobs= True, overwrite= True, qsub_args=f'-q {machineName} -o /ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/Exta_temp_files/Logs/ -j y'))

    else:
        for trait_i in trait_studies_input:# HPC benefits (Put a qsub job on each trait)
            metal_file_gz=metal_improved_execution_function(trait_i,"abcd123__".join(trait_studies_input[trait_i]),meta_type)
            General_Munge(metal_file_gz)
            # metal_output_file_paths.append(extras_metal_output_files+trait_i_name+"1RE_"+trait_i_name+".tbl.gz")
        os.removedirs(extras_metal_qced_files) # To delete the QCED temporary files
    
    file_metal_results_path=extras_metal_input_files+"metal_results_gz_files.txt"
    # with open(file_metal_results_path, 'w') as met_out_paths:
    #     for line in metal_output_file_paths:
    #         met_out_paths.write(line+"\n")
    print("Metal Completed")
    return "file_metal_results_path"


def metal_Analysis_Module(studies_list,meta_type):
    print("META Analysis Module")
    dict_trait_pathStudy=metal_checker(studies_list)
    return metal_main_program(dict_trait_pathStudy,hpc_bool=True,meta_type=meta_type,machineName="runnow.q")
    