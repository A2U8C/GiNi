
def LAVA_TWAS_run_function(tissue_sqtl:str,ref_file:str,TWAS_case_control:str,chr_i:str):

    import os
    import sys
    import subprocess
        
    sys.path.append('/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/')
    import CONSTANTS
     
    # r_command=f'Rscript /ifs/loni/faculty/njahansh/GAMBIT/genome_WebApplication/static/LAVA/LAVA_TWAS_Grid_sQTL.R "/ifs/loni/faculty/njahansh/nerds/ravi/genetics/lava-mdd-endo-2023/data/g1000_eur" "/ifs/loni/faculty/njahansh/nerds/ravi/genetics/LAVA-main/LAVA_info.txt" ${{tissue_type}} "/ifs/loni/faculty/njahansh/nerds/ravi/genetics/LAVA-main/sqtl" ${{chr}} "LAVA_TWAS_sQTL_${{tissue_type}}_${{chr}}"'

    output_file_sqtl=CONSTANTS.Extra_temp_files_dict["extra_LAVA_TWAS_Results"]+"/"+tissue_sqtl+"/"+"LAVA_TWAS_sQTL_"+tissue_sqtl+"_"+str(chr_i)
    r_command=f'Rscript /ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/LAVA_TWAS_Module/LAVA_TWAS_Grid_sQTL.R "{ref_file}" "{TWAS_case_control}" "{tissue_sqtl}" "/ifs/loni/faculty/njahansh/nerds/ravi/genetics/LAVA-main/sqtl" {chr_i} "{output_file_sqtl}"'
    print(r_command)

    # Define the environment variables specific to your R environment
    r_environment_path = "/ifs/loni/faculty/njahansh/nerds/shruti/software/python/anaconda3/envs/r_env/bin"

    r_environment = os.environ.copy()
    r_environment["PATH"] = f"{r_environment_path}:{r_environment['PATH']}"

    # Execute the command
    subprocess.run(r_command, env=r_environment,shell=True)

    return output_folder


def Lava_TWAS_input_file(study_files:list,TWAS_case_control):
    import pandas as pd
    import sys
    import time
    sys.path.append('/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/')
    import CONSTANTS 

    def func_Pheno_Ext(row):
        return row.split("/")[-1].split(".")[0].split(CONSTANTS.sep_btw_study_trait)[1].lstrip("_")

    
    if TWAS_case_control==None:
        dataframe_data=pd.DataFrame(columns=["phenotype","cases","filename"])
        for el in study_files:
            data_i={"phenotype":el.split("/")[-1].split(".")[0].split(CONSTANTS.sep_btw_study_trait)[1].lstrip("_"),"cases":"NA","controls":"NA","filename":el} 
            new_row_df = pd.DataFrame([data_i])
            dataframe_data = pd.concat([dataframe_data, new_row_df], ignore_index=True)
        df_case_controls=pd.DataFrame(data=dataframe_data)
    else:
        df_case_controls=pd.read_csv(TWAS_case_control,sep="\t")
        df_case_controls["phenotype"]=df_case_controls["file_location"].apply(func_Pheno_Ext) 
        df_case_controls.rename(columns={"file_location":"filename","Cases":"cases","Controls":"controls"},inplace=True)
    
    f_name=CONSTANTS.Extra_temp_files_dict["extra_LAVA_TWAS_CaseControl"]+"/"+"LAVA_TWAS_"+str(time.time()).split(".")[0]+".txt"
    print(df_case_controls)
    df_case_controls[["phenotype","cases","controls","filename"]].to_csv(f_name, na_rep="NA",sep="\t",index=False)


    LAVA_input_out=f_name
    return LAVA_input_out


def LAVA_TWAS_script(ref_file:str,study_files:list,TWAS_case_control,machineName:str="runnow.q") -> None:
    import os
    import sys
    import time
    sys.path.append('/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/')
    import CONSTANTS

    import nipype
    import nipype.pipeline.engine as pe
    import nipype.interfaces.utility as niu
    from nipype.interfaces.utility import Function
    import shutil
    from LAVA_TWAS_Module import LAVA_TWAS_manager
    from nipype import config

    LAVA_TWAS_CaseControl=Lava_TWAS_input_file(study_files,TWAS_case_control)
    wf_name='LAVA_TWAS_exec_'+str(time.time()).split(".")[0]
    workflow_metal_mungng = pe.Workflow(name=wf_name)
    workflow_metal_mungng.base_dir = CONSTANTS.Extra_temp_files_dict["extra_Nipype_Wf"]+"/" #"/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/Exta_temp_files/"
    workflow_metal_mungng.config["job_finished_timeout"]="60"

    
    config_dict={'execution': {'job_finished_timeout': '60'}}
    config.update_config(config_dict)
    # output_folder="/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/Exta_temp_files/LAVA_Out/LAVA_Results"+"/"+trait_name+"/"+wSA_wTHICK

    all_tissues=["brain_amygdala","brain_anterior_cingulate_ba24","brain_caudate","brain_cerebellar","brain_cerebellum","brain_cortex","brain_frontal_cortex_ba9","brain_hippocampus","brain_hypothalamus","brain_nucleus_accumbens","brain_putamen","brain_spinal_cord","brain_substantia_nigra","cells_fibroblasts","cells_lymphocytes","whole_blood"]


    iterable_new_list=[]
    for tissue_sqtl in all_tissues:
        output_file_sqtl=CONSTANTS.Extra_temp_files_dict["extra_LAVA_TWAS_Results"]+"/"+tissue_sqtl
        if not os.path.exists(output_file_sqtl):
            os.makedirs(output_file_sqtl)
            print("Nested folders created successfully.")

        iterable_new_list+=[[tissue_sqtl,chr_el_i] for chr_el_i in range(1,23)]

    print(iterable_new_list)

    # iterable_new_list=list(range(1,18381))
    inputnode = pe.Node(
        niu.IdentityInterface(
            fields=['tissue_sqtl_input','chr_i_input']
            ),
            synchronize=True,
            iterables=[('tissue_sqtl_input','chr_i_input'),iterable_new_list],
            name='inputnode'
            )
        
    node1_parallel_LAVA_TWAS = pe.Node(Function(input_names=['tissue_sqtl','ref_file','TWAS_case_control','chr_i'],
                    output_names=['out_file_comp'],
                    function=LAVA_TWAS_manager.LAVA_TWAS_run_function),
            name='node1_parallel_LAVA_TWAS')
    node1_parallel_LAVA_TWAS.inputs.ref_file=ref_file
    node1_parallel_LAVA_TWAS.inputs.TWAS_case_control=LAVA_TWAS_CaseControl

    
    connections=[
                    (inputnode, node1_parallel_LAVA_TWAS, [('tissue_sqtl_input', 'tissue_sqtl')]),
                    (inputnode, node1_parallel_LAVA_TWAS, [('chr_i_input', 'chr_i')])
                 ] 
        

    workflow_metal_mungng.connect(connections)    
    res=workflow_metal_mungng.run('SGE',plugin_args=dict(dont_resubmit_completed_jobs= True, overwrite= True, qsub_args=f'-q {machineName} -o /ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/Exta_temp_files/Logs/ -j y'))

    return "Completed"
