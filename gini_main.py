import os
import sys
import click


import nipype
import nipype.pipeline.engine as pe
import nipype.interfaces.utility as niu
from nipype.interfaces.utility import Function
from nipype import config
import shutil

sys.path.append('/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/')

from METAL_Module.metal_script import *
from CONSTANTS import *
import CONSTANTS

from LDSC_Module.ldsc_manager import * ##General_Munge
from LAVA_Module.LAVA_script import * ##Lava_FilePrep Lava_input_file
#python gini_main.py input-module-wrapper --input_txt /ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/temp_file.txt --n_studies 2 --ethnicity European --analysis_list Heritability --tissues_cells Astrocytes --meta random

#python gini_main.py input-module-wrapper --input_txt /ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/temp_file_nonMetal.txt --n_studies 1 --ethnicity European --analysis_list Heritability --tissues_cells Astrocytes --meta random


########################################################Review:
### Check if Click handle the inout 


##############################################################



@click.group()
def gini_gwas():
    pass


def file_checker(file_path_in)->bool:
    if not (file_path_in and os.path.exists(file_path_in)):
        # print(file_path_in)
        print(f"The file {file_path_in} does not exist or the path is wrong")
        sys.exit()
    return True


def input_check(in_txt1,n_studies, ethnicity,analysis_list,tissues_cells,meta_type)->tuple([bool,list]):

    analysis_global_set={ 'ALL','Global_Correlation', 'Local_Correlation', 'Heritability', 'GSMR','TWAS'}
    ethnicity_global_list=['European', 'African','East_Asian', 'South_Asian', 'Middle_South_American', 'custom_entry']
    tissues_cells_global_list={'All_Cells', 'Astrocytes', 'Oligodendrocytes', 'Neurons', 'Amygdala', 'Anterior_cingulate', 'Caudate', 'Cerebellar', 'Cerebellum', 'Cortex', 'Frontal_cortex', 'Hippocampus', 'Hypothalamus', 'Nucleus_accumbens', 'Putamen', 'Spinal_cord', 'Substantia_nigra', 'All_Brain_tissue', 'Everything'}
    meta_types_global=['random','fixed']

    ##### Check the text file and the files mentioned in the text file, if it exists###
    file_checker(in_txt1)
    
    with open(in_txt1, 'r') as file:
        n_files = file.read().splitlines()

    if len(n_files)==0:
        print('The input file is empty')
        sys.exit()

    if len(n_files)!=n_studies:
        print('The number of studies in the input file does not match with the number of studies mentioned')
        sys.exit()

    if len(n_files)>1:
        meta_Analysis_op=True
        if meta_type.lower() not in meta_types_global:
            print('Invalid Meta Analysis Format')
            print('Avaliable Meta Analysis Formats are: ', ", ".join(meta_types_global))
            sys.exit()
    else:
        meta_Analysis_op=False
    #############################################################################
    #### Changer each element to lower case to check if the current entry is present in the list
    if not (ethnicity and ethnicity in ethnicity_global_list):
        print(ethnicity not in ethnicity_global_list)
        print('Please enter a valid ethnicity')
        print('Aavliable ethnicity options: ',", ".join(ethnicity_global_list))
        sys.exit()
    
    if not (analysis_list and set(analysis_list).issubset(analysis_global_set)):
        print('Please enter valid analysis list. Following analysis from your list is/are not available: (', ", ".join(list(set(analysis_list).difference( analysis_global_set)))+')')
        print('Aavliable Analysis list: ',", ".join(['ALL']+[i for i in list(analysis_global_set) if i!='ALL']))
        sys.exit()

    if not (tissues_cells and set(tissues_cells).issubset(tissues_cells_global_list)):
        print('Please enter valid tiisues and cells list. Following tissues and cells from your list is/are not available: (', ", ".join(list(set(tissues_cells).difference(tissues_cells_global_list)))+')')
        print('Aavliable Tissues and cells list: ',", ".join(['All_Cells']+[i for i in list(tissues_cells_global_list) if i!='All_Cells']))
        sys.exit()

    return meta_Analysis_op, n_files


def Non_Metal_manager(in_txt_file:list):
    with open(in_txt_file[0], 'r') as file:
        study_files = file.read().splitlines()
        
    heritability_res_files=[]
    
    workflow_metal_mungng = pe.Workflow(name='Non_METAL_Munging_wf2')
    workflow_metal_mungng.base_dir = "/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/Exta_temp_files/"



    for study_i_inp_file in study_files:
        heritability_res_files.append(CONSTANTS.Extra_temp_files_dict["extras_LDSC_Heri_files"]+"/Heritability_"+study_i_inp_file.replace(".tbl","").replace(".csv","").replace(".txt","").replace(".gz","").replace(".tsv","").replace(".vcf","").split("/")[-1].split(".gz")[0]+".sumstats.log")

    iterable_new_list=study_files
    
    print(iterable_new_list)
    inputnode = pe.Node(
        niu.IdentityInterface(
            fields=['filePathInp']
            ),
        synchronize=True,
        iterables=[('filePathInp',iterable_new_list)],
        name='inputnode'
    )
    
    node1_LDSC_Munge = pe.Node(Function(input_names=['filePathInp','kwargs'],
                    output_names=['LDSC_Munge_out'],
                    function=General_Munge),
            name='node1_LDSC_Munge')
    node1_LDSC_Munge.inputs.kwargs=dict(N_col="TotalN",snp="MarkerName",frq="Freq1",signed_sumstats="EffectARE",out=CONSTANTS.Extra_temp_files_dict["extras_LDSC_Munge_files"]+"/",a1="Allele1",a2="Allele2",p="PvalueARE")
    
    
    node2_LDSC_Heritability = pe.Node(Function(input_names=['filePathInp_Munged'],
                    output_names=['LDSC_Heritability_out'],
                    function=HeritabilityLDSC),
            name='node2_LDSC_Heritability')
    
    
    node3_LDSC_rG_manager = pe.Node(Function(input_names=['Munged_file','hpc_bool','machineName'],
                    output_names=['LDSC_Heritability_out'],
                    function=rG_Node_Manager),
            name='node3_LDSC_rG_manager')
    
    node3_LDSC_rG_manager.inputs.hpc_bool=True
    node3_LDSC_rG_manager.inputs.machineName="runnow.q"

    connections=[(inputnode, node1_LDSC_Munge, [('filePathInp', 'filePathInp')]),
                    (node1_LDSC_Munge, node2_LDSC_Heritability, [('LDSC_Munge_out', 'filePathInp_Munged')]),
                    # (node1_LDSC_Munge, node3_LDSC_rG_manager, [('LDSC_Munge_out', 'Munged_file')])
                    ] 
    
    workflow_metal_mungng.connect(connections)    
    workflow_metal_mungng.run('SGE',plugin_args=dict(dont_resubmit_completed_jobs= True, overwrite= True, qsub_args=f'-q runnow.q -o /ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/Exta_temp_files/Logs/ -j y'))
    
    print("All Heritability results: /*-/*--*/-*/-*/-/-*//-*/*-/*/-/*-/*-*")
    print(heritability_res_files)
    Heritability_Log_Extraction(heritability_res_files)
    pass



def Non_Metal_manager_ver2(trait_study_dict:dict,gwas_format:str,lava_control_cases,hpc_bool:bool=False,machineName:str="iniadmin7.q"):
    study_files = sum(trait_study_dict.values(), [])
    
        
    heritability_res_files=[]

    all_ENIGMA__munged_cortical=[os.path.join(CONSTANTS.rG_folder["wSA"], f) for f in os.listdir(CONSTANTS.rG_folder["wSA"]) if (os.path.isfile(os.path.join(CONSTANTS.rG_folder["wSA"], f)) and f.split(".")[-1].lower()=="gz")]+[os.path.join(CONSTANTS.rG_folder["wTHICK"], f) for f in os.listdir(CONSTANTS.rG_folder["wTHICK"]) if (os.path.isfile(os.path.join(CONSTANTS.rG_folder["wTHICK"], f)) and f.split(".")[-1].lower()=="gz")]

    global_rG_res_files=[CONSTANTS.Extra_temp_files_dict["extras_LDSC_rG_results"]+"/"+ trait_files.split(",")[0].split("/")[-1].split(".")[0]+"___"+trait_files.split(",")[1].split("/")[-1].split(".")[0]+".log" for trait_files in CONSTANTS.trait_Combinations_for_rG(study_files,all_ENIGMA__munged_cortical)]


    





    for study_i_inp_file in study_files:
        heritability_res_files.append(CONSTANTS.Extra_temp_files_dict["extras_LDSC_Heri_files"]+"/Heritability_"+study_i_inp_file.replace(".tbl","").replace(".csv","").replace(".txt","").replace(".gz","").replace(".tsv","").replace(".vcf","").split("/")[-1].split(".gz")[0]+".sumstats.log")
        


    '''
    Change the way the relative path is being set using the config file
    '''

    if hpc_bool:
        config_dict={'execution': {'job_finished_timeout': '60'}}
        config.update_config(config_dict)


        wf_name='Non_METAL_Munging_wf3'
        workflow_metal_mungng = pe.Workflow(name=wf_name)
        workflow_metal_mungng.base_dir = "/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/Exta_temp_files/"


        if lava_control_cases==None or gwas_format=='metal':
            iterable_new_list=[[el,"NA","NA"] for el in study_files] #trait,cases,control
        
        inputnode = pe.Node(
            niu.IdentityInterface(
                fields=['filePathInp','Controls_N','Cases_N']
                ),
            synchronize=True,
            iterables=[('filePathInp','Controls_N','Cases_N'),iterable_new_list],
            # iterables=[('filePathInp',iterable_new_list)],
            name='inputnode'
        )
        
        node1_LDSC_Munge = pe.Node(Function(input_names=['filePathInp','kwargs'],
                        output_names=['LDSC_Munge_out'],
                        function=General_Munge),
                name='node1_LDSC_Munge')
        node1_LDSC_Munge.inputs.kwargs=dict(N_col="TotalN",snp="MarkerName",frq="Freq1",signed_sumstats="EffectARE",out=CONSTANTS.Extra_temp_files_dict["extras_LDSC_Munge_files"]+"/",a1="Allele1",a2="Allele2",p="PvalueARE")
        
        
        node2_LDSC_Heritability = pe.Node(Function(input_names=['filePathInp_Munged'],
                        output_names=['LDSC_Heritability_out'],
                        function=HeritabilityLDSC),
                name='node2_LDSC_Heritability')
        
        
        node3_LDSC_Cell_h2 = pe.Node(Function(input_names=['files_Munged','CellAnalysisName','LDCTSFile'],
                        output_names=['LDSC_Cell_H2_out'],
                        function=CellTypeLDSC),
                name='node3_LDSC_Cell_h2')
        node3_LDSC_Cell_h2.inputs.CellAnalysisName="Corces"
        node3_LDSC_Cell_h2.inputs.LDCTSFile="Corces_ATAC.ldcts"
        
        node4_LDSC_rG_manager = pe.Node(Function(input_names=['Munged_file','hpc_bool','machineName'],
                                                 output_names=['LDSC_trait_out'],
                                                 function=rG_Node_Manager),
                                                 name='node4_LDSC_rG_manager')
        node4_LDSC_rG_manager.inputs.hpc_bool=True
        node4_LDSC_rG_manager.inputs.machineName=machineName


        
        node_LAVA_FilePrep = pe.Node(Function(input_names=['filename','Gwas_format'],
                                                 output_names=['LAVA_prep_out'],
                                                 function=Lava_FilePrep),
                                                 name='node_LAVA_FilePrep')
        node_LAVA_FilePrep.inputs.Gwas_format=gwas_format

        

        
        node_LAVA_input_file = pe.Node(Function(input_names=['trait_path', 'Enigma_input','lava_control','lava_cases'],
                                                 output_names=['LAVA_input_out'],
                                                 function=Lava_input_file),
                                                 name='node_LAVA_input_file')
        node_LAVA_input_file.inputs.Enigma_input='cortical'


        node_LAVA_Matrix_formation = pe.Node(Function(input_names=['rG_log_files'],
                                                 output_names=['LAVA_matrix_file_out'],
                                                 function=LAVA_Matrix_Formation),
                                                 name='node_LAVA_Matrix_formation')
        

        
        node_LAVA_Shell_Creation = pe.Node(Function(input_names=['lava_matrix','lava_case_control'],
                                                 output_names=['LAVA_shell_file_out'],
                                                 function=LAVA_shell_call_script),
                                                 name='node_LAVA_Shell_Creation')
        

        
##################################################################################################
        # # iterable_new_list_rG_Fixer=["/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/Exta_temp_files/LDSC_Out/Munged_results/ukbb_white_abcd_white____Witelson5_Genu_Area1.sumstats.gz"]

        # # iterable_new_list_rG_Fixer=["/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/Exta_temp_files/LDSC_Out/rG_Results/ukbb_white_abcd_white____Witelson5_Genu_Area1___ENIGMA3_mixed_se_wSA_Mean_transversetemporal_surfavg_20190429.logabcdFileX123__/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/Exta_temp_files/LDSC_Out/rG_Results/ukbb_white_abcd_white____Witelson5_Genu_Area1___ENIGMA3_mixed_se_wSA_Mean_lingual_surfavg_20190429.logabcdFileX123__/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/Exta_temp_files/LDSC_Out/rG_Results/ukbb_white_abcd_white____Witelson5_Genu_Area1___ENIGMA3_mixed_se_wSA_Mean_postcentral_surfavg_20190429.logabcdFileX123__/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/Exta_temp_files/LDSC_Out/rG_Results/ukbb_white_abcd_white____Witelson5_Genu_Area1___ENIGMA3_mixed_se_wSA_Mean_superiorparietal_surfavg_20190429.logabcdFileX123__/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/Exta_temp_files/LDSC_Out/rG_Results/ukbb_white_abcd_white____Witelson5_Genu_Area1___ENIGMA3_mixed_se_wSA_Mean_parsorbitalis_surfavg_20190429.logabcdFileX123__/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/Exta_temp_files/LDSC_Out/rG_Results/ukbb_white_abcd_white____Witelson5_Genu_Area1___ENIGMA3_mixed_se_wSA_Mean_superiorfrontal_surfavg_20190429.logabcdFileX123__/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/Exta_temp_files/LDSC_Out/rG_Results/ukbb_white_abcd_white____Witelson5_Genu_Area1___ENIGMA3_mixed_se_wSA_Mean_precentral_surfavg_20190429.logabcdFileX123__/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/Exta_temp_files/LDSC_Out/rG_Results/ukbb_white_abcd_white____Witelson5_Genu_Area1___ENIGMA3_mixed_se_wSA_Mean_caudalmiddlefrontal_surfavg_20190429.logabcdFileX123__/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/Exta_temp_files/LDSC_Out/rG_Results/ukbb_white_abcd_white____Witelson5_Genu_Area1___ENIGMA3_mixed_se_wSA_Mean_middletemporal_surfavg_20190429.logabcdFileX123__/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/Exta_temp_files/LDSC_Out/rG_Results/ukbb_white_abcd_white____Witelson5_Genu_Area1___ENIGMA3_mixed_se_wSA_Mean_parsopercularis_surfavg_20190429.logabcdFileX123__/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/Exta_temp_files/LDSC_Out/rG_Results/ukbb_white_abcd_white____Witelson5_Genu_Area1___ENIGMA3_mixed_se_wSA_Mean_parahippocampal_surfavg_20190429.logabcdFileX123__/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/Exta_temp_files/LDSC_Out/rG_Results/ukbb_white_abcd_white____Witelson5_Genu_Area1___ENIGMA3_mixed_se_wSA_Mean_cuneus_surfavg_20190429.logabcdFileX123__/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/Exta_temp_files/LDSC_Out/rG_Results/ukbb_white_abcd_white____Witelson5_Genu_Area1___ENIGMA3_mixed_se_wSA_Mean_parstriangularis_surfavg_20190429.logabcdFileX123__/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/Exta_temp_files/LDSC_Out/rG_Results/ukbb_white_abcd_white____Witelson5_Genu_Area1___ENIGMA3_mixed_se_wSA_Mean_inferiortemporal_surfavg_20190429.logabcdFileX123__/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/Exta_temp_files/LDSC_Out/rG_Results/ukbb_white_abcd_white____Witelson5_Genu_Area1___ENIGMA3_mixed_se_wSA_Mean_superiortemporal_surfavg_20190429.logabcdFileX123__/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/Exta_temp_files/LDSC_Out/rG_Results/ukbb_white_abcd_white____Witelson5_Genu_Area1___ENIGMA3_mixed_se_wSA_Mean_lateralorbitofrontal_surfavg_20190429.logabcdFileX123__/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/Exta_temp_files/LDSC_Out/rG_Results/ukbb_white_abcd_white____Witelson5_Genu_Area1___ENIGMA3_mixed_se_wSA_Mean_isthmuscingulate_surfavg_20190429.logabcdFileX123__/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/Exta_temp_files/LDSC_Out/rG_Results/ukbb_white_abcd_white____Witelson5_Genu_Area1___ENIGMA3_mixed_se_wSA_Mean_lateraloccipital_surfavg_20190429.logabcdFileX123__/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/Exta_temp_files/LDSC_Out/rG_Results/ukbb_white_abcd_white____Witelson5_Genu_Area1___ENIGMA3_mixed_se_wSA_Mean_insula_surfavg_20190429.logabcdFileX123__/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/Exta_temp_files/LDSC_Out/rG_Results/ukbb_white_abcd_white____Witelson5_Genu_Area1___ENIGMA3_mixed_se_wSA_Mean_paracentral_surfavg_20190429.logabcdFileX123__/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/Exta_temp_files/LDSC_Out/rG_Results/ukbb_white_abcd_white____Witelson5_Genu_Area1___ENIGMA3_mixed_se_wSA_Mean_supramarginal_surfavg_20190429.logabcdFileX123__/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/Exta_temp_files/LDSC_Out/rG_Results/ukbb_white_abcd_white____Witelson5_Genu_Area1___ENIGMA3_mixed_se_wSA_Mean_rostralanteriorcingulate_surfavg_20190429.logabcdFileX123__/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/Exta_temp_files/LDSC_Out/rG_Results/ukbb_white_abcd_white____Witelson5_Genu_Area1___ENIGMA3_mixed_se_wSA_Mean_temporalpole_surfavg_20190429.logabcdFileX123__/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/Exta_temp_files/LDSC_Out/rG_Results/ukbb_white_abcd_white____Witelson5_Genu_Area1___ENIGMA3_mixed_se_wSA_Mean_caudalanteriorcingulate_surfavg_20190429.logabcdFileX123__/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/Exta_temp_files/LDSC_Out/rG_Results/ukbb_white_abcd_white____Witelson5_Genu_Area1___ENIGMA3_mixed_se_wSA_Mean_posteriorcingulate_surfavg_20190429.logabcdFileX123__/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/Exta_temp_files/LDSC_Out/rG_Results/ukbb_white_abcd_white____Witelson5_Genu_Area1___ENIGMA3_mixed_se_wSA_Mean_inferiorparietal_surfavg_20190429.logabcdFileX123__/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/Exta_temp_files/LDSC_Out/rG_Results/ukbb_white_abcd_white____Witelson5_Genu_Area1___ENIGMA3_mixed_se_wSA_Mean_frontalpole_surfavg_20190429.logabcdFileX123__/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/Exta_temp_files/LDSC_Out/rG_Results/ukbb_white_abcd_white____Witelson5_Genu_Area1___ENIGMA3_mixed_se_wSA_Mean_medialorbitofrontal_surfavg_20190429.logabcdFileX123__/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/Exta_temp_files/LDSC_Out/rG_Results/ukbb_white_abcd_white____Witelson5_Genu_Area1___ENIGMA3_mixed_se_wSA_Mean_rostralmiddlefrontal_surfavg_20190429.logabcdFileX123__/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/Exta_temp_files/LDSC_Out/rG_Results/ukbb_white_abcd_white____Witelson5_Genu_Area1___ENIGMA3_mixed_se_wSA_Mean_entorhinal_surfavg_20190429.logabcdFileX123__/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/Exta_temp_files/LDSC_Out/rG_Results/ukbb_white_abcd_white____Witelson5_Genu_Area1___ENIGMA3_mixed_se_wSA_Mean_precuneus_surfavg_20190429.logabcdFileX123__/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/Exta_temp_files/LDSC_Out/rG_Results/ukbb_white_abcd_white____Witelson5_Genu_Area1___ENIGMA3_mixed_se_wSA_Mean_fusiform_surfavg_20190429.logabcdFileX123__/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/Exta_temp_files/LDSC_Out/rG_Results/ukbb_white_abcd_white____Witelson5_Genu_Area1___ENIGMA3_mixed_se_wSA_Mean_bankssts_surfavg_20190429.logabcdFileX123__/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/Exta_temp_files/LDSC_Out/rG_Results/ukbb_white_abcd_white____Witelson5_Genu_Area1___ENIGMA3_mixed_se_wSA_Mean_pericalcarine_surfavg_20190429.logabcdFileX123__/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/Exta_temp_files/LDSC_Out/rG_Results/ukbb_white_abcd_white____Witelson5_Genu_Area1___ENIGMA3_mixed_se_wTHICK_Mean_rostralmiddlefrontal_thickavg_20200522.logabcdFileX123__/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/Exta_temp_files/LDSC_Out/rG_Results/ukbb_white_abcd_white____Witelson5_Genu_Area1___ENIGMA3_mixed_se_wTHICK_Mean_inferiorparietal_thickavg_20200522.logabcdFileX123__/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/Exta_temp_files/LDSC_Out/rG_Results/ukbb_white_abcd_white____Witelson5_Genu_Area1___ENIGMA3_mixed_se_wTHICK_Mean_cuneus_thickavg_20200522.logabcdFileX123__/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/Exta_temp_files/LDSC_Out/rG_Results/ukbb_white_abcd_white____Witelson5_Genu_Area1___ENIGMA3_mixed_se_wTHICK_Mean_caudalanteriorcingulate_thickavg_20200522.logabcdFileX123__/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/Exta_temp_files/LDSC_Out/rG_Results/ukbb_white_abcd_white____Witelson5_Genu_Area1___ENIGMA3_mixed_se_wTHICK_Mean_bankssts_thickavg_20200522.logabcdFileX123__/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/Exta_temp_files/LDSC_Out/rG_Results/ukbb_white_abcd_white____Witelson5_Genu_Area1___ENIGMA3_mixed_se_wTHICK_Mean_rostralanteriorcingulate_thickavg_20200522.logabcdFileX123__/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/Exta_temp_files/LDSC_Out/rG_Results/ukbb_white_abcd_white____Witelson5_Genu_Area1___ENIGMA3_mixed_se_wTHICK_Mean_entorhinal_thickavg_20200522.logabcdFileX123__/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/Exta_temp_files/LDSC_Out/rG_Results/ukbb_white_abcd_white____Witelson5_Genu_Area1___ENIGMA3_mixed_se_wTHICK_Mean_caudalmiddlefrontal_thickavg_20200522.logabcdFileX123__/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/Exta_temp_files/LDSC_Out/rG_Results/ukbb_white_abcd_white____Witelson5_Genu_Area1___ENIGMA3_mixed_se_wTHICK_Mean_parahippocampal_thickavg_20200522.logabcdFileX123__/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/Exta_temp_files/LDSC_Out/rG_Results/ukbb_white_abcd_white____Witelson5_Genu_Area1___ENIGMA3_mixed_se_wTHICK_Mean_precuneus_thickavg_20200522.logabcdFileX123__/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/Exta_temp_files/LDSC_Out/rG_Results/ukbb_white_abcd_white____Witelson5_Genu_Area1___ENIGMA3_mixed_se_wTHICK_Mean_supramarginal_thickavg_20200522.logabcdFileX123__/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/Exta_temp_files/LDSC_Out/rG_Results/ukbb_white_abcd_white____Witelson5_Genu_Area1___ENIGMA3_mixed_se_wTHICK_Mean_lateralorbitofrontal_thickavg_20200522.logabcdFileX123__/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/Exta_temp_files/LDSC_Out/rG_Results/ukbb_white_abcd_white____Witelson5_Genu_Area1___ENIGMA3_mixed_se_wTHICK_Mean_frontalpole_thickavg_20200522.logabcdFileX123__/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/Exta_temp_files/LDSC_Out/rG_Results/ukbb_white_abcd_white____Witelson5_Genu_Area1___ENIGMA3_mixed_se_wTHICK_Mean_posteriorcingulate_thickavg_20200522.logabcdFileX123__/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/Exta_temp_files/LDSC_Out/rG_Results/ukbb_white_abcd_white____Witelson5_Genu_Area1___ENIGMA3_mixed_se_wTHICK_Mean_lingual_thickavg_20200522.logabcdFileX123__/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/Exta_temp_files/LDSC_Out/rG_Results/ukbb_white_abcd_white____Witelson5_Genu_Area1___ENIGMA3_mixed_se_wTHICK_Mean_transversetemporal_thickavg_20200522.logabcdFileX123__/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/Exta_temp_files/LDSC_Out/rG_Results/ukbb_white_abcd_white____Witelson5_Genu_Area1___ENIGMA3_mixed_se_wTHICK_Mean_middletemporal_thickavg_20200522.logabcdFileX123__/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/Exta_temp_files/LDSC_Out/rG_Results/ukbb_white_abcd_white____Witelson5_Genu_Area1___ENIGMA3_mixed_se_wTHICK_Mean_superiorparietal_thickavg_20200522.logabcdFileX123__/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/Exta_temp_files/LDSC_Out/rG_Results/ukbb_white_abcd_white____Witelson5_Genu_Area1___ENIGMA3_mixed_se_wTHICK_Mean_isthmuscingulate_thickavg_20200522.logabcdFileX123__/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/Exta_temp_files/LDSC_Out/rG_Results/ukbb_white_abcd_white____Witelson5_Genu_Area1___ENIGMA3_mixed_se_wTHICK_Mean_paracentral_thickavg_20200522.logabcdFileX123__/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/Exta_temp_files/LDSC_Out/rG_Results/ukbb_white_abcd_white____Witelson5_Genu_Area1___ENIGMA3_mixed_se_wTHICK_Mean_fusiform_thickavg_20200522.logabcdFileX123__/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/Exta_temp_files/LDSC_Out/rG_Results/ukbb_white_abcd_white____Witelson5_Genu_Area1___ENIGMA3_mixed_se_wTHICK_Mean_insula_thickavg_20200522.logabcdFileX123__/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/Exta_temp_files/LDSC_Out/rG_Results/ukbb_white_abcd_white____Witelson5_Genu_Area1___ENIGMA3_mixed_se_wTHICK_Mean_precentral_thickavg_20200522.logabcdFileX123__/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/Exta_temp_files/LDSC_Out/rG_Results/ukbb_white_abcd_white____Witelson5_Genu_Area1___ENIGMA3_mixed_se_wTHICK_Mean_medialorbitofrontal_thickavg_20200522.logabcdFileX123__/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/Exta_temp_files/LDSC_Out/rG_Results/ukbb_white_abcd_white____Witelson5_Genu_Area1___ENIGMA3_mixed_se_wTHICK_Mean_parsorbitalis_thickavg_20200522.logabcdFileX123__/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/Exta_temp_files/LDSC_Out/rG_Results/ukbb_white_abcd_white____Witelson5_Genu_Area1___ENIGMA3_mixed_se_wTHICK_Mean_superiortemporal_thickavg_20200522.logabcdFileX123__/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/Exta_temp_files/LDSC_Out/rG_Results/ukbb_white_abcd_white____Witelson5_Genu_Area1___ENIGMA3_mixed_se_wTHICK_Mean_superiorfrontal_thickavg_20200522.logabcdFileX123__/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/Exta_temp_files/LDSC_Out/rG_Results/ukbb_white_abcd_white____Witelson5_Genu_Area1___ENIGMA3_mixed_se_wTHICK_Mean_temporalpole_thickavg_20200522.logabcdFileX123__/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/Exta_temp_files/LDSC_Out/rG_Results/ukbb_white_abcd_white____Witelson5_Genu_Area1___ENIGMA3_mixed_se_wTHICK_Mean_parstriangularis_thickavg_20200522.logabcdFileX123__/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/Exta_temp_files/LDSC_Out/rG_Results/ukbb_white_abcd_white____Witelson5_Genu_Area1___ENIGMA3_mixed_se_wTHICK_Mean_parsopercularis_thickavg_20200522.logabcdFileX123__/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/Exta_temp_files/LDSC_Out/rG_Results/ukbb_white_abcd_white____Witelson5_Genu_Area1___ENIGMA3_mixed_se_wTHICK_Mean_inferiortemporal_thickavg_20200522.logabcdFileX123__/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/Exta_temp_files/LDSC_Out/rG_Results/ukbb_white_abcd_white____Witelson5_Genu_Area1___ENIGMA3_mixed_se_wTHICK_Mean_postcentral_thickavg_20200522.logabcdFileX123__/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/Exta_temp_files/LDSC_Out/rG_Results/ukbb_white_abcd_white____Witelson5_Genu_Area1___ENIGMA3_mixed_se_wTHICK_Mean_lateraloccipital_thickavg_20200522.logabcdFileX123__/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/Exta_temp_files/LDSC_Out/rG_Results/ukbb_white_abcd_white____Witelson5_Genu_Area1___ENIGMA3_mixed_se_wTHICK_Mean_pericalcarine_thickavg_20200522.log"]



        # inputnode_rG_Fixer=pe.Node(
        #     niu.IdentityInterface(
        #         fields=['filePathInp']
        #         ),
        #     synchronize=True,
        #     iterables=[('filePathInp',iterable_new_list_rG_Fixer)],
        #     # iterables=[('filePathInp',iterable_new_list)],
        #     name='inputnode'
        # )

        # connections=[
        #             # (inputnode_rG_Fixer, node4_LDSC_rG_manager, [('filePathInp', 'Munged_file')]),
        #             # (node4_LDSC_rG_manager, node_LAVA_Matrix_formation, [('LDSC_trait_out', 'rG_log_files')]),
                    
        #             (inputnode_rG_Fixer, node_LAVA_Matrix_formation, [('filePathInp', 'rG_log_files')])
        #                 ] 
        
##################################################################################################
        connections=[
                    (inputnode, node_LAVA_input_file, [('filePathInp', 'trait_path')]),
                    (inputnode, node_LAVA_input_file, [('Controls_N', 'lava_control')]),
                    (inputnode, node_LAVA_input_file, [('Cases_N', 'lava_cases')]),
                    (inputnode, node1_LDSC_Munge, [('filePathInp', 'filePathInp')]),
                    (inputnode, node_LAVA_FilePrep, [('filePathInp', 'filename')]),
                    (node1_LDSC_Munge, node4_LDSC_rG_manager, [('LDSC_Munge_out', 'Munged_file')]),
                    (node4_LDSC_rG_manager, node_LAVA_Matrix_formation, [('LDSC_trait_out', 'rG_log_files')]),
                    (node_LAVA_input_file, node_LAVA_Shell_Creation, [('LAVA_input_out', 'lava_case_control')]),
                    (node_LAVA_Matrix_formation, node_LAVA_Shell_Creation, [('LAVA_matrix_file_out', 'lava_matrix')])
                        ] 

        workflow_metal_mungng.connect(connections)    
        workflow_metal_mungng.run('SGE',plugin_args=dict(dont_resubmit_completed_jobs= True, overwrite= True, qsub_args=f'-q {machineName} -o /ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/Exta_temp_files/Logs/ -j y'))
        # workflow_metal_mungng.run('SGE',plugin_args=dict(qsub_args=f'-q {machineName} -o /ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/Exta_temp_files/Logs/ -j y'))

        # workflow_metal_mungng.run()

        # shutil.rmtree(workflow_metal_mungng.base_dir+wf_name+"/") # To remove the NiPype LOGS
        #  (node4_LDSC_rG_manager, node_LAVA_Matrix_formation, [('LDSC_trait_out', 'rG_log_files')])
        #  (node1_LDSC_Munge, node2_LDSC_Heritability, [('LDSC_Munge_out', 'filePathInp_Munged')]),
        # (node1_LDSC_Munge, node3_LDSC_Cell_h2, [('LDSC_Munge_out', 'files_Munged')]),
        

    
    else:
        for study_file in study_files:
            HeritabilityLDSC(General_Munge(study_file,dict(N_col="TotalN",snp="MarkerName",frq="Freq1",signed_sumstats="EffectARE",out=CONSTANTS.Extra_temp_files_dict["extras_LDSC_Munge_files"]+"/",a1="Allele1",a2="Allele2",p="PvalueARE")
        ))
    
    # print(f"Heritability results are extracted saved in a CSV file in {CONSTANTS.Extra_temp_files_dict['extras_LDSC_Heri__results']}")
    # Heritability_Log_Extraction(heritability_res_files)
    rG_Log_Extraction(global_rG_res_files)
    # for el in global_rG_res_files:
    #     print(el)


def rG_Node_Manager(Munged_file:str,hpc_bool:bool=False,machineName:str="runnow.q"):

    import os
    import sys
    import click

    sys.path.append('/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/')
    
    from LDSC_Module import ldsc_manager ##General_Munge
    import CONSTANTS 

    import nipype
    import nipype.pipeline.engine as pe
    import nipype.interfaces.utility as niu
    from nipype.interfaces.utility import Function
    import shutil
        

    # All combinations with the Enigma GWAS
    iterable_new_list=[Munged_file+","+os.path.join(CONSTANTS.rG_folder["wSA"], f) for f in os.listdir(CONSTANTS.rG_folder["wSA"]) if (os.path.isfile(os.path.join(CONSTANTS.rG_folder["wSA"], f)) and f.split(".")[-1].lower()=="gz")]+[Munged_file+","+os.path.join(CONSTANTS.rG_folder["wTHICK"], f) for f in os.listdir(CONSTANTS.rG_folder["wTHICK"]) if (os.path.isfile(os.path.join(CONSTANTS.rG_folder["wTHICK"], f)) and f.split(".")[-1].lower()=="gz")]


    wf_name='Non_METAL_rG_wf_'+Munged_file.split("/")[-1].split(".")[0]
    workflow_metal_mungng = pe.Workflow(name=wf_name)
    workflow_metal_mungng.base_dir = "/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/Exta_temp_files/"
    workflow_metal_mungng.config["job_finished_timeout"]=30

    inputnode = pe.Node(
            niu.IdentityInterface(
                fields=['trait_files']
                ),
            synchronize=True,
            iterables=[('trait_files',iterable_new_list)],
            name='inputnode'
        )
        
    node1_LDSC_rG = pe.Node(Function(input_names=['trait_files'],
                    output_names=['out_file_comp'],
                    function=ldsc_manager.rG_LDSC),
            name='node1_LDSC_rG')
    
    
    connections=[(inputnode, node1_LDSC_rG, [('trait_files', 'trait_files')])
                 ] 
        

    workflow_metal_mungng.connect(connections)    
    res=workflow_metal_mungng.run('SGE',plugin_args=dict(dont_resubmit_completed_jobs= True, overwrite= True, qsub_args=f'-q {machineName} -o /ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/Exta_temp_files/Logs/ -j y'))


    print(res.nodes())
    outputs = CONSTANTS.file_joiner_str.join([x.result.outputs.out_file_comp for x in res.nodes() if x.name == 'node1_LDSC_rG'])
    print(outputs)
    # shutil.rmtree(workflow_metal_mungng.base_dir+wf_name+"/") # To remove the NiPype LOGS
    return outputs





def input_check_ver_2(in_txt1,n_studies, ethnicity,analysis_list,tissues_cells,meta_type)->tuple([bool,list]):
    analysis_global_set={ 'ALL','Global_Correlation', 'Local_Correlation', 'Heritability', 'GSMR','TWAS'}
    ethnicity_global_list=['European', 'African','East_Asian', 'South_Asian', 'Middle_South_American', 'custom_entry']
    tissues_cells_global_list={'All_Cells', 'Astrocytes', 'Oligodendrocytes', 'Neurons', 'Amygdala', 'Anterior_cingulate', 'Caudate', 'Cerebellar', 'Cerebellum', 'Cortex', 'Frontal_cortex', 'Hippocampus', 'Hypothalamus', 'Nucleus_accumbens', 'Putamen', 'Spinal_cord', 'Substantia_nigra', 'All_Brain_tissue', 'Everything'}
    meta_types_global=['random','fixed']

    ##### Check the text file and the files mentioned in the text file, if it exists###
    file_checker(in_txt1)
    
    with open(in_txt1, 'r') as file:
        n_files = file.read().splitlines()

    if len(n_files)==0:
        print('The input file is empty')
        sys.exit()

    #############################################################################
    #### Changer each element to lower case to check if the current entry is present in the list
    if not (ethnicity and ethnicity in ethnicity_global_list):
        print(ethnicity not in ethnicity_global_list)
        print('Please enter a valid ethnicity')
        print('Aavliable ethnicity options: ',", ".join(ethnicity_global_list))
        sys.exit()
    
    if not (analysis_list and set(analysis_list).issubset(analysis_global_set)):
        print('Please enter valid analysis list. Following analysis from your list is/are not available: (', ", ".join(list(set(analysis_list).difference( analysis_global_set)))+')')
        print('Aavliable Analysis list: ',", ".join(['ALL']+[i for i in list(analysis_global_set) if i!='ALL']))
        sys.exit()

    if not (tissues_cells and set(tissues_cells).issubset(tissues_cells_global_list)):
        print('Please enter valid tiisues and cells list. Following tissues and cells from your list is/are not available: (', ", ".join(list(set(tissues_cells).difference(tissues_cells_global_list)))+')')
        print('Aavliable Tissues and cells list: ',", ".join(['All_Cells']+[i for i in list(tissues_cells_global_list) if i!='All_Cells']))
        sys.exit()

    set_all_traits=set() # reference traits list; To cross check that trait_list are same for all the studies
    
    trait_study_dict={}#Create a dict for all traits, trait_i->[list of traits_paths for all studies ]
    trait_study_name={}#Create a dict for all traits, trait_i->{name of studies}

    

    for line_i in n_files:
        line_i_proc=line_i.strip("\n").strip("\r")
        file_checker(line_i_proc)# To check if the file path exists
        study_i_trait_j=line_i_proc.split("/")[-1].split(".")[0].split(sep_btw_study_trait)
        if len(study_i_trait_j)!=2:
            print('The following file name is not in the mentioned format: study_name___trait_name.file_type for file:', line_i_proc)
            sys.exit()
        study_i=study_i_trait_j[0]
        trait_j=study_i_trait_j[1]

        trait_study_name[trait_j]=trait_study_name.get(trait_j,set())
        trait_study_name[trait_j].add(study_i)

        trait_study_dict[trait_j]=trait_study_dict.get(trait_j,[])+[line_i_proc]
    
    tester_no_studies_ref=len(trait_study_dict[trait_j]) #To check that all the traits has same number of studies

    tester_study_name_set_ref=trait_study_name[trait_j] #To check that all the traits has same studies


    if tester_no_studies_ref!=n_studies:
        print('The number of studies in the input file does not match with the number of studies mentioned')
        sys.exit()


    for trait_el in trait_study_name.keys():
        trait_study_el_paths=trait_study_dict[trait_el]
        if len(trait_study_el_paths)!=n_studies:
            print('The number of studies for the trait: ',trait_el,' does not match with the number of studies mentioned')
            sys.exit()

        trait_study_el_set=trait_study_name[trait_el]

        if len(trait_study_el_set.symmetric_difference(tester_study_name_set_ref))!=0:
            print('The studies for the trait: ',trait_el,' does not match with other traits')
            sys.exit()
        

    if tester_no_studies_ref>1:
        meta_Analysis_op=True
        if meta_type.lower() not in meta_types_global:
            print('Invalid Meta Analysis Format')
            print('Avaliable Meta Analysis Formats are: ', ", ".join(meta_types_global))
            sys.exit()
    else:
        meta_Analysis_op=False

    ################################################ Temporary Folders Section ################################################################
    for tp_fold_i in CONSTANTS.Extra_temp_files_dict.keys():
        if not (os.path.exists(CONSTANTS.Extra_temp_files_dict[tp_fold_i]) and os.path.isdir(CONSTANTS.Extra_temp_files_dict[tp_fold_i])): # Check if folder exists; to add temporary folder in it
            os.makedirs(CONSTANTS.Extra_temp_files_dict[tp_fold_i])
    ################################################ Temporary Folders Section ################################################################
        
    return meta_Analysis_op, trait_study_dict



def METAL_NonMETAL_Manager(Analysis_ops:list, trait_study_dict:list):
    pass


def pain_manager(machineName="runnow.q"):
    trait_file_1="/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/TEMP_inputs/Pain_traits_txt1.txt"
    trait_file_2="/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/TEMP_inputs/Pain_traits_txt2.txt"

    workflow_metal_mungng = pe.Workflow(name='Pain_Munging_wf2')
    workflow_metal_mungng.base_dir = "/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/Exta_temp_files/"

    with open(trait_file_1, 'r') as file1:
        trait_file_1_files = file1.read().splitlines()

    with open(trait_file_2, 'r') as file2:
        trait_file_2_files = file2.read().splitlines()

    iterable_new_list=trait_file_1_files+trait_file_2_files
    Munging1_res_files=[]
    for study_i_inp_file in trait_file_1_files:
        Munging1_res_files.append(CONSTANTS.Extra_temp_files_dict["extras_LDSC_Munge_files"]+"/"+CONSTANTS.file_name_process(study_i_inp_file.split("/")[-1])+".sumstats.gz")

    Munging2_res_files=[]
    for study_i_inp_file in trait_file_2_files:
        Munging2_res_files.append(CONSTANTS.Extra_temp_files_dict["extras_LDSC_Munge_files"]+"/"+CONSTANTS.file_name_process(study_i_inp_file.split("/")[-1])+".sumstats.gz")

    all_combinations_for_rG=trait_Combinations_for_rG(Munging1_res_files,Munging2_res_files).to_list()

    print(all_combinations_for_rG)

    print(all_combinations_for_rG)
    inputnode = pe.Node(
        niu.IdentityInterface(
            fields=['filePathInp']
            ),
        synchronize=True,
        iterables=[('filePathInp',iterable_new_list)],
        name='inputnode'
    )
    
    node1_LDSC_Munge = pe.Node(Function(input_names=['filePathInp','kwargs'],
                    output_names=['LDSC_Munge_out'],
                    function=General_Munge),
            name='node1_LDSC_Munge')
    node1_LDSC_Munge.inputs.kwargs=dict(N_col="N",snp="SNP",frq="freq",OR=False,signed_sumstats="b",out=CONSTANTS.Extra_temp_files_dict["extras_LDSC_Munge_files"]+"/",a1="A1",a2="A2",p="p")
    
    connections=[(inputnode, node1_LDSC_Munge, [('filePathInp', 'filePathInp')])] 
    
    workflow_metal_mungng.connect(connections)    
    workflow_metal_mungng.run('SGE',plugin_args=dict(dont_resubmit_completed_jobs= True, overwrite= True, qsub_args=f'-q {machineName} -o /ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/Exta_temp_files/Logs/ -j y'))





    workflow_metal_rG = pe.Workflow(name='Pain_Munging_rG_wf2')
    workflow_metal_rG.base_dir = "/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/Exta_temp_files/"

    inputnode1 = pe.Node(
        niu.IdentityInterface(
            fields=['filePathInp']
            ),
        synchronize=True,
        iterables=[('filePathInp',all_combinations_for_rG)],
        name='inputnode1'
    )
    
    node3_LDSC_rG = pe.Node(Function(input_names=['trait_files'],
                    output_names=['LDSC_rG_out'],
                    function=rG_LDSC),
            name='node3_LDSC_rG')
    
    connections2=[(inputnode1, node3_LDSC_rG, [('filePathInp', 'trait_files')])] 
        
    
    workflow_metal_rG.connect(connections2)    
    workflow_metal_rG.run('SGE',plugin_args=dict(dont_resubmit_completed_jobs= True, overwrite= True, qsub_args=f'-q {machineName} -o /ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/Exta_temp_files/Logs/ -j y'))

    
    rg_Extractor_list=[]
    for trait_file_i in all_combinations_for_rG:
        trait_file_i=trait_file_i.replace("\n","").replace("\r","")
        trait_file_1=trait_file_i.split(",")[0] #
        trait_file_2=trait_file_i.split(",")[1] #
        
        rg_Extractor_list.append("/ifs/loni/faculty/njahansh/nerds/ankush/for_Ravi/Pain_Genetics/Pain_rG_NonDuplicates/"+trait_file_1.split("/")[-1].split(".")[0]+"___"+trait_file_2.split("/")[-1].split(".")[0]+".log")

    rG_Log_Extraction(rg_Extractor_list)


    pass






@gini_gwas.command()
@click.option('--input_txt', help='input file. Each line contains the path location of the txt file for each study, which contains file path for each trait. Trait files should be named according to: study_name___trait_name... MAKE SURE THE STUDY NAME AND THE TRAIT NAME ARE SEPERATED BY 3 UNDERSCORES. THERE SHOULD BE ONLY ONE OCCURENCE OF 3 UNDERSCORES TOGETHER')
@click.option('--n_studies', type=int, help='Number of studies. To decide if Meta Analysis is needed or not')
@click.option('--ethnicity', metavar='| For Reference panel selection', help="""Ethnicity for selecting appropriate Reference panel. Select one of the available options: 'European' or 'African' or 'East_Asian' or 'South_Asian' or 'Middle_South_American' or 'custom_entry'""")
@click.option('--analysis_list', metavar='| Selection of Analyses and tests', help="""Specify one or more of the analyses which you want to be executed. Following are the available options: 'Global_Correlation' or 'Local_Correlation' or 'Heritability' or 'GSMR' or 'TWAS' or 'ALL'""")
@click.option('--tissues_cells',metavar=' | Selection of Brain Tissue and Cells',help="""Specify one or more of the Brain Tissue and Cells for which you want to performe TWAS and heritability. Following are the available options: 'Astrocytes' or 'Oligodendrocytes' or 'Neurons' or 'All_Cells' or 'Amygdala' or 'Anterior_cingulate' or 'Caudate' or 'Cerebellar' or 'Cerebellum' or 'Cortex' or 'Frontal_cortex' or 'Hippocampus' or 'Hypothalamus' or 'Nucleus_accumbens' or 'Putamen' or 'Spinal_cord' or 'Substantia_nigra' or 'All_Brain_tissue' or 'Everything'""")
@click.option('--meta', metavar='| Analysis Type', help="""Analysis type for Meta-Analysis; Please enter 'fixed' for the fixed-effects or 'random' for random-effects""")
@click.option('--lava_control_cases', metavar='| # of Control/cases', help="""Path to the text file containing the count of control/cases for each trait. Format with the column heading: 'file_location' \t '#Controls' \t '#Cases' """, default=None)
@click.option('--gwas_format', metavar='| GWAS format', help="""Format of the GWASs which are given as input""", default='Regenie')
def input_module_wrapper(input_txt,n_studies, ethnicity,analysis_list,tissues_cells,meta,lava_control_cases,gwas_format):

    analysis_list=analysis_list.split(',')
    tissues_cells=tissues_cells.split(',')
    
    # meta_Analysis_op, study_files = input_check(input_txt,n_studies, ethnicity,analysis_list,tissues_cells,meta)
    meta_Analysis_op, study_files = input_check_ver_2(input_txt,n_studies, ethnicity,analysis_list,tissues_cells,meta)


    # #check if the system on which the program is running is a HPC: 
    # # os.popen('hostname').read().encode('utf-8')
    # #  psutil.virtual_memory()
    
    print(meta_Analysis_op)
    if meta_Analysis_op:
        metal_path_gz_files=metal_Analysis_Module(study_files,meta.upper())
    else:
        Non_Metal_manager_ver2(study_files,gwas_format=gwas_format,lava_control_cases=lava_control_cases,hpc_bool=True)

    # rG_Node_Manager("")
    # pain_manager()


if __name__ == "__main__":
    gini_gwas()
