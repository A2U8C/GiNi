import os
import sys
import click


import nipype
import nipype.pipeline.engine as pe
import nipype.interfaces.utility as niu
from nipype.interfaces.utility import Function


sys.path.append('/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/')

from METAL_Module.metal_script import *

from LDSC_Module.ldsc_manager import * ##General_Munge
#python gini_main.py input-module-wrapper --input_txt /ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/temp_file.txt --n_studies 2 --ethnicity European --analysis_list Heritability --tissues_cells Astrocytes --meta random


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


def Non_Metal_manager(study_files:list):
    
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


    workflow_metal_mungng = pe.Workflow(name='Non_METAL_Munging_wf2')
    workflow_metal_mungng.base_dir = "/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/Exta_temp_files/"


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
    node1_LDSC_Munge.inputs.kwargs=dict(N_col="TotalN",snp="MarkerName",frq="Freq1",signed_sumstats="EffectARE",out=extras_LDSC_Munge_files+"/",a1="Allele1",a2="Allele2",p="PvalueARE")
    
    
    node2_LDSC_Heritability = pe.Node(Function(input_names=['filePathInp_Munged'],
                    output_names=['LDSC_Heritability_out'],
                    function=HeritabilityLDSC),
            name='node2_LDSC_Heritability')
    

    connections=[(inputnode, node1_LDSC_Munge, [('filePathInp', 'filePathInp')]),
                    (node1_LDSC_Munge, node2_LDSC_Heritability, [('LDSC_Munge_out', 'filePathInp_Munged')])] 
    
    workflow_metal_mungng.connect(connections)    
    workflow_metal_mungng.run('SGE',plugin_args=dict(dont_resubmit_completed_jobs= True, overwrite= True, qsub_args=f'-q runnow.q -o /ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/Exta_temp_files/Logs/ -j y'))
    
    pass

@gini_gwas.command()
@click.option('--input_txt', help='input file. Each line contains the path location of the txt file for each study, which contains file path for each trait. Trait files should be named according to: study_name___trait_name... MAKE SURE THE STUDY NAME AND THE TRAIT NAME ARE SEPERATED BY 3 UNDERSCORES. THERE SHOULD BE ONLY ONE OCCURENCE OF 3 UNDERSCORES TOGETHER')
@click.option('--n_studies', type=int, help='Number of studies. To decide if Meta Analysis is needed or not')
@click.option('--ethnicity', metavar='| For Reference panel selection', help="""Ethnicity for selecting appropriate Reference panel. Select one of the available options: 'European' or 'African' or 'East_Asian' or 'South_Asian' or 'Middle_South_American' or 'custom_entry'""")
@click.option('--analysis_list', metavar='| Selection of Analyses and tests', help="""Specify one or more of the analyses which you want to be executed. Following are the available options: 'Global_Correlation' or 'Local_Correlation' or 'Heritability' or 'GSMR' or 'TWAS' or 'ALL'""")
@click.option('--tissues_cells',metavar=' | Selection of Brain Tissue and Cells',help="""Specify one or more of the Brain Tissue and Cells for which you want to performe TWAS and heritability. Following are the available options: 'Astrocytes' or 'Oligodendrocytes' or 'Neurons' or 'All_Cells' or 'Amygdala' or 'Anterior_cingulate' or 'Caudate' or 'Cerebellar' or 'Cerebellum' or 'Cortex' or 'Frontal_cortex' or 'Hippocampus' or 'Hypothalamus' or 'Nucleus_accumbens' or 'Putamen' or 'Spinal_cord' or 'Substantia_nigra' or 'All_Brain_tissue' or 'Everything'""")
@click.option('--meta', metavar='| Analysis Type', help="""Analysis type for Meta-Analysis; Please enter 'fixed' for the fixed-effects or 'random' for random-effects""")
def input_module_wrapper(input_txt,n_studies, ethnicity,analysis_list,tissues_cells,meta):

    analysis_list=analysis_list.split(',')
    tissues_cells=tissues_cells.split(',')
    
    meta_Analysis_op, study_files = input_check(input_txt,n_studies, ethnicity,analysis_list,tissues_cells,meta)


    #check if the system on which the program is running is a HPC: 
    # os.popen('hostname').read().encode('utf-8')
    #  psutil.virtual_memory()
    
    
    if meta_Analysis_op:
        metal_path_gz_files=metal_Analysis_Module(study_files,meta.upper())
    else:
        Non_Metal_manager(study_files)



if __name__ == "__main__":
    gini_gwas()
