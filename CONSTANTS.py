import os
import sys

import nipype
import nipype.pipeline.engine as pe
import nipype.interfaces.utility as niu
from nipype.interfaces.utility import Function


sys.path.append('/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/')
from METAL_Module.metal_script import *
from LDSC_Module.ldsc_manager import * ##General_Munge



cell_type_ldsc_ref="/ifs/loni/faculty/njahansh/nerds/ravi/genetics/ldsc/ldsc_seg_ldscores/"
#change it to the "-"; Mention  in the documentation
sep_btw_study_trait="___" #Seperator between study_name___trait_name.regenie 

os.chdir('/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/')
root_dir_to_Gini=os.getcwd()#+"/"
# root_dir_to_Gini="/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/"

Extra_temp_files_dict=dict( extras_folder_path = os.path.join(root_dir_to_Gini, "Exta_temp_files"),# Check if folder exists; to add temporary folder in it
extras_metal_files = os.path.join(root_dir_to_Gini, "Exta_temp_files/METAL_Out"),# Check if Metal extras folder exists; to add temporary files in it
extras_metal_qced_files = os.path.join(root_dir_to_Gini, "Exta_temp_files/METAL_Out/METAL_QCED_file"),# to add QCED files in it
extras_metal_input_files = os.path.join(root_dir_to_Gini, "Exta_temp_files/METAL_Out/METAL_input_files"),# to add .metal input files in it
extras_metal_output_files = os.path.join(root_dir_to_Gini, "Exta_temp_files/METAL_Out/METAL_output_files"),# For METAL results
extras_LDSC_Munge_files = os.path.join(root_dir_to_Gini, "Exta_temp_files/LDSC_Out/Munged_results"),# to add LDSC Mumnge results
extras_LDSC_Heri_files = os.path.join(root_dir_to_Gini, "Exta_temp_files/LDSC_Out/Heritability_Results"),# to add Heritability resutls
extras_LDSC_Heri__results = os.path.join(root_dir_to_Gini, "Exta_temp_files/LDSC_Out/Heritability_CSV_files"),# to add Heritability tables from logs
extras_LDSC_rG_results    =   os.path.join(root_dir_to_Gini, "Exta_temp_files/LDSC_Out/rG_Results"),
rG_CSV_files    =   os.path.join(root_dir_to_Gini, "Exta_temp_files/LDSC_Out/rG_CSV_files"),
extra_cell_h2=os.path.join(root_dir_to_Gini, "Exta_temp_files/LDSC_Out/Cell_Type_Heritability_Out")
                                      )

rG_folder="/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/Exta_temp_files/INPUT_TESTING_PACKAGE/Enigma_SE_GC/"
def file_name_process(fileName:str):
    return fileName.replace(".tbl","").replace(".csv","").replace(".txt","").replace(".gz","").replace(".tsv","").replace(".vcf","")

# node_METAL_Results = pe.Node(Function(input_names=['trait_name', 'trait_study_loc_paths_str', 'meta_type'],
#                       output_names=['out_zipped_metal'],
#                       function=metal_improved_execution_function),
#              name='node_METAL_Results')

# node_LDSC_Munge = pe.Node(Function(input_names=['filePathInp','kwargs'],
#                     output_names=['LDSC_Munge_out'],
#                     function=General_Munge),
#             name='node_LDSC_Munge')


# node_LDSC_Heritability = pe.Node(Function(input_names=['filePathInp_Munged'],
#                     output_names=['LDSC_Heritability_out'],
#                     function=HeritabilityLDSC),
#             name='node_LDSC_Heritability')

# Nipype_Nodes_dict=dict(node_METAL_Results=node_METAL_Results,
#                        node_LDSC_Munge=node_LDSC_Munge,
#                        node_LDSC_Heritability=node_LDSC_Heritability
#                        )