import os
import sys
import json
'''
Update the paths for the downloaded (sQTLS, eQTLS, ENIGMA GWAS):    66GB
OPTIONS:
    Dropbox (Ask to download and store at specific location)

'''

sys.path.append('/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/')
from METAL_Module.metal_script import *
from LDSC_Module.ldsc_manager import * ##General_Munge

timeout_wait_setter='3600'
locii_seperator="Num_seperator_el"

cell_type_ldsc_ref="/ifs/loni/faculty/njahansh/nerds/ravi/genetics/ldsc/ldsc_seg_ldscores/"
#change it to the "-"; Mention  in the documentation
sep_btw_study_trait="___" #Seperator between study_name___trait_name.regenie 

os.chdir('/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/')
# root_dir_to_Gini=os.getcwd()#+"/"
root_dir_to_Gini_temp_files=os.getcwd()#+"/"

# root_dir_to_Gini="/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/"


# root_dir_to_Gini_temp_files="/ifs/loni/faculty/njahansh/nerds/ankush/for_Ravi/Pain_Genetics/Pain_ABCD_X_ENIGMA/"
root_dir_to_Gini=os.getcwd()#+"/"

Extra_temp_files_dict=dict( extras_folder_path = os.path.join(root_dir_to_Gini_temp_files, "Exta_temp_files"),# Check if folder exists; to add temporary folder in it
extras_metal_files = os.path.join(root_dir_to_Gini_temp_files, "Exta_temp_files/METAL_Out"),# Check if Metal extras folder exists; to add temporary files in it
extras_metal_qced_files = os.path.join(root_dir_to_Gini_temp_files, "Exta_temp_files/METAL_Out/METAL_QCED_file"),# to add QCED files in it
extras_metal_input_files = os.path.join(root_dir_to_Gini_temp_files, "Exta_temp_files/METAL_Out/METAL_input_files"),# to add .metal input files in it
extras_metal_output_files = os.path.join(root_dir_to_Gini_temp_files, "Exta_temp_files/METAL_Out/METAL_output_files"),# For METAL results
extras_LDSC_Munge_files = os.path.join(root_dir_to_Gini_temp_files, "Exta_temp_files/LDSC_Out/Munged_results"),# to add LDSC Mumnge results
extras_LDSC_Heri_files = os.path.join(root_dir_to_Gini_temp_files, "Exta_temp_files/LDSC_Out/Heritability_Results"),# to add Heritability resutls
extras_LDSC_Heri__results = os.path.join(root_dir_to_Gini_temp_files, "Exta_temp_files/LDSC_Out/Heritability_CSV_files"),# to add Heritability tables from logs
extras_LDSC_rG_results    =   os.path.join(root_dir_to_Gini_temp_files, "Exta_temp_files/LDSC_Out/rG_Results"),
rG_CSV_files    =   os.path.join(root_dir_to_Gini_temp_files, "Exta_temp_files/LDSC_Out/rG_CSV_files"),
extra_cell_h2=os.path.join(root_dir_to_Gini_temp_files, "Exta_temp_files/LDSC_Out/Cell_Type_Heritability_Out"),
extra_LAVA_filePrep=os.path.join(root_dir_to_Gini_temp_files, "Exta_temp_files/LAVA_Out/LAVA_data"),
extra_LAVA_input_files=os.path.join(root_dir_to_Gini_temp_files, "Exta_temp_files/LAVA_Out/LAVA_input"),
extra_LAVA_Matrix_files=os.path.join(root_dir_to_Gini_temp_files, "Exta_temp_files/LAVA_Out/LAVA_Matrix"),
extra_LAVA_Shell_files=os.path.join(root_dir_to_Gini_temp_files, "Exta_temp_files/LAVA_Out/LAVA_shell"),
extra_LAVA_output=os.path.join(root_dir_to_Gini_temp_files, "Exta_temp_files/LAVA_Out/LAVA_Results"),
extra_Nipype_Wf=os.path.join(root_dir_to_Gini_temp_files, "Exta_temp_files/Nipype_Wf"),
extra_CAUSE_output=os.path.join(root_dir_to_Gini_temp_files, "Exta_temp_files/CAUSE_Out"),
extra_LAVA_TWAS_output=os.path.join(root_dir_to_Gini_temp_files, "Exta_temp_files/LAVA_TWAS_Out"),
extra_LAVA_TWAS_CaseControl=os.path.join(root_dir_to_Gini_temp_files, "Exta_temp_files/LAVA_TWAS_Out/LAVA_TWAS_input"),
extra_LAVA_TWAS_Results=os.path.join(root_dir_to_Gini_temp_files, "Exta_temp_files/LAVA_TWAS_Out/LAVA_TWAS_Results")
)
                                      

rG_folder={"wSA":os.path.join(root_dir_to_Gini, "LDSC_Module/Enigma_GC_Munged/wSA/"),
           "wTHICK":os.path.join(root_dir_to_Gini, "LDSC_Module/Enigma_GC_Munged/wTHICK/"),
           "Subcortical":os.path.join(root_dir_to_Gini, "LDSC_Module/Enigma_Subcortical_Munged/")}

# Lava_cortical_subcortical={"cortical":[
#                         os.path.join(root_dir_to_Gini, "LAVA_Module/Enigma/ENIGMA3_mixed_se_wSA_Mean_bankssts_surfavg_20190429.txt.gz"),
#                         os.path.join(root_dir_to_Gini, "LAVA_Module/Enigma/ENIGMA3_mixed_se_wSA_Mean_caudalanteriorcingulate_surfavg_20190429.txt.gz"),
#                         os.path.join(root_dir_to_Gini, "LAVA_Module/Enigma/ENIGMA3_mixed_se_wSA_Mean_caudalmiddlefrontal_surfavg_20190429.txt.gz"),
#                         os.path.join(root_dir_to_Gini, "LAVA_Module/Enigma/ENIGMA3_mixed_se_wSA_Mean_cuneus_surfavg_20190429.txt.gz"),
#                         os.path.join(root_dir_to_Gini, "LAVA_Module/Enigma/ENIGMA3_mixed_se_wSA_Mean_entorhinal_surfavg_20190429.txt.gz"),
#                         os.path.join(root_dir_to_Gini, "LAVA_Module/Enigma/ENIGMA3_mixed_se_wSA_Mean_frontalpole_surfavg_20190429.txt.gz"),
#                         os.path.join(root_dir_to_Gini, "LAVA_Module/Enigma/ENIGMA3_mixed_se_wSA_Mean_fusiform_surfavg_20190429.txt.gz"),
#                         os.path.join(root_dir_to_Gini, "LAVA_Module/Enigma/ENIGMA3_mixed_se_wSA_Mean_inferiorparietal_surfavg_20190429.txt.gz"),
#                         os.path.join(root_dir_to_Gini, "LAVA_Module/Enigma/ENIGMA3_mixed_se_wSA_Mean_inferiortemporal_surfavg_20190429.txt.gz"),
#                         os.path.join(root_dir_to_Gini, "LAVA_Module/Enigma/ENIGMA3_mixed_se_wSA_Mean_insula_surfavg_20190429.txt.gz"),
#                         os.path.join(root_dir_to_Gini, "LAVA_Module/Enigma/ENIGMA3_mixed_se_wSA_Mean_isthmuscingulate_surfavg_20190429.txt.gz"),
#                         os.path.join(root_dir_to_Gini, "LAVA_Module/Enigma/ENIGMA3_mixed_se_wSA_Mean_lateraloccipital_surfavg_20190429.txt.gz"),
#                         os.path.join(root_dir_to_Gini, "LAVA_Module/Enigma/ENIGMA3_mixed_se_wSA_Mean_lateralorbitofrontal_surfavg_20190429.txt.gz"),
#                         os.path.join(root_dir_to_Gini, "LAVA_Module/Enigma/ENIGMA3_mixed_se_wSA_Mean_lingual_surfavg_20190429.txt.gz"),
#                         os.path.join(root_dir_to_Gini, "LAVA_Module/Enigma/ENIGMA3_mixed_se_wSA_Mean_medialorbitofrontal_surfavg_20190429.txt.gz"),
#                         os.path.join(root_dir_to_Gini, "LAVA_Module/Enigma/ENIGMA3_mixed_se_wSA_Mean_middletemporal_surfavg_20190429.txt.gz"),
#                         os.path.join(root_dir_to_Gini, "LAVA_Module/Enigma/ENIGMA3_mixed_se_wSA_Mean_paracentral_surfavg_20190429.txt.gz"),
#                         os.path.join(root_dir_to_Gini, "LAVA_Module/Enigma/ENIGMA3_mixed_se_wSA_Mean_parahippocampal_surfavg_20190429.txt.gz"),
#                         os.path.join(root_dir_to_Gini, "LAVA_Module/Enigma/ENIGMA3_mixed_se_wSA_Mean_parsopercularis_surfavg_20190429.txt.gz"),
#                         os.path.join(root_dir_to_Gini, "LAVA_Module/Enigma/ENIGMA3_mixed_se_wSA_Mean_parsorbitalis_surfavg_20190429.txt.gz"),
#                         os.path.join(root_dir_to_Gini, "LAVA_Module/Enigma/ENIGMA3_mixed_se_wSA_Mean_parstriangularis_surfavg_20190429.txt.gz"),
#                         os.path.join(root_dir_to_Gini, "LAVA_Module/Enigma/ENIGMA3_mixed_se_wSA_Mean_pericalcarine_surfavg_20190429.txt.gz"),
#                         os.path.join(root_dir_to_Gini, "LAVA_Module/Enigma/ENIGMA3_mixed_se_wSA_Mean_postcentral_surfavg_20190429.txt.gz"),
#                         os.path.join(root_dir_to_Gini, "LAVA_Module/Enigma/ENIGMA3_mixed_se_wSA_Mean_posteriorcingulate_surfavg_20190429.txt.gz"),
#                         os.path.join(root_dir_to_Gini, "LAVA_Module/Enigma/ENIGMA3_mixed_se_wSA_Mean_precentral_surfavg_20190429.txt.gz"),
#                         os.path.join(root_dir_to_Gini, "LAVA_Module/Enigma/ENIGMA3_mixed_se_wSA_Mean_precuneus_surfavg_20190429.txt.gz"),
#                         os.path.join(root_dir_to_Gini, "LAVA_Module/Enigma/ENIGMA3_mixed_se_wSA_Mean_rostralanteriorcingulate_surfavg_20190429.txt.gz"),
#                         os.path.join(root_dir_to_Gini, "LAVA_Module/Enigma/ENIGMA3_mixed_se_wSA_Mean_rostralmiddlefrontal_surfavg_20190429.txt.gz"),
#                         os.path.join(root_dir_to_Gini, "LAVA_Module/Enigma/ENIGMA3_mixed_se_wSA_Mean_superiorfrontal_surfavg_20190429.txt.gz"),
#                         os.path.join(root_dir_to_Gini, "LAVA_Module/Enigma/ENIGMA3_mixed_se_wSA_Mean_superiorparietal_surfavg_20190429.txt.gz"),
#                         os.path.join(root_dir_to_Gini, "LAVA_Module/Enigma/ENIGMA3_mixed_se_wSA_Mean_superiortemporal_surfavg_20190429.txt.gz"),
#                         os.path.join(root_dir_to_Gini, "LAVA_Module/Enigma/ENIGMA3_mixed_se_wSA_Mean_supramarginal_surfavg_20190429.txt.gz"),
#                         os.path.join(root_dir_to_Gini, "LAVA_Module/Enigma/ENIGMA3_mixed_se_wSA_Mean_temporalpole_surfavg_20190429.txt.gz"),
#                         os.path.join(root_dir_to_Gini, "LAVA_Module/Enigma/ENIGMA3_mixed_se_wSA_Mean_transversetemporal_surfavg_20190429.txt.gz"),
#                         os.path.join(root_dir_to_Gini, "LAVA_Module/Enigma/ENIGMA3_mixed_se_wTHICK_Mean_bankssts_thickavg_20200522.txt.gz"),
#                         os.path.join(root_dir_to_Gini, "LAVA_Module/Enigma/ENIGMA3_mixed_se_wTHICK_Mean_caudalanteriorcingulate_thickavg_20200522.txt.gz"),
#                         os.path.join(root_dir_to_Gini, "LAVA_Module/Enigma/ENIGMA3_mixed_se_wTHICK_Mean_caudalmiddlefrontal_thickavg_20200522.txt.gz"),
#                         os.path.join(root_dir_to_Gini, "LAVA_Module/Enigma/ENIGMA3_mixed_se_wTHICK_Mean_cuneus_thickavg_20200522.txt.gz"),
#                         os.path.join(root_dir_to_Gini, "LAVA_Module/Enigma/ENIGMA3_mixed_se_wTHICK_Mean_entorhinal_thickavg_20200522.txt.gz"),
#                         os.path.join(root_dir_to_Gini, "LAVA_Module/Enigma/ENIGMA3_mixed_se_wTHICK_Mean_frontalpole_thickavg_20200522.txt.gz"),
#                         os.path.join(root_dir_to_Gini, "LAVA_Module/Enigma/ENIGMA3_mixed_se_wTHICK_Mean_fusiform_thickavg_20200522.txt.gz"),
#                         os.path.join(root_dir_to_Gini, "LAVA_Module/Enigma/ENIGMA3_mixed_se_wTHICK_Mean_inferiorparietal_thickavg_20200522.txt.gz"),
#                         os.path.join(root_dir_to_Gini, "LAVA_Module/Enigma/ENIGMA3_mixed_se_wTHICK_Mean_inferiortemporal_thickavg_20200522.txt.gz"),
#                         os.path.join(root_dir_to_Gini, "LAVA_Module/Enigma/ENIGMA3_mixed_se_wTHICK_Mean_insula_thickavg_20200522.txt.gz"),
#                         os.path.join(root_dir_to_Gini, "LAVA_Module/Enigma/ENIGMA3_mixed_se_wTHICK_Mean_isthmuscingulate_thickavg_20200522.txt.gz"),
#                         os.path.join(root_dir_to_Gini, "LAVA_Module/Enigma/ENIGMA3_mixed_se_wTHICK_Mean_lateraloccipital_thickavg_20200522.txt.gz"),
#                         os.path.join(root_dir_to_Gini, "LAVA_Module/Enigma/ENIGMA3_mixed_se_wTHICK_Mean_lateralorbitofrontal_thickavg_20200522.txt.gz"),
#                         os.path.join(root_dir_to_Gini, "LAVA_Module/Enigma/ENIGMA3_mixed_se_wTHICK_Mean_lingual_thickavg_20200522.txt.gz"),
#                         os.path.join(root_dir_to_Gini, "LAVA_Module/Enigma/ENIGMA3_mixed_se_wTHICK_Mean_medialorbitofrontal_thickavg_20200522.txt.gz"),
#                         os.path.join(root_dir_to_Gini, "LAVA_Module/Enigma/ENIGMA3_mixed_se_wTHICK_Mean_middletemporal_thickavg_20200522.txt.gz"),
#                         os.path.join(root_dir_to_Gini, "LAVA_Module/Enigma/ENIGMA3_mixed_se_wTHICK_Mean_paracentral_thickavg_20200522.txt.gz"),
#                         os.path.join(root_dir_to_Gini, "LAVA_Module/Enigma/ENIGMA3_mixed_se_wTHICK_Mean_parahippocampal_thickavg_20200522.txt.gz"),
#                         os.path.join(root_dir_to_Gini, "LAVA_Module/Enigma/ENIGMA3_mixed_se_wTHICK_Mean_parsopercularis_thickavg_20200522.txt.gz"),
#                         os.path.join(root_dir_to_Gini, "LAVA_Module/Enigma/ENIGMA3_mixed_se_wTHICK_Mean_parsorbitalis_thickavg_20200522.txt.gz"),
#                         os.path.join(root_dir_to_Gini, "LAVA_Module/Enigma/ENIGMA3_mixed_se_wTHICK_Mean_parstriangularis_thickavg_20200522.txt.gz"),
#                         os.path.join(root_dir_to_Gini, "LAVA_Module/Enigma/ENIGMA3_mixed_se_wTHICK_Mean_pericalcarine_thickavg_20200522.txt.gz"),
#                         os.path.join(root_dir_to_Gini, "LAVA_Module/Enigma/ENIGMA3_mixed_se_wTHICK_Mean_postcentral_thickavg_20200522.txt.gz"),
#                         os.path.join(root_dir_to_Gini, "LAVA_Module/Enigma/ENIGMA3_mixed_se_wTHICK_Mean_posteriorcingulate_thickavg_20200522.txt.gz"),
#                         os.path.join(root_dir_to_Gini, "LAVA_Module/Enigma/ENIGMA3_mixed_se_wTHICK_Mean_precentral_thickavg_20200522.txt.gz"),
#                         os.path.join(root_dir_to_Gini, "LAVA_Module/Enigma/ENIGMA3_mixed_se_wTHICK_Mean_precuneus_thickavg_20200522.txt.gz"),
#                         os.path.join(root_dir_to_Gini, "LAVA_Module/Enigma/ENIGMA3_mixed_se_wTHICK_Mean_rostralanteriorcingulate_thickavg_20200522.txt.gz"),
#                         os.path.join(root_dir_to_Gini, "LAVA_Module/Enigma/ENIGMA3_mixed_se_wTHICK_Mean_rostralmiddlefrontal_thickavg_20200522.txt.gz"),
#                         os.path.join(root_dir_to_Gini, "LAVA_Module/Enigma/ENIGMA3_mixed_se_wTHICK_Mean_superiorfrontal_thickavg_20200522.txt.gz"),
#                         os.path.join(root_dir_to_Gini, "LAVA_Module/Enigma/ENIGMA3_mixed_se_wTHICK_Mean_superiorparietal_thickavg_20200522.txt.gz"),
#                         os.path.join(root_dir_to_Gini, "LAVA_Module/Enigma/ENIGMA3_mixed_se_wTHICK_Mean_superiortemporal_thickavg_20200522.txt.gz"),
#                         os.path.join(root_dir_to_Gini, "LAVA_Module/Enigma/ENIGMA3_mixed_se_wTHICK_Mean_supramarginal_thickavg_20200522.txt.gz"),
#                         os.path.join(root_dir_to_Gini, "LAVA_Module/Enigma/ENIGMA3_mixed_se_wTHICK_Mean_temporalpole_thickavg_20200522.txt.gz"),
#                         os.path.join(root_dir_to_Gini, "LAVA_Module/Enigma/ENIGMA3_mixed_se_wTHICK_Mean_transversetemporal_thickavg_20200522.txt.gz")
#                         ],
#                            "subcortical":"/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/Exta_temp_files/INPUT_TESTING_PACKAGE/Enigma_SE_GC/"}


association_mapped_CAUSE=["SNPID","BETA","SE", "EA","NEA","P"]

Cortical_directory=os.path.join(root_dir_to_Gini, "LAVA_Module/Enigma/")

Lava_cortical_subcortical={"cortical":{
                            "wSA":[
                                os.path.join(root_dir_to_Gini, "LAVA_Module/Enigma/ENIGMA3_mixed_se_wSA_Mean_bankssts_surfavg_20190429.txt.gz"),
                                os.path.join(root_dir_to_Gini, "LAVA_Module/Enigma/ENIGMA3_mixed_se_wSA_Mean_caudalanteriorcingulate_surfavg_20190429.txt.gz"),
                                os.path.join(root_dir_to_Gini, "LAVA_Module/Enigma/ENIGMA3_mixed_se_wSA_Mean_caudalmiddlefrontal_surfavg_20190429.txt.gz"),
                                os.path.join(root_dir_to_Gini, "LAVA_Module/Enigma/ENIGMA3_mixed_se_wSA_Mean_cuneus_surfavg_20190429.txt.gz"),
                                os.path.join(root_dir_to_Gini, "LAVA_Module/Enigma/ENIGMA3_mixed_se_wSA_Mean_entorhinal_surfavg_20190429.txt.gz"),
                                os.path.join(root_dir_to_Gini, "LAVA_Module/Enigma/ENIGMA3_mixed_se_wSA_Mean_frontalpole_surfavg_20190429.txt.gz"),
                                os.path.join(root_dir_to_Gini, "LAVA_Module/Enigma/ENIGMA3_mixed_se_wSA_Mean_fusiform_surfavg_20190429.txt.gz"),
                                os.path.join(root_dir_to_Gini, "LAVA_Module/Enigma/ENIGMA3_mixed_se_wSA_Mean_inferiorparietal_surfavg_20190429.txt.gz"),
                                os.path.join(root_dir_to_Gini, "LAVA_Module/Enigma/ENIGMA3_mixed_se_wSA_Mean_inferiortemporal_surfavg_20190429.txt.gz"),
                                os.path.join(root_dir_to_Gini, "LAVA_Module/Enigma/ENIGMA3_mixed_se_wSA_Mean_insula_surfavg_20190429.txt.gz"),
                                os.path.join(root_dir_to_Gini, "LAVA_Module/Enigma/ENIGMA3_mixed_se_wSA_Mean_isthmuscingulate_surfavg_20190429.txt.gz"),
                                os.path.join(root_dir_to_Gini, "LAVA_Module/Enigma/ENIGMA3_mixed_se_wSA_Mean_lateraloccipital_surfavg_20190429.txt.gz"),
                                os.path.join(root_dir_to_Gini, "LAVA_Module/Enigma/ENIGMA3_mixed_se_wSA_Mean_lateralorbitofrontal_surfavg_20190429.txt.gz"),
                                os.path.join(root_dir_to_Gini, "LAVA_Module/Enigma/ENIGMA3_mixed_se_wSA_Mean_lingual_surfavg_20190429.txt.gz"),
                                os.path.join(root_dir_to_Gini, "LAVA_Module/Enigma/ENIGMA3_mixed_se_wSA_Mean_medialorbitofrontal_surfavg_20190429.txt.gz"),
                                os.path.join(root_dir_to_Gini, "LAVA_Module/Enigma/ENIGMA3_mixed_se_wSA_Mean_middletemporal_surfavg_20190429.txt.gz"),
                                os.path.join(root_dir_to_Gini, "LAVA_Module/Enigma/ENIGMA3_mixed_se_wSA_Mean_paracentral_surfavg_20190429.txt.gz"),
                                os.path.join(root_dir_to_Gini, "LAVA_Module/Enigma/ENIGMA3_mixed_se_wSA_Mean_parahippocampal_surfavg_20190429.txt.gz"),
                                os.path.join(root_dir_to_Gini, "LAVA_Module/Enigma/ENIGMA3_mixed_se_wSA_Mean_parsopercularis_surfavg_20190429.txt.gz"),
                                os.path.join(root_dir_to_Gini, "LAVA_Module/Enigma/ENIGMA3_mixed_se_wSA_Mean_parsorbitalis_surfavg_20190429.txt.gz"),
                                os.path.join(root_dir_to_Gini, "LAVA_Module/Enigma/ENIGMA3_mixed_se_wSA_Mean_parstriangularis_surfavg_20190429.txt.gz"),
                                os.path.join(root_dir_to_Gini, "LAVA_Module/Enigma/ENIGMA3_mixed_se_wSA_Mean_pericalcarine_surfavg_20190429.txt.gz"),
                                os.path.join(root_dir_to_Gini, "LAVA_Module/Enigma/ENIGMA3_mixed_se_wSA_Mean_postcentral_surfavg_20190429.txt.gz"),
                                os.path.join(root_dir_to_Gini, "LAVA_Module/Enigma/ENIGMA3_mixed_se_wSA_Mean_posteriorcingulate_surfavg_20190429.txt.gz"),
                                os.path.join(root_dir_to_Gini, "LAVA_Module/Enigma/ENIGMA3_mixed_se_wSA_Mean_precentral_surfavg_20190429.txt.gz"),
                                os.path.join(root_dir_to_Gini, "LAVA_Module/Enigma/ENIGMA3_mixed_se_wSA_Mean_precuneus_surfavg_20190429.txt.gz"),
                                os.path.join(root_dir_to_Gini, "LAVA_Module/Enigma/ENIGMA3_mixed_se_wSA_Mean_rostralanteriorcingulate_surfavg_20190429.txt.gz"),
                                os.path.join(root_dir_to_Gini, "LAVA_Module/Enigma/ENIGMA3_mixed_se_wSA_Mean_rostralmiddlefrontal_surfavg_20190429.txt.gz"),
                                os.path.join(root_dir_to_Gini, "LAVA_Module/Enigma/ENIGMA3_mixed_se_wSA_Mean_superiorfrontal_surfavg_20190429.txt.gz"),
                                os.path.join(root_dir_to_Gini, "LAVA_Module/Enigma/ENIGMA3_mixed_se_wSA_Mean_superiorparietal_surfavg_20190429.txt.gz"),
                                os.path.join(root_dir_to_Gini, "LAVA_Module/Enigma/ENIGMA3_mixed_se_wSA_Mean_superiortemporal_surfavg_20190429.txt.gz"),
                                os.path.join(root_dir_to_Gini, "LAVA_Module/Enigma/ENIGMA3_mixed_se_wSA_Mean_supramarginal_surfavg_20190429.txt.gz"),
                                os.path.join(root_dir_to_Gini, "LAVA_Module/Enigma/ENIGMA3_mixed_se_wSA_Mean_temporalpole_surfavg_20190429.txt.gz"),
                                os.path.join(root_dir_to_Gini, "LAVA_Module/Enigma/ENIGMA3_mixed_se_wSA_Mean_transversetemporal_surfavg_20190429.txt.gz")
                            ], 
                            "wTHICK":[
                                os.path.join(root_dir_to_Gini, "LAVA_Module/Enigma/ENIGMA3_mixed_se_wTHICK_Mean_bankssts_thickavg_20200522.txt.gz"),
                                os.path.join(root_dir_to_Gini, "LAVA_Module/Enigma/ENIGMA3_mixed_se_wTHICK_Mean_caudalanteriorcingulate_thickavg_20200522.txt.gz"),
                                os.path.join(root_dir_to_Gini, "LAVA_Module/Enigma/ENIGMA3_mixed_se_wTHICK_Mean_caudalmiddlefrontal_thickavg_20200522.txt.gz"),
                                os.path.join(root_dir_to_Gini, "LAVA_Module/Enigma/ENIGMA3_mixed_se_wTHICK_Mean_cuneus_thickavg_20200522.txt.gz"),
                                os.path.join(root_dir_to_Gini, "LAVA_Module/Enigma/ENIGMA3_mixed_se_wTHICK_Mean_entorhinal_thickavg_20200522.txt.gz"),
                                os.path.join(root_dir_to_Gini, "LAVA_Module/Enigma/ENIGMA3_mixed_se_wTHICK_Mean_frontalpole_thickavg_20200522.txt.gz"),
                                os.path.join(root_dir_to_Gini, "LAVA_Module/Enigma/ENIGMA3_mixed_se_wTHICK_Mean_fusiform_thickavg_20200522.txt.gz"),
                                os.path.join(root_dir_to_Gini, "LAVA_Module/Enigma/ENIGMA3_mixed_se_wTHICK_Mean_inferiorparietal_thickavg_20200522.txt.gz"),
                                os.path.join(root_dir_to_Gini, "LAVA_Module/Enigma/ENIGMA3_mixed_se_wTHICK_Mean_inferiortemporal_thickavg_20200522.txt.gz"),
                                os.path.join(root_dir_to_Gini, "LAVA_Module/Enigma/ENIGMA3_mixed_se_wTHICK_Mean_insula_thickavg_20200522.txt.gz"),
                                os.path.join(root_dir_to_Gini, "LAVA_Module/Enigma/ENIGMA3_mixed_se_wTHICK_Mean_isthmuscingulate_thickavg_20200522.txt.gz"),
                                os.path.join(root_dir_to_Gini, "LAVA_Module/Enigma/ENIGMA3_mixed_se_wTHICK_Mean_lateraloccipital_thickavg_20200522.txt.gz"),
                                os.path.join(root_dir_to_Gini, "LAVA_Module/Enigma/ENIGMA3_mixed_se_wTHICK_Mean_lateralorbitofrontal_thickavg_20200522.txt.gz"),
                                os.path.join(root_dir_to_Gini, "LAVA_Module/Enigma/ENIGMA3_mixed_se_wTHICK_Mean_lingual_thickavg_20200522.txt.gz"),
                                os.path.join(root_dir_to_Gini, "LAVA_Module/Enigma/ENIGMA3_mixed_se_wTHICK_Mean_medialorbitofrontal_thickavg_20200522.txt.gz"),
                                os.path.join(root_dir_to_Gini, "LAVA_Module/Enigma/ENIGMA3_mixed_se_wTHICK_Mean_middletemporal_thickavg_20200522.txt.gz"),
                                os.path.join(root_dir_to_Gini, "LAVA_Module/Enigma/ENIGMA3_mixed_se_wTHICK_Mean_paracentral_thickavg_20200522.txt.gz"),
                                os.path.join(root_dir_to_Gini, "LAVA_Module/Enigma/ENIGMA3_mixed_se_wTHICK_Mean_parahippocampal_thickavg_20200522.txt.gz"),
                                os.path.join(root_dir_to_Gini, "LAVA_Module/Enigma/ENIGMA3_mixed_se_wTHICK_Mean_parsopercularis_thickavg_20200522.txt.gz"),
                                os.path.join(root_dir_to_Gini, "LAVA_Module/Enigma/ENIGMA3_mixed_se_wTHICK_Mean_parsorbitalis_thickavg_20200522.txt.gz"),
                                os.path.join(root_dir_to_Gini, "LAVA_Module/Enigma/ENIGMA3_mixed_se_wTHICK_Mean_parstriangularis_thickavg_20200522.txt.gz"),
                                os.path.join(root_dir_to_Gini, "LAVA_Module/Enigma/ENIGMA3_mixed_se_wTHICK_Mean_pericalcarine_thickavg_20200522.txt.gz"),
                                os.path.join(root_dir_to_Gini, "LAVA_Module/Enigma/ENIGMA3_mixed_se_wTHICK_Mean_postcentral_thickavg_20200522.txt.gz"),
                                os.path.join(root_dir_to_Gini, "LAVA_Module/Enigma/ENIGMA3_mixed_se_wTHICK_Mean_posteriorcingulate_thickavg_20200522.txt.gz"),
                                os.path.join(root_dir_to_Gini, "LAVA_Module/Enigma/ENIGMA3_mixed_se_wTHICK_Mean_precentral_thickavg_20200522.txt.gz"),
                                os.path.join(root_dir_to_Gini, "LAVA_Module/Enigma/ENIGMA3_mixed_se_wTHICK_Mean_precuneus_thickavg_20200522.txt.gz"),
                                os.path.join(root_dir_to_Gini, "LAVA_Module/Enigma/ENIGMA3_mixed_se_wTHICK_Mean_rostralanteriorcingulate_thickavg_20200522.txt.gz"),
                                os.path.join(root_dir_to_Gini, "LAVA_Module/Enigma/ENIGMA3_mixed_se_wTHICK_Mean_rostralmiddlefrontal_thickavg_20200522.txt.gz"),
                                os.path.join(root_dir_to_Gini, "LAVA_Module/Enigma/ENIGMA3_mixed_se_wTHICK_Mean_superiorfrontal_thickavg_20200522.txt.gz"),
                                os.path.join(root_dir_to_Gini, "LAVA_Module/Enigma/ENIGMA3_mixed_se_wTHICK_Mean_superiorparietal_thickavg_20200522.txt.gz"),
                                os.path.join(root_dir_to_Gini, "LAVA_Module/Enigma/ENIGMA3_mixed_se_wTHICK_Mean_superiortemporal_thickavg_20200522.txt.gz"),
                                os.path.join(root_dir_to_Gini, "LAVA_Module/Enigma/ENIGMA3_mixed_se_wTHICK_Mean_supramarginal_thickavg_20200522.txt.gz"),
                                os.path.join(root_dir_to_Gini, "LAVA_Module/Enigma/ENIGMA3_mixed_se_wTHICK_Mean_temporalpole_thickavg_20200522.txt.gz"),
                                os.path.join(root_dir_to_Gini, "LAVA_Module/Enigma/ENIGMA3_mixed_se_wTHICK_Mean_transversetemporal_thickavg_20200522.txt.gz")
                        ]},
                           "subcortical":"/ifs/loni/faculty/njahansh/nerds/ravi/genetics/enigma_subcortical_gwas"}


def file_name_process(fileName:str):
    return fileName.replace(".tbl","").replace(".csv","").replace(".txt","").replace(".gz","").replace(".tsv","").replace(".vcf","").split("/")[-1].split(".")[0]

# with open('/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/association_format.json', 'r') as json_file:
#     json_association = json.load(json_file)

with open('/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/LAVA_formatted_Association.json', 'r') as json_file:
    json_association = json.load(json_file)

file_joiner_str="abcdFileX123__"