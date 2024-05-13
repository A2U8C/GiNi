import pandas as pd



def CAUSE_run(input_trait:str, Enigma:str, bfile:str, input_trait_cols:str):
    import os
    import sys
    import subprocess
        
    sys.path.append('/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/')
    import CONSTANTS
     
    trait_name=input_trait.split("___")[1].split(".")[0].lstrip("_")
    # trait_name=input_trait.split(CONSTANTS.sep_btw_study_trait)[1].split(".")[0].lstrip("_")
    Enigma_trait=Enigma.split("ENIGMA3_mixed_se_")[-1].split("_")[2]
    
    wSA_wTHICK=Enigma.split("ENIGMA3_mixed_se_")[-1].split("_")[0]
    # output_folder="/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/Exta_temp_files/CAUSE_Out"+"/"+trait_name+"/"+wSA_wTHICK
    # output_folder="/ifs/loni/faculty/njahansh/nerds/ankush/for_Ravi/Pain_Genetics/Pain_ABCD_X_ENIGMA/Exta_temp_files/CAUSE_Out"+"/"+trait_name+"/"+wSA_wTHICK
    output_folder=CONSTANTS.Extra_temp_files_dict["extra_CAUSE_output"]+"/"+trait_name+"/"+wSA_wTHICK

    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
        print("Nested folders created successfully.")

    p_val_threshold=0.001#0.05/68;
    r_command=f'Rscript /ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/CAUSE_Module/Cause.R "{input_trait}" "{Enigma}" "{bfile}" "{input_trait_cols}" "{trait_name}__{Enigma_trait}" "{output_folder}" {p_val_threshold}'


    print(r_command)

    
    # Define the environment variables specific to your R environment
    r_environment_path = "/ifs/loni/faculty/njahansh/nerds/ankush/software/Anaconda/envs/r_env_ankush/bin"#"/ifs/loni/faculty/njahansh/nerds/shruti/software/python/anaconda3/envs/r_env/bin"

    r_environment = os.environ.copy()
    r_environment["PATH"] = f"{r_environment_path}:{r_environment['PATH']}"

    # Execute the command
    subprocess.run(r_command, env=r_environment,shell=True)

    return output_folder



def CAUSE_prep_function(rG_csv_file_path:str, gwasFormat:str,dict_munge_raw:dict,machineName):
    import pandas as pd
    import sys
    
    sys.path.append('/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/')
    import CONSTANTS
    import nipype
    import nipype.pipeline.engine as pe
    import nipype.interfaces.utility as niu
    from nipype.interfaces.utility import Function
    from nipype import config

    df_rG_csv_file = pd.read_csv(rG_csv_file_path)

    interate_list=[] #[input_trait,Enigma]
    thresh_freq=0.05/68#df_rG_csv_file.shape[0]#0.05
    column_names_cause=[]

    dict_col_header_mapper={v:k for k,v in CONSTANTS.json_association[gwasFormat]['format_dict'].items()}


    for el in CONSTANTS.association_mapped_CAUSE:
        column_names_cause.append(dict_col_header_mapper[el])



    for index, row in df_rG_csv_file.iterrows():
        if row['p']<thresh_freq:
            interate_list.append([dict_munge_raw[row["p1"]],CONSTANTS.Cortical_directory+"/"+row["p2"].split("/")[-1].split(".")[0]+".txt.gz"])

    

    wf_name='Cause_wf_Significants'
    workflow_metal_mungng = pe.Workflow(name=wf_name)
    workflow_metal_mungng.base_dir = CONSTANTS.Extra_temp_files_dict["extra_Nipype_Wf"]+"/"#"/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/Exta_temp_files/"
    workflow_metal_mungng.config["job_finished_timeout"]="60"

        
    config_dict={'execution': {'job_finished_timeout': '60'}}
    config.update_config(config_dict)
    print(interate_list)
    inputnode = pe.Node(
        niu.IdentityInterface(
            fields=['input_trait', 'Enigma']
            ),
        synchronize=True,
        iterables=[('input_trait', 'Enigma'),interate_list],
        name='inputnode')
        
    node1_LDSC_rG = pe.Node(Function(input_names=['input_trait', 'Enigma','bfile', 'input_trait_cols'],
                    output_names=['out_file_comp'],
                    function=CAUSE_run),
            name='node1_LDSC_rG')
    node1_LDSC_rG.inputs.bfile="/ifs/loni/faculty/njahansh/nerds/ravi/genetics/LAVA-main/support_data/eur/g1000_eur"
    node1_LDSC_rG.inputs.input_trait_cols=" ".join(column_names_cause)

    
    connections=[(inputnode, node1_LDSC_rG, [('input_trait', 'input_trait')]),
                 (inputnode, node1_LDSC_rG, [('Enigma', 'Enigma')])
                 ] 
        

    workflow_metal_mungng.connect(connections)    
    res=workflow_metal_mungng.run('SGE',plugin_args=dict(dont_resubmit_completed_jobs= True, overwrite= True, qsub_args=f'-q {machineName} -o /ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/Exta_temp_files/Logs/ -j y'))


# CAUSE_run("/ifs/loni/faculty/njahansh/nerds/ankush/for_Ravi/Pain_Genetics/Pain_Inputs_GiNi_Formatted/UKBB_HUNT___COPC.csv","/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/LAVA_Module/Enigma//ENIGMA3_mixed_se_wTHICK_Mean_superiorfrontal_thickavg_20200522.txt.gz","/ifs/loni/faculty/njahansh/nerds/ravi/genetics/LAVA-main/support_data/eur/g1000_eur","SNP BETA SE A1 A2 PVAL")


# all_arr=[['/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/Exta_temp_files/METAL_Out/METAL_output_files/ukbb_white_abcd_white____Total_MeanThickness1.tbl.gz', '/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/LAVA_Module/Enigma//ENIGMA3_mixed_se_wSA_Mean_inferiortemporal_surfavg_20190429.txt.gz'], ['/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/Exta_temp_files/METAL_Out/METAL_output_files/ukbb_white_abcd_white____Total_MeanThickness1.tbl.gz', '/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/LAVA_Module/Enigma//ENIGMA3_mixed_se_wSA_Mean_lateraloccipital_surfavg_20190429.txt.gz'], ['/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/Exta_temp_files/METAL_Out/METAL_output_files/ukbb_white_abcd_white____Total_MeanThickness1.tbl.gz', '/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/LAVA_Module/Enigma//ENIGMA3_mixed_se_wSA_Mean_paracentral_surfavg_20190429.txt.gz'], ['/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/Exta_temp_files/METAL_Out/METAL_output_files/ukbb_white_abcd_white____Total_MeanThickness1.tbl.gz', '/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/LAVA_Module/Enigma//ENIGMA3_mixed_se_wSA_Mean_rostralanteriorcingulate_surfavg_20190429.txt.gz'], ['/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/Exta_temp_files/METAL_Out/METAL_output_files/ukbb_white_abcd_white____Total_MeanThickness1.tbl.gz', '/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/LAVA_Module/Enigma//ENIGMA3_mixed_se_wSA_Mean_caudalanteriorcingulate_surfavg_20190429.txt.gz'], ['/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/Exta_temp_files/METAL_Out/METAL_output_files/ukbb_white_abcd_white____Total_MeanThickness1.tbl.gz', '/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/LAVA_Module/Enigma//ENIGMA3_mixed_se_wSA_Mean_medialorbitofrontal_surfavg_20190429.txt.gz'], ['/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/Exta_temp_files/METAL_Out/METAL_output_files/ukbb_white_abcd_white____Total_MeanThickness1.tbl.gz', '/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/LAVA_Module/Enigma//ENIGMA3_mixed_se_wSA_Mean_precuneus_surfavg_20190429.txt.gz'], ['/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/Exta_temp_files/METAL_Out/METAL_output_files/ukbb_white_abcd_white____Total_MeanThickness1.tbl.gz', '/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/LAVA_Module/Enigma//ENIGMA3_mixed_se_wTHICK_Mean_inferiorparietal_thickavg_20200522.txt.gz'], ['/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/Exta_temp_files/METAL_Out/METAL_output_files/ukbb_white_abcd_white____Total_MeanThickness1.tbl.gz', '/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/LAVA_Module/Enigma//ENIGMA3_mixed_se_wTHICK_Mean_cuneus_thickavg_20200522.txt.gz'], ['/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/Exta_temp_files/METAL_Out/METAL_output_files/ukbb_white_abcd_white____Total_MeanThickness1.tbl.gz', '/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/LAVA_Module/Enigma//ENIGMA3_mixed_se_wTHICK_Mean_caudalanteriorcingulate_thickavg_20200522.txt.gz'], ['/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/Exta_temp_files/METAL_Out/METAL_output_files/ukbb_white_abcd_white____Total_MeanThickness1.tbl.gz', '/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/LAVA_Module/Enigma//ENIGMA3_mixed_se_wTHICK_Mean_rostralanteriorcingulate_thickavg_20200522.txt.gz'], ['/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/Exta_temp_files/METAL_Out/METAL_output_files/ukbb_white_abcd_white____Total_MeanThickness1.tbl.gz', '/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/LAVA_Module/Enigma//ENIGMA3_mixed_se_wTHICK_Mean_posteriorcingulate_thickavg_20200522.txt.gz'], ['/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/Exta_temp_files/METAL_Out/METAL_output_files/ukbb_white_abcd_white____Total_MeanThickness1.tbl.gz', '/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/LAVA_Module/Enigma//ENIGMA3_mixed_se_wTHICK_Mean_lingual_thickavg_20200522.txt.gz'], ['/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/Exta_temp_files/METAL_Out/METAL_output_files/ukbb_white_abcd_white____Total_MeanThickness1.tbl.gz', '/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/LAVA_Module/Enigma//ENIGMA3_mixed_se_wTHICK_Mean_fusiform_thickavg_20200522.txt.gz'], ['/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/Exta_temp_files/METAL_Out/METAL_output_files/ukbb_white_abcd_white____Total_MeanThickness1.tbl.gz', '/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/LAVA_Module/Enigma//ENIGMA3_mixed_se_wTHICK_Mean_precentral_thickavg_20200522.txt.gz'], ['/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/Exta_temp_files/METAL_Out/METAL_output_files/ukbb_white_abcd_white____Total_MeanThickness1.tbl.gz', '/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/LAVA_Module/Enigma//ENIGMA3_mixed_se_wTHICK_Mean_medialorbitofrontal_thickavg_20200522.txt.gz'], ['/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/Exta_temp_files/METAL_Out/METAL_output_files/ukbb_white_abcd_white____Witelson5_Genu_Area1.tbl.gz', '/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/LAVA_Module/Enigma//ENIGMA3_mixed_se_wSA_Mean_frontalpole_surfavg_20190429.txt.gz'], ['/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/Exta_temp_files/METAL_Out/METAL_output_files/ukbb_white_abcd_white____Witelson5_Genu_Area1.tbl.gz', '/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/LAVA_Module/Enigma//ENIGMA3_mixed_se_wSA_Mean_rostralmiddlefrontal_surfavg_20190429.txt.gz'], ['/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/Exta_temp_files/METAL_Out/METAL_output_files/ukbb_white_abcd_white____Witelson5_Genu_Area1.tbl.gz', '/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/LAVA_Module/Enigma//ENIGMA3_mixed_se_wSA_Mean_precuneus_surfavg_20190429.txt.gz'], ['/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/Exta_temp_files/METAL_Out/METAL_output_files/ukbb_white_abcd_white____Witelson5_Genu_Area1.tbl.gz', '/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/LAVA_Module/Enigma//ENIGMA3_mixed_se_wTHICK_Mean_cuneus_thickavg_20200522.txt.gz'], ['/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/Exta_temp_files/METAL_Out/METAL_output_files/ukbb_white_abcd_white____Witelson5_Genu_Area1.tbl.gz', '/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/LAVA_Module/Enigma//ENIGMA3_mixed_se_wTHICK_Mean_rostralanteriorcingulate_thickavg_20200522.txt.gz'], ['/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/Exta_temp_files/METAL_Out/METAL_output_files/ukbb_white_abcd_white____Witelson5_Genu_Area1.tbl.gz', '/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/LAVA_Module/Enigma//ENIGMA3_mixed_se_wTHICK_Mean_posteriorcingulate_thickavg_20200522.txt.gz'], ['/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/Exta_temp_files/METAL_Out/METAL_output_files/ukbb_white_abcd_white____Witelson5_Genu_Area1.tbl.gz', '/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/LAVA_Module/Enigma//ENIGMA3_mixed_se_wTHICK_Mean_lingual_thickavg_20200522.txt.gz'], ['/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/Exta_temp_files/METAL_Out/METAL_output_files/ukbb_white_abcd_white____Witelson5_Genu_Area1.tbl.gz', '/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/LAVA_Module/Enigma//ENIGMA3_mixed_se_wTHICK_Mean_middletemporal_thickavg_20200522.txt.gz'], ['/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/Exta_temp_files/METAL_Out/METAL_output_files/ukbb_white_abcd_white____Witelson5_Genu_Area1.tbl.gz', '/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/LAVA_Module/Enigma//ENIGMA3_mixed_se_wTHICK_Mean_medialorbitofrontal_thickavg_20200522.txt.gz'], ['/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/Exta_temp_files/METAL_Out/METAL_output_files/ukbb_white_abcd_white____Witelson5_Isthmus_MeanThickness1.tbl.gz', '/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/LAVA_Module/Enigma//ENIGMA3_mixed_se_wSA_Mean_superiorfrontal_surfavg_20190429.txt.gz'], ['/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/Exta_temp_files/METAL_Out/METAL_output_files/ukbb_white_abcd_white____Witelson5_Isthmus_MeanThickness1.tbl.gz', '/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/LAVA_Module/Enigma//ENIGMA3_mixed_se_wSA_Mean_middletemporal_surfavg_20190429.txt.gz'], ['/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/Exta_temp_files/METAL_Out/METAL_output_files/ukbb_white_abcd_white____Witelson5_Isthmus_MeanThickness1.tbl.gz', '/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/LAVA_Module/Enigma//ENIGMA3_mixed_se_wSA_Mean_cuneus_surfavg_20190429.txt.gz'], ['/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/Exta_temp_files/METAL_Out/METAL_output_files/ukbb_white_abcd_white____Witelson5_Isthmus_MeanThickness1.tbl.gz', '/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/LAVA_Module/Enigma//ENIGMA3_mixed_se_wSA_Mean_inferiortemporal_surfavg_20190429.txt.gz'], ['/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/Exta_temp_files/METAL_Out/METAL_output_files/ukbb_white_abcd_white____Witelson5_Isthmus_MeanThickness1.tbl.gz', '/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/LAVA_Module/Enigma//ENIGMA3_mixed_se_wSA_Mean_lateralorbitofrontal_surfavg_20190429.txt.gz'], ['/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/Exta_temp_files/METAL_Out/METAL_output_files/ukbb_white_abcd_white____Witelson5_Isthmus_MeanThickness1.tbl.gz', '/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/LAVA_Module/Enigma//ENIGMA3_mixed_se_wSA_Mean_lateraloccipital_surfavg_20190429.txt.gz'], ['/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/Exta_temp_files/METAL_Out/METAL_output_files/ukbb_white_abcd_white____Witelson5_Isthmus_MeanThickness1.tbl.gz', '/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/LAVA_Module/Enigma//ENIGMA3_mixed_se_wSA_Mean_paracentral_surfavg_20190429.txt.gz'], ['/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/Exta_temp_files/METAL_Out/METAL_output_files/ukbb_white_abcd_white____Witelson5_Isthmus_MeanThickness1.tbl.gz', '/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/LAVA_Module/Enigma//ENIGMA3_mixed_se_wSA_Mean_rostralanteriorcingulate_surfavg_20190429.txt.gz'], ['/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/Exta_temp_files/METAL_Out/METAL_output_files/ukbb_white_abcd_white____Witelson5_Isthmus_MeanThickness1.tbl.gz', '/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/LAVA_Module/Enigma//ENIGMA3_mixed_se_wSA_Mean_medialorbitofrontal_surfavg_20190429.txt.gz'], ['/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/Exta_temp_files/METAL_Out/METAL_output_files/ukbb_white_abcd_white____Witelson5_Isthmus_MeanThickness1.tbl.gz', '/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/LAVA_Module/Enigma//ENIGMA3_mixed_se_wSA_Mean_precuneus_surfavg_20190429.txt.gz'], ['/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/Exta_temp_files/METAL_Out/METAL_output_files/ukbb_white_abcd_white____Witelson5_Isthmus_MeanThickness1.tbl.gz', '/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/LAVA_Module/Enigma//ENIGMA3_mixed_se_wTHICK_Mean_rostralanteriorcingulate_thickavg_20200522.txt.gz'], ['/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/Exta_temp_files/METAL_Out/METAL_output_files/ukbb_white_abcd_white____Witelson5_Isthmus_MeanThickness1.tbl.gz', '/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/LAVA_Module/Enigma//ENIGMA3_mixed_se_wTHICK_Mean_caudalmiddlefrontal_thickavg_20200522.txt.gz'], ['/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/Exta_temp_files/METAL_Out/METAL_output_files/ukbb_white_abcd_white____Witelson5_Isthmus_MeanThickness1.tbl.gz', '/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/LAVA_Module/Enigma//ENIGMA3_mixed_se_wTHICK_Mean_posteriorcingulate_thickavg_20200522.txt.gz'], ['/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/Exta_temp_files/METAL_Out/METAL_output_files/ukbb_white_abcd_white____Witelson5_Isthmus_MeanThickness1.tbl.gz', '/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/LAVA_Module/Enigma//ENIGMA3_mixed_se_wTHICK_Mean_lingual_thickavg_20200522.txt.gz'], ['/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/Exta_temp_files/METAL_Out/METAL_output_files/ukbb_white_abcd_white____Witelson5_Isthmus_MeanThickness1.tbl.gz', '/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/LAVA_Module/Enigma//ENIGMA3_mixed_se_wTHICK_Mean_precentral_thickavg_20200522.txt.gz'], ['/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/Exta_temp_files/METAL_Out/METAL_output_files/ukbb_white_abcd_white____Witelson5_Isthmus_MeanThickness1.tbl.gz', '/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/LAVA_Module/Enigma//ENIGMA3_mixed_se_wTHICK_Mean_postcentral_thickavg_20200522.txt.gz']]
# CAUSE_run(all_arr[0][0],all_arr[0][1],"/ifs/loni/faculty/njahansh/nerds/ravi/genetics/LAVA-main/support_data/eur/g1000_eur","MarkerName Effect StdErr Allele1 Allele2 PvalueARE ")


