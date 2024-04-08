

def Lava_FilePrep(filename:str,Gwas_format:str):
    import pandas as pd
    import numpy as np
    import re
    import subprocess
    import os
    import sys

    Gwas_format=Gwas_format.lower()
    
    sys.path.append('/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/')
    import CONSTANTS

    abs_path_data=CONSTANTS.Extra_temp_files_dict["extra_LAVA_filePrep"]

    df=pd.read_csv(filename,delim_whitespace=True)
    
    
    dict_col_header_mapper=CONSTANTS.json_association[Gwas_format]['LAVA_Formatted_key']
    print(Gwas_format)
    df.rename(columns=dict_col_header_mapper,inplace=True)

    if Gwas_format=="regenie":    
        df["P"]= (1/10)**df['LOG10P']

    df["A2"]=df["A2"].str.upper()
    df["A1"]=df["A1"].str.upper()

    if "BETA" not in df.columns and "OR" in df.columns:
        df["BETA"]=np.log10(df["OR"])


    
    f_file_name=abs_path_data+"/"+filename.split("/")[-1].split(".")[0]
    
    df[["MarkerName","A1","A2","N","BETA","P"]].to_csv(f_file_name+".tbl",sep="\t",index=False,header=True)
    ukbb_gzip_qced_cmd=f"""gzip -9 {f_file_name}.tbl"""
    subprocess.call(ukbb_gzip_qced_cmd, shell=True)

    LAVA_prep_out=f_file_name+".tbl.gz"
    return LAVA_prep_out


def Lava_input_file(trait_path:str, Enigma_input,lava_control,lava_cases):
    import pandas as pd
    import sys
    sys.path.append('/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/')
    import CONSTANTS 

    trait_name=trait_path.split("/")[-1].split(".")[0]
    trait_loc=CONSTANTS.Extra_temp_files_dict["extra_LAVA_filePrep"]+"/"+trait_name+".tbl.gz"

    
    
    LAVA_input_out_arr=[]#wsa_file+CONSTANTS.file_joiner_str+wTHICK_file
    for Enigma_wsa_wthick in ["wSA","wTHICK"]:
        input_list=CONSTANTS.Lava_cortical_subcortical[Enigma_input][Enigma_wsa_wthick]
        dataframe_data=pd.DataFrame([{"phenotype":trait_name.split(CONSTANTS.sep_btw_study_trait)[1],"cases":lava_cases,"controls":lava_control,"filename":trait_loc}])
        
        for line in input_list:
            line=line.rstrip("\n")
            line=line.rstrip("\r")
            Enigma_trait=line.split("ENIGMA3_mixed_se_")[1].split("_")[2]
            data_i={"phenotype":Enigma_trait,"cases":"NA","controls":"NA","filename":line}
            new_row_df = pd.DataFrame([data_i])
            dataframe_data = pd.concat([dataframe_data, new_row_df], ignore_index=True)
        df=pd.DataFrame(data=dataframe_data)
        f_name=CONSTANTS.Extra_temp_files_dict["extra_LAVA_input_files"]+"/"+Enigma_wsa_wthick+"_"+trait_name+".txt"
        df.to_csv(f_name, na_rep="NA",sep="\t",index=False)
        LAVA_input_out_arr.append(f_name)


    LAVA_input_out=CONSTANTS.file_joiner_str.join(LAVA_input_out_arr)
    print(LAVA_input_out)
    return LAVA_input_out

def LAVA_Matrix_Formation(rG_log_files:str):
        import re
        import sys
        
        sys.path.append('/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/')
        import CONSTANTS 


        print(rG_log_files)
        print(type(rG_log_files))
        rG_log_files_list=rG_log_files.split(CONSTANTS.file_joiner_str)
        

        trait_name=CONSTANTS.file_name_process(rG_log_files_list[0]).split(CONSTANTS.sep_btw_study_trait)[1][1:]

        print(trait_name)
        LAVA_matrix_file_out_arr=[]
        for Enigma_wsa_wthick in ["wSA","wTHICK"]:

            matrix=[[0]*(len(rG_log_files_list)+1) for _ in range(len(rG_log_files_list)+1)]
            matrix_id={}#trait_parcellation->index
            trait_parcellation_array=[]
            i=0
            
            
            for line in rG_log_files_list:
                line=line.strip("\n")
                line=line.strip("\r")
                

                Enigma_parcellation_arr=line.split("ENIGMA3_mixed_se_")[1].split("_")
                Enigma_parcellation=Enigma_parcellation_arr[2]

                if Enigma_parcellation_arr[0]!=Enigma_wsa_wthick:
                     continue
                trait_parcellation=[trait_name,Enigma_parcellation]

                print(trait_parcellation,i)



                with open(line, "r") as file:
                            content = file.read()
                        
                # Extract information from "Heritability of phenotype 1"
                heritability1 = re.search(r"Heritability of phenotype 1\n-+\nTotal Observed scale h2: (.*?) \((.*?)\)\nLambda GC: (.*?)\nMean Chi\^2: (.*?)\nIntercept: (.*?) \((.*?)\)\n", content)
                
                if heritability1:
                    total_h2_1 = heritability1.group(1)
                    total_h2_1_se = heritability1.group(2)
                    lambda_gc_1 = heritability1.group(3)
                    mean_chi2_1 = heritability1.group(4)
                    intercept_h1 = heritability1.group(5)
                    intercept_1_se = heritability1.group(6)



                heritability2 = re.search(r"Heritability of phenotype 2/2\n-+\nTotal Observed scale h2: (.*?) \((.*?)\)\nLambda GC: (.*?)\nMean Chi\^2: (.*?)\nIntercept: (.*?) \((.*?)\)\n", content)
                
                if heritability2:
                    total_h2_2 = heritability2.group(1)
                    total_h2_2_in_bracket = heritability2.group(2)
                    lambda_gc_2 = heritability2.group(3)
                    mean_chi2_2 = heritability2.group(4)
                    intercept_h2 = heritability2.group(5)
                    intercept_2_se = heritability2.group(6)

                genetic_covariance = re.search(r"Genetic Covariance\n-+\nTotal Observed scale gencov: (.*?) \((.*?)\)\nMean z1\*z2: (.*?)\nIntercept: (.*?) \((.*?)\)", content)
                if genetic_covariance:
                    total_gencov = genetic_covariance.group(1)
                    total_gencov_in_bracket = genetic_covariance.group(2)
                    mean_z1z2 = genetic_covariance.group(3)
                    intercept_gencov = genetic_covariance.group(4)
                    intercept_gencov_in_bracket = genetic_covariance.group(5)
                
                hypothesis_1_intercept=float(intercept_h1)
                hypothesis_2_intercept=float(intercept_h2)
                covariance_intercept=float(intercept_gencov)
                
                if trait_parcellation[0] not in matrix_id:
                    matrix_id[trait_parcellation[0]]=i
                    trait_parcellation_array.append(trait_parcellation[0])
                    i+=1
                matrix_id[trait_parcellation[1]]=i
                i+=1
                trait_parcellation_array.append(trait_parcellation[1])
                matrix[matrix_id[trait_parcellation[0]]][matrix_id[trait_parcellation[1]]]=covariance_intercept
                matrix[matrix_id[trait_parcellation[1]]][matrix_id[trait_parcellation[0]]]=covariance_intercept
                matrix[matrix_id[trait_parcellation[0]]][matrix_id[trait_parcellation[0]]]=hypothesis_1_intercept
                matrix[matrix_id[trait_parcellation[1]]][matrix_id[trait_parcellation[1]]]=hypothesis_2_intercept



            swapped_matrix_id = {value: key for key, value in matrix_id.items()} #index -> trait_parcellation
            csv_file = CONSTANTS.Extra_temp_files_dict["extra_LAVA_Matrix_files"]+"/"+Enigma_wsa_wthick+"_matrix_LAVA_Input_"+trait_name+".txt"
            LAVA_matrix_file_out_arr.append(csv_file)

            arr_index_list=[el_i for el_i in range(len(matrix_id.keys()))]
            print("Matrix for loop done",arr_index_list)
            all_parcellation_names=[swapped_matrix_id[el_i] for el_i in arr_index_list]

            with open(csv_file, 'w') as file:
                # Write column indices in the first row
                file.write("" + " ".join(all_parcellation_names) + "\n")
                
                # Write rows with row index and data
                for i, key2 in enumerate(all_parcellation_names):
                    file.write(key2 + " " + " ".join(map(str, matrix[i])) + "\n")

        print("Table saved to 'table.txt'")
        LAVA_matrix_file_out=CONSTANTS.file_joiner_str.join(LAVA_matrix_file_out_arr)
        return LAVA_matrix_file_out
        

def LAVA_shell_call_script(lava_matrix:str,lava_case_control:str) -> None:
    import os
    import sys
        
    sys.path.append('/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/')
    import CONSTANTS
    abs_LAVA=CONSTANTS.Extra_temp_files_dict["extra_LAVA_Shell_files"]+"/"

    lava_matrices=lava_matrix.split(CONSTANTS.file_joiner_str) ##for wSA and wTHICK
    lava_case_control_files=lava_case_control.split(CONSTANTS.file_joiner_str) ##for wSA and wTHICK
    

    trait_name=lava_case_control_files[0].split(CONSTANTS.sep_btw_study_trait)[1].split(".")[0].lstrip("_")
    

    shell_content='''#!/bin/bash
#$ -cwd
#$ -o /ifs/loni/faculty/njahansh/GAMBIT/genome_WebApplication/static/LAVA/logs/ -j y
#$ -N {LAVA}
#$ -t 1:18380
#$ -q iniadmin7.q


export PATH="/ifs/loni/faculty/njahansh/nerds/shruti/software/python/anaconda3/envs/r_env/bin:$PATH"

LOCI=($(seq 1 1 18380))
locus=${{LOCI[${{SGE_TASK_ID}}-1]}}

cd /ifs/loni/faculty/njahansh/nerds/ankush/for_Ravi/Pain_Genetics/Pain_LAVA_NonDuplicates/Pain_Lava_Results/


Rscript /ifs/loni/faculty/njahansh/GAMBIT/genome_WebApplication/static/LAVA/LAVA.R "/ifs/loni/faculty/njahansh/nerds/ravi/genetics/lava-mdd-endo-2023/data/g1000_eur" "/ifs/loni/faculty/njahansh/nerds/ravi/genetics/lava-mdd-endo-2023/data/analysis_files/gencode.v26.GRCh38.protein_coding.1Mb-cis.SNP-IDs.dbsnp-g1000-subset.loci" "{lava_case_control}" "{lava_matrix}" "{trait_name} bankssts caudalanteriorcingulate caudalmiddlefrontal cuneus entorhinal frontalpole fusiform inferiorparietal inferiortemporal insula isthmuscingulate lateraloccipital lateralorbitofrontal lingual medialorbitofrontal middletemporal paracentral parahippocampal parsopercularis parsorbitalis parstriangularis pericalcarine postcentral posteriorcingulate precentral precuneus rostralanteriorcingulate rostralmiddlefrontal superiorfrontal superiorparietal superiortemporal supramarginal temporalpole transversetemporal" "LAVA_{trait_name}_{wSA_wTHICK}" ${{locus}} "{trait_name}"
'''
    
    
    for lava_mat_i, case_control_i in zip(lava_matrices,lava_case_control_files):
        wSA_wTHICK=case_control_i.split("/")[-1].split("_")[0]
        formatted_content = shell_content.format(trait_name=trait_name,
                                                lava_case_control=case_control_i,
                                                lava_matrix=lava_mat_i,
                                                wSA_wTHICK=wSA_wTHICK,
                                                LAVA=trait_name)
    # Specify the file name for the shell script
        shell_script_file = abs_LAVA+wSA_wTHICK+"_"+trait_name+"_LAVA_script.sh"

        # Write the formatted content to the shell script file
        with open(shell_script_file, "w") as file:
            file.write(formatted_content)
        os.chmod(shell_script_file, 0o777)


