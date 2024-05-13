

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
        dataframe_data=pd.DataFrame([{"phenotype":trait_name.split(CONSTANTS.sep_btw_study_trait)[1].lstrip("_"),"cases":lava_cases,"controls":lava_control,"filename":trait_loc}])
        
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



def Lava_input_file_Enigma(trait_path:str, Enigma_input:str,Enigma_wsa_wthick:str,lava_control:str,lava_cases:str):
    import pandas as pd
    import sys
    sys.path.append('/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/')
    import CONSTANTS 

    trait_name=trait_path.split("/")[-1].split(".")[0]
    trait_loc=CONSTANTS.Extra_temp_files_dict["extra_LAVA_filePrep"]+"/"+trait_name+".tbl.gz"

    
    
    input_list=CONSTANTS.Lava_cortical_subcortical[Enigma_input][Enigma_wsa_wthick]
    dataframe_data=pd.DataFrame([{"phenotype":trait_name.split(CONSTANTS.sep_btw_study_trait)[1].lstrip("_"),"cases":lava_cases,"controls":lava_control,"filename":trait_loc}])
    
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


    LAVA_input_out=f_name
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
        

        trait_name=CONSTANTS.file_name_process(rG_log_files_list[0]).split(CONSTANTS.sep_btw_study_trait)[1].lstrip("_")

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
        

def LAVA_Matrix_Formation_Enigma(rG_log_files:str,Enigma_wsa_wthick:str):
    import re
    import sys
    
    sys.path.append('/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/')
    import CONSTANTS 


    print(rG_log_files)
    print(type(rG_log_files))
    rG_log_files_list=rG_log_files.split(CONSTANTS.file_joiner_str)
    

    trait_name=CONSTANTS.file_name_process(rG_log_files_list[0]).split(CONSTANTS.sep_btw_study_trait)[1].lstrip("_")

    print(len(rG_log_files_list))

    matrix=[[0]*(35) for _ in range(35)]
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
    LAVA_matrix_file_out=csv_file
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


def LAVA_shell_call_script_Enigma(lava_matrix:str,lava_case_control:str) -> None:
    import os
    import sys
        
    sys.path.append('/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/')
    import CONSTANTS
    abs_LAVA=CONSTANTS.Extra_temp_files_dict["extra_LAVA_Shell_files"]+"/"

    

    trait_name=lava_case_control.split(CONSTANTS.sep_btw_study_trait)[1].split(".")[0].lstrip("_")
    

    shell_content='''#!/bin/bash
#$ -cwd
#$ -o /ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/Exta_temp_files/Logs/ -j y
#$ -N {LAVA}
#$ -t 1:10
#$ -q iniadmin7.q


LOCI=($(seq 1 1 18380))
locus=${{LOCI[${{SGE_TASK_ID}}-1]}}

cd {output_folder}


Rscript /ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/LAVA_Module/LAVA.R "/ifs/loni/faculty/njahansh/nerds/ravi/genetics/lava-mdd-endo-2023/data/g1000_eur" "/ifs/loni/faculty/njahansh/nerds/ravi/genetics/lava-mdd-endo-2023/data/analysis_files/gencode.v26.GRCh38.protein_coding.1Mb-cis.SNP-IDs.dbsnp-g1000-subset.loci" "{lava_case_control}" "{lava_matrix}" "{trait_name} bankssts caudalanteriorcingulate caudalmiddlefrontal cuneus entorhinal frontalpole fusiform inferiorparietal inferiortemporal insula isthmuscingulate lateraloccipital lateralorbitofrontal lingual medialorbitofrontal middletemporal paracentral parahippocampal parsopercularis parsorbitalis parstriangularis pericalcarine postcentral posteriorcingulate precentral precuneus rostralanteriorcingulate rostralmiddlefrontal superiorfrontal superiorparietal superiortemporal supramarginal temporalpole transversetemporal" "LAVA_{trait_name}_{wSA_wTHICK}" ${{locus}} "{trait_name}" "{output_folder}"
'''
    

    wSA_wTHICK=lava_case_control.split("/")[-1].split("_")[0]
    output_folder=CONSTANTS.Extra_temp_files_dict["extra_LAVA_output"]+"/"+trait_name+"/"+wSA_wTHICK
    formatted_content = shell_content.format(trait_name=trait_name,
                                            lava_case_control=lava_case_control,
                                            lava_matrix=lava_matrix,
                                            wSA_wTHICK=wSA_wTHICK,
                                            output_folder=output_folder,
                                            LAVA=trait_name)
# Specify the file name for the shell script
    shell_script_file = abs_LAVA+wSA_wTHICK+"_"+trait_name+"_LAVA_script.sh"

    # Write the formatted content to the shell script file
    with open(shell_script_file, "w") as file:
        file.write(formatted_content)
    os.chmod(shell_script_file, 0o777)





def LAVA_run_function(lava_matrix:str,lava_case_control:str,locus):

    import os
    import sys
    import subprocess
        
    sys.path.append('/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/')
    import CONSTANTS
     
    trait_name=lava_case_control.split(CONSTANTS.sep_btw_study_trait)[1].split(".")[0].lstrip("_")
    # trait_name=lava_case_control.split("___")[1].split(".")[0].lstrip("_")
    
    wSA_wTHICK=lava_case_control.split("/")[-1].split("_")[0]
    output_folder=CONSTANTS.Extra_temp_files_dict["extra_LAVA_output"]+"/"+trait_name+"/"+wSA_wTHICK
    # output_folder="/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/Exta_temp_files/LAVA_Out/LAVA_Results"+"/"+trait_name+"/"+wSA_wTHICK

    r_command=f'Rscript /ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/LAVA_Module/LAVA.R "/ifs/loni/faculty/njahansh/nerds/ravi/genetics/lava-mdd-endo-2023/data/g1000_eur" "/ifs/loni/faculty/njahansh/nerds/ravi/genetics/lava-mdd-endo-2023/data/analysis_files/gencode.v26.GRCh38.protein_coding.1Mb-cis.SNP-IDs.dbsnp-g1000-subset.loci" "{lava_case_control}" "{lava_matrix}" "{trait_name} bankssts caudalanteriorcingulate caudalmiddlefrontal cuneus entorhinal frontalpole fusiform inferiorparietal inferiortemporal insula isthmuscingulate lateraloccipital lateralorbitofrontal lingual medialorbitofrontal middletemporal paracentral parahippocampal parsopercularis parsorbitalis parstriangularis pericalcarine postcentral posteriorcingulate precentral precuneus rostralanteriorcingulate rostralmiddlefrontal superiorfrontal superiorparietal superiortemporal supramarginal temporalpole transversetemporal" "LAVA_{trait_name}_{wSA_wTHICK}" {locus} "{trait_name}" {output_folder}'


    print(r_command)

#     r_command = [
#     "Rscript",
#     "/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/LAVA_Module/LAVA.R",
#     "/ifs/loni/faculty/njahansh/nerds/ravi/genetics/lava-mdd-endo-2023/data/g1000_eur",
#     "/ifs/loni/faculty/njahansh/nerds/ravi/genetics/lava-mdd-endo-2023/data/analysis_files/gencode.v26.GRCh38.protein_coding.1Mb-cis.SNP-IDs.dbsnp-g1000-subset.loci",
#     lava_case_control,
#     lava_matrix,
#     trait_name+" bankssts caudalanteriorcingulate caudalmiddlefrontal cuneus entorhinal frontalpole fusiform inferiorparietal inferiortemporal insula isthmuscingulate lateraloccipital lateralorbitofrontal lingual medialorbitofrontal middletemporal paracentral parahippocampal parsopercularis parsorbitalis parstriangularis pericalcarine postcentral posteriorcingulate precentral precuneus rostralanteriorcingulate rostralmiddlefrontal superiorfrontal superiorparietal superiortemporal supramarginal temporalpole transversetemporal",
#     "LAVA_"+trait_name+"_"+wSA_wTHICK,
#     locus,
#     trait_name,
#     output_folder
# ]
    
    # Define the environment variables specific to your R environment
    r_environment_path = "/ifs/loni/faculty/njahansh/nerds/shruti/software/python/anaconda3/envs/r_env/bin"

    r_environment = os.environ.copy()
    r_environment["PATH"] = f"{r_environment_path}:{r_environment['PATH']}"

    # Execute the command
    subprocess.run(r_command, env=r_environment,shell=True)

    return output_folder





def LAVA_run_function_split_locus(lava_matrix:str,lava_case_control:str,locii):

    import os
    import sys
    import subprocess
        
    sys.path.append('/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/')
    import CONSTANTS
     
    trait_name=lava_case_control.split(CONSTANTS.sep_btw_study_trait)[1].split(".")[0].lstrip("_")
    # trait_name=lava_case_control.split("___")[1].split(".")[0].lstrip("_")
    
    wSA_wTHICK=lava_case_control.split("/")[-1].split("_")[0]
    output_folder=CONSTANTS.Extra_temp_files_dict["extra_LAVA_output"]+"/"+trait_name+"/"+wSA_wTHICK
    # output_folder="/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/Exta_temp_files/LAVA_Out/LAVA_Results"+"/"+trait_name+"/"+wSA_wTHICK

    locii_list=locii.split(CONSTANTS.locii_seperator)
    r_environment_path = "/ifs/loni/faculty/njahansh/nerds/shruti/software/python/anaconda3/envs/r_env/bin"

    r_environment = os.environ.copy()
    r_environment["PATH"] = f"{r_environment_path}:{r_environment['PATH']}"
    for locus in locii_list:
        r_command=f'Rscript /ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/LAVA_Module/LAVA.R "/ifs/loni/faculty/njahansh/nerds/ravi/genetics/lava-mdd-endo-2023/data/g1000_eur" "/ifs/loni/faculty/njahansh/nerds/ravi/genetics/lava-mdd-endo-2023/data/analysis_files/gencode.v26.GRCh38.protein_coding.1Mb-cis.SNP-IDs.dbsnp-g1000-subset.loci" "{lava_case_control}" "{lava_matrix}" "{trait_name} bankssts caudalanteriorcingulate caudalmiddlefrontal cuneus entorhinal frontalpole fusiform inferiorparietal inferiortemporal insula isthmuscingulate lateraloccipital lateralorbitofrontal lingual medialorbitofrontal middletemporal paracentral parahippocampal parsopercularis parsorbitalis parstriangularis pericalcarine postcentral posteriorcingulate precentral precuneus rostralanteriorcingulate rostralmiddlefrontal superiorfrontal superiorparietal superiortemporal supramarginal temporalpole transversetemporal" "LAVA_{trait_name}_{wSA_wTHICK}" {int(locus)} "{trait_name}" {output_folder}'


        print(r_command)
        subprocess.run(r_command, env=r_environment,shell=True)

    return output_folder



def LAVA_shell_call_script_Enigma_Nipype(lava_matrix:str,lava_case_control:str,trait_LAVA_prepd:str,machineName:str="runnow.q"):
    import os
    import sys
    import math
    sys.path.append('/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/')
    import CONSTANTS

    import nipype
    import nipype.pipeline.engine as pe
    import nipype.interfaces.utility as niu
    from nipype.interfaces.utility import Function
    import shutil
    from LAVA_Module import LAVA_script
    from nipype import config

    trait_name=lava_case_control.split(CONSTANTS.sep_btw_study_trait)[1].split(".")[0].lstrip("_")

    wSA_wTHICK=lava_case_control.split("/")[-1].split("_")[0]
    
    wf_name='LAVA_exec_'+wSA_wTHICK+"_"+trait_name
    workflow_metal_mungng = pe.Workflow(name=wf_name)
    workflow_metal_mungng.base_dir = CONSTANTS.Extra_temp_files_dict["extra_Nipype_Wf"]+"/" #"/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/Exta_temp_files/"
    workflow_metal_mungng.config["job_finished_timeout"]=CONSTANTS.timeout_wait_setter

    
    config_dict={'execution': {'job_finished_timeout': CONSTANTS.timeout_wait_setter}}
    config.update_config(config_dict)
    # config.enable_resource_monitor()
    # output_folder="/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/Exta_temp_files/LAVA_Out/LAVA_Results"+"/"+trait_name+"/"+wSA_wTHICK
    output_folder=CONSTANTS.Extra_temp_files_dict["extra_LAVA_output"]+"/"+trait_name+"/"+wSA_wTHICK
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
        print("Nested folders created successfully.")

    num_of_sub_jobs=400
    num_of_locii_in_each_job=math.ceil(18380.0/num_of_sub_jobs)
    # locus_comp_list=list(range(1,18381))
    locus_comp_list = [str(i) for i in range(1, 18381)]
    iterable_new_list = [CONSTANTS.locii_seperator.join(locus_comp_list[i:i+num_of_locii_in_each_job]) for i in range(0, len(locus_comp_list), num_of_locii_in_each_job)]

    # iterable_new_list=list(range(1,11))
    print(len(iterable_new_list))
    inputnode = pe.Node(
            niu.IdentityInterface(
                fields=['locus_i']
                ),
            synchronize=True,
            iterables=[('locus_i',iterable_new_list)],
            name='inputnode'
        )
        
    node1_parallel_LAVA = pe.Node(Function(input_names=['lava_matrix','lava_case_control','locii'],
                    output_names=['out_file_comp'],
                    function=LAVA_script.LAVA_run_function_split_locus),
            name='node1_parallel_LAVA')
    node1_parallel_LAVA.inputs.lava_matrix=lava_matrix
    node1_parallel_LAVA.inputs.lava_case_control=lava_case_control

    
    connections=[(inputnode, node1_parallel_LAVA, [('locus_i', 'locii')])
                 ] 
        

    workflow_metal_mungng.connect(connections)    
    res=workflow_metal_mungng.run('SGE',plugin_args=dict(
        max_jobs =num_of_sub_jobs,
        overwrite= True, 
        qsub_args=f'-q {machineName} -o /ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/Exta_temp_files/Logs/ -j y'))

    # res=workflow_metal_mungng.run('SGE',plugin_args=dict(dont_resubmit_completed_jobs= True, overwrite= True, qsub_args=f'-q {machineName} -o /ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/Exta_temp_files/Logs/ -j y'))
    
    outputs = CONSTANTS.file_joiner_str.join([x.result.outputs.out_file_comp for x in res.nodes() if x.name == 'node1_parallel_LAVA'])
    return outputs






def LAVA_shell_call_script_Nipype_without_split(lava_matrix:str,lava_case_control:str,trait_LAVA_prepd:str,machineName:str="runnow.q"):
    import os
    import sys
        
    sys.path.append('/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/')
    import CONSTANTS

    import nipype
    import nipype.pipeline.engine as pe
    import nipype.interfaces.utility as niu
    from nipype.interfaces.utility import Function
    import shutil
    from LAVA_Module import LAVA_script
    from nipype import config

    trait_name=lava_case_control.split(CONSTANTS.sep_btw_study_trait)[1].split(".")[0].lstrip("_")

    wSA_wTHICK=lava_case_control.split("/")[-1].split("_")[0]
    
    wf_name='LAVA_exec_'+wSA_wTHICK+"_"+trait_name
    workflow_metal_mungng = pe.Workflow(name=wf_name)
    workflow_metal_mungng.base_dir = CONSTANTS.Extra_temp_files_dict["extra_Nipype_Wf"]+"/" #"/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/Exta_temp_files/"
    workflow_metal_mungng.config["job_finished_timeout"]=CONSTANTS.timeout_wait_setter

    
    config_dict={'execution': {'job_finished_timeout': CONSTANTS.timeout_wait_setter}}
    config.update_config(config_dict)
    # config.enable_resource_monitor()
    # output_folder="/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/Exta_temp_files/LAVA_Out/LAVA_Results"+"/"+trait_name+"/"+wSA_wTHICK
    output_folder=CONSTANTS.Extra_temp_files_dict["extra_LAVA_output"]+"/"+trait_name+"/"+wSA_wTHICK
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
        print("Nested folders created successfully.")

    iterable_new_list=list(range(1,18381))
    # iterable_new_list=list(range(1,11))
    print(len(iterable_new_list))
    inputnode = pe.Node(
            niu.IdentityInterface(
                fields=['locus_i']
                ),
            synchronize=True,
            iterables=[('locus_i',iterable_new_list)],
            name='inputnode'
        )
        
    node1_parallel_LAVA = pe.Node(Function(input_names=['lava_matrix','lava_case_control','locus'],
                    output_names=['out_file_comp'],
                    function=LAVA_script.LAVA_run_function),
            name='node1_parallel_LAVA')
    node1_parallel_LAVA.inputs.lava_matrix=lava_matrix
    node1_parallel_LAVA.inputs.lava_case_control=lava_case_control

    
    connections=[(inputnode, node1_parallel_LAVA, [('locus_i', 'locus')])
                 ] 
        

    workflow_metal_mungng.connect(connections)    
    res=workflow_metal_mungng.run('SGE',plugin_args=dict(
        max_jobs =400,
        overwrite= True, 
        qsub_args=f'-q {machineName} -o /ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/Exta_temp_files/Logs/ -j y'))

    # res=workflow_metal_mungng.run('SGE',plugin_args=dict(dont_resubmit_completed_jobs= True, overwrite= True, qsub_args=f'-q {machineName} -o /ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/Exta_temp_files/Logs/ -j y'))
    
    outputs = CONSTANTS.file_joiner_str.join([x.result.outputs.out_file_comp for x in res.nodes() if x.name == 'node1_parallel_LAVA'])
    return outputs



