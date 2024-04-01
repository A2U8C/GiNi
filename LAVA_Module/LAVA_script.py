import pandas as pd
import numpy as np
import re
import os

def Lava_Meta_FilePrep(filename:str):
    import subprocess
    import pandas as pd
    abs_path="/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/Exta_temp_files/LAVA_Out/LAVA_data/"
    df=pd.read_csv(filename,delim_whitespace=True)


    # df["A2"]=df["A2"].str.upper()
    # df["A1"]=df["A1"].str.upper()

    # if "b" not in df.columns and "OR" in df.columns:
    #     df["b"]=np.log10(df["OR"])


    disease_name_map={"GCST90016564_buildGRCh37": "ibs",}
    
    ukb_f_name=filename.split("/")[-1].split("1RE")[0]
    print(ukb_f_name)
    ukb_f_name=disease_name_map[ukb_f_name] if ukb_f_name in disease_name_map else ukb_f_name

    df.rename(columns={"SNP":"MarkerName",
                        "Allele1":"A1" , 
                        "Allele2":"A2" , 
                        "TotalN":"N" , 
                        "EffectARE":"BETA","b":"BETA", 
                        "PvalueARE":"P","p":"P"
                        },inplace=True)
    df[["MarkerName","A1","A2","N","BETA","P"]].to_csv(abs_path+ukb_f_name+".tbl",sep="\t",index=False,header=True)
    ukbb_gzip_qced_cmd=f"""gzip -9 {abs_path+ukb_f_name}.tbl"""
    subprocess.call(ukbb_gzip_qced_cmd, shell=True)


# Lava_Meta_FilePrep("/ifs/loni/faculty/njahansh/nerds/ravi/genetics/li_anxiety_2023/Anxiety-Disorders-plink.meta.P_Sorted.2021.3.31_Final.txt")

class Lava_input_file:
    def __init__(self, trait_name, input_list):

        trait_loc="/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/Exta_temp_files/LAVA_Out/LAVA_data/"+trait_name+".tbl.gz"

        cases_control={
            "COPC":{
                "cases":82812,
                "controls":81966 
            },
            "neuroticism":{
                "cases":None,
                "controls":None
            },
            "insomnia":{
                "cases":593724,
                "controls":1771286 
            },
            "adhd":{
                "cases":38691 ,
                "controls":186843
            },
            "suicide":{
                "cases":26590,
                "controls":492022
            },
            "intelligence":{
                "cases":None,
                "controls":None
            },
            "asd":{
                "cases":18381,
                "controls":27969
            },
            "ibs":{
                "cases":53400,
                "controls":433201
            },
            "anxiety":{
                "cases":74973,
                "controls":400243
            },
        }
        
        dataframe_data=pd.DataFrame([{"phenotype":trait_name,"cases":cases_control[trait_name]['cases'],"controls":cases_control[trait_name]['controls'],"filename":trait_loc}])
        
        
        for line_i in input_list:
            line="/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/Exta_temp_files/LAVA_Out/LAVA_data/"+line_i+".tbl.gz"
            line=line.rstrip("\n")
            line=line.rstrip("\r")
            disease_file=line.split("/")[-1].split(".")[0]
            disease_name=disease_file

            data_i={"phenotype":disease_name,"cases":cases_control[disease_name]['cases'],"controls":cases_control[disease_name]['controls'],"filename":line}

            new_row_df = pd.DataFrame([data_i])

            dataframe_data = pd.concat([dataframe_data, new_row_df], ignore_index=True)

        df=pd.DataFrame(data=dataframe_data)
        f_name="/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/Exta_temp_files/LAVA_Out/LAVA_input/"+trait_name+".txt"
        df.to_csv(f_name, na_rep="NA",sep="\t",index=False)
        print(f_name)

      
# Lava_input_file("COPC",["suicide","neuroticism","adhd","asd","insomnia","intelligence","ibs","anxiety"])

class LAVA_Matrix_Formation:
    def __init__(self) -> None:

        trait_name="COPC"
        file_out_name="/ifs/loni/faculty/njahansh/nerds/ankush/webApplication_Genome/Git_Genome/GenomeAPI_2/File_Path_text_Files/painCorrelationLogs.txt"
        

        matrix=[[0]*9 for _ in range(9)]
        matrix_id={}#trait_parcellation->index
        trait_parcellation_array=[]
        i=0
        
        with open(file_out_name) as f:
            for line in f:
                line=line.strip("\n")
                line=line.strip("\r")
                #Total_Area_pgc-panic2019.log
                if line.split("/")[-1].split(".log")[0].split("_")[-1]=="noGC":
                    continue
                
                disease_parcellation=line.split("___")[-1].split("_")[0].lower()

                trait_parcellation=[trait_name,disease_parcellation]




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
        csv_file = "/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/Exta_temp_files/LAVA_Out/LAVA_Matrix/matrix_LAVA_Input_"+trait_name+".txt"

        arr_index_list=[el_i for el_i in range(len(matrix_id.keys()))]
        all_parcellation_names=[swapped_matrix_id[el_i] for el_i in arr_index_list]

        print(matrix)

        with open(csv_file, 'w') as file:
            # Write column indices in the first row
            file.write("" + " ".join(all_parcellation_names) + "\n")
            
            # Write rows with row index and data
            for i, key2 in enumerate(all_parcellation_names):
                file.write(key2 + " " + " ".join(map(str, matrix[i])) + "\n")

        print("Table saved to 'table.txt'")
        





class LAVA_shell_call_script:
    def __init__(self) -> None:
        ## Neuro
        # file_12_traits="/ifs/loni/faculty/njahansh/nerds/ankush/webApplication_Genome/Git_Genome/GenomeAPI_2/File_Path_text_Files/12_trait_files_output.txt"
        abs_LAVA="/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/Exta_temp_files/LAVA_Out/LAVA_shell/"
        # lava_job_name=
        shell_content='''#!/bin/bash
#$ -cwd
#$ -o /ifs/loni/faculty/njahansh/GAMBIT/genome_WebApplication/static/LAVA/logs/ -j y
#$ -N {LAVA}
#$ -t 1:18380
#$ -q iniadmin7.q


export PATH="/ifs/loni/faculty/njahansh/nerds/shruti/software/python/anaconda3/envs/r_env/bin:$PATH"

LOCI=($(seq 1 1 18380))
locus=${{LOCI[${{SGE_TASK_ID}}-1]}}

cd /ifs/loni/faculty/njahansh/nerds/ankush/for_Ravi/Pain_Genetics/Pain_Lava/

Rscript /ifs/loni/faculty/njahansh/GAMBIT/genome_WebApplication/static/LAVA/LAVA.R "/ifs/loni/faculty/njahansh/nerds/ravi/genetics/lava-mdd-endo-2023/data/g1000_eur" "/ifs/loni/faculty/njahansh/nerds/ravi/genetics/lava-mdd-endo-2023/data/analysis_files/gencode.v26.GRCh38.protein_coding.1Mb-cis.SNP-IDs.dbsnp-g1000-subset.loci" "/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/Exta_temp_files/LAVA_Out/LAVA_input/{trait_name}.txt" "/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/Exta_temp_files/LAVA_Out/LAVA_Matrix/matrix_LAVA_Input_{trait_name}.txt" "{trait_name} suicide neuroticism adhd asd insomnia intelligence ibs anxiety" "LAVA_{trait_name}" ${{locus}} "{trait_name}"

'''
        
        line="COPC"
        formatted_content = shell_content.format(trait_name=line,LAVA=line)
        # Specify the file name for the shell script
        shell_script_file = abs_LAVA+line+"_LAVA_script.sh"

        # Write the formatted content to the shell script file
        with open(shell_script_file, "w") as file:
            file.write(formatted_content)
        os.chmod(shell_script_file, 0o777)


        # with open(file_12_traits) as f:
            # for line in f:
            #     line=line.strip("\n")
            #     line=line.strip("\r")
            #     formatted_content = shell_content.format(trait_name=line,LAVA=line)
            #     # Specify the file name for the shell script
            #     shell_script_file = abs_LAVA+line+"_LAVA_script.sh"

            #     # Write the formatted content to the shell script file
            #     with open(shell_script_file, "w") as file:
            #         file.write(formatted_content)
            #     os.chmod(shell_script_file, 0o777)


LAVA_shell_call_script()
# LAVA_Matrix_Formation()