import subprocess

import sys
sys.path.append('/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/')

from gini_main import * ##file_checker
import CONSTANTS
def General_Munge(filePathInp:str,kwargs):
    import subprocess
    import pandas as pd
    import sys
    sys.path.append('/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/')
    import CONSTANTS 
        
    fileName=""
    filePathInp=filePathInp.replace("\n","")
    filePathInp=filePathInp.replace("\r","")
    
    df=pd.read_csv(filePathInp,sep="\t",nrows=15)
    
    fileName=filePathInp.split("/")[-1]#.split(".")[0]
    output_file_name=CONSTANTS.file_name_process(fileName)
    new_el=""
    if kwargs.pop("OR",False):
        new_el+=f"""--signed-sumstats {kwargs.pop("signed_sumstats","OR")},1 \
            """
    else:
        new_el+=f"""--signed-sumstats {kwargs.pop("signed_sumstats","BETA")},0 \
            """
    # with open('/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/association_format.json', 'r') as json_file:
    #     json_association = json.load(json_file)

    # gwas_format_name=filePathInp.split("/")[-1] if filePathInp.split("/")[-1].split(".")[-1].lower()!="gz" else filePathInp.split("/")[-1].split(".")[-2]
    # dict_col_header_mapper=json_association[gwas_format_name]['format_dict']

    
    # output_file_name=fileName.replace(".tbl","").replace(".csv","").replace(".txt","").replace(".gz","").replace(".tsv","").replace(".vcf","")
    
    out_file_name=kwargs.pop('out',CONSTANTS.Extra_temp_files_dict["extras_LDSC_Munge_files"])+"/"+output_file_name
    new_el+=f"""--sumstats {filePathInp} \
        --out {out_file_name} \
        --a1 {kwargs.pop('a1',"A1")} \
        --a2 {kwargs.pop('a2',"A2")} \
            """

    # if "signed_sumstats" in kwargs:
    #     new_el+=f"""--signed-sumstats {kwargs.pop("signed_sumstats","OR")},{kwargs.pop("signed_sumstats_int","0")} \
    #         """

    if "N_col" in kwargs:
        new_el+=f"""--N-col {kwargs.pop("N_col","N")} \
            """
            
    for k,v in kwargs.items():
        k=k.replace("_","-")
        new_el+=f"""--{k} {v} \
            """



    cmd_line=f"/ifs/loni/faculty/njahansh/nerds/ankush/software/Anaconda/envs/ldsc/bin/python2.7 \
            /ifs/loni/faculty/njahansh/nerds/ankush/webApplication_Genome/Git_Genome/GenomeAPI/ldsc/munge_sumstats_Meta.py  \
        --chunksize 500000 \
        --merge-alleles /ifs/loni/faculty/njahansh/nerds/ravi/genetics/ldsc/eur_w_ld_chr/w_hm3.snplist \
        "+new_el


    print(cmd_line)
    subprocess.call(cmd_line, shell=True)
    
    LDSC_Munge_out=out_file_name+".sumstats.gz"
    
    return LDSC_Munge_out


def HeritabilityLDSC(filePathInp_Munged):
    import subprocess

    
    import sys
    sys.path.append('/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/')
    import CONSTANTS


    kwargs_Munged=dict(out=CONSTANTS.Extra_temp_files_dict["extras_LDSC_Heri_files"]+"/",
                       ref_ld_chr="/ifs/loni/faculty/njahansh/nerds/ravi/genetics/ldsc/1000G_EUR_Phase3_baseline/baseline.",
                       frqfile_chr="/ifs/loni/faculty/njahansh/nerds/ravi/genetics/ldsc/1000G_Phase3_frq/1000G.EUR.QC.",
                       w_ld_chr="/ifs/loni/faculty/njahansh/nerds/ravi/genetics/ldsc/weights_hm3_no_hla/weights.")
    
    filePathInp=filePathInp_Munged
    kwargs=kwargs_Munged

    fileNamesList=filePathInp.split(" ")
    all_files=",".join(fileNamesList)

    if len(fileNamesList)>1:
        fileOutName="MultiFileOut"
    else:
        fileOutName="Heritability_"+fileNamesList[0].split("/")[-1].split(".gz")[0]
    

    new_el=""
    out_heritable_file=kwargs.pop('out',CONSTANTS.Extra_temp_files_dict["extras_LDSC_Heri_files"]+"/")+fileOutName

    new_el=f"""--out {out_heritable_file} \
        --ref-ld-chr {kwargs.pop('ref_ld_chr',"/ifs/loni/faculty/njahansh/nerds/ravi/genetics/ldsc/1000G_EUR_Phase3_baseline/baseline.")} \
            --frqfile-chr {kwargs.pop('frqfile_chr',"/ifs/loni/faculty/njahansh/nerds/ravi/genetics/ldsc/1000G_Phase3_frq/1000G.EUR.QC.")} \
                --w-ld-chr {kwargs.pop('w_ld_chr',"/ifs/loni/faculty/njahansh/nerds/ravi/genetics/ldsc/weights_hm3_no_hla/weights.")} \
"""
            
    for k,v in kwargs.items():
        k=k.replace("_","-")
        new_el+=f"""--{k} {v} \
            """


    cmd_line=f"/ifs/loni/faculty/njahansh/nerds/ankush/software/Anaconda/envs/ldsc/bin/python2.7 \
        /ifs/loni/faculty/njahansh/nerds/ankush/webApplication_Genome/Git_Genome/GenomeAPI/ldsc/ldsc.py \
            --h2 "+all_files+" \
            --overlap-annot \
            "+new_el
    
    subprocess.call(cmd_line, shell=True)
    return out_heritable_file



def Heritability_Log_Extraction(input_file_arr:list):
    import pandas as pd
    import re
    
    import sys
    sys.path.append('/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/')
    import CONSTANTS

    # input_file="/ifs/loni/faculty/njahansh/nerds/ankush/webApplication_Genome/Git_Genome/GenomeAPI_2/File_Path_text_Files/metaAnalysisHeritabilityLog.txt"
    df_final=pd.DataFrame(columns = ["Trait Name","h\u00b2", "h\u00b2 SE", "Mean Chi\u00b2","LambdaGC", "Intercept","Intercept SE","Ratio"])
    out_file_name="_".join(input_file_arr[0].split("/")[-1].split(CONSTANTS.sep_btw_study_trait)[0].split("_")[1:])+"_"+str(len(input_file_arr))
    # with open(input_file,"r") as lines:
    for line in input_file_arr:
        line=line.strip("\r")
        line=line.strip("\n")
        trait_name="_".join(line.split("/")[-1].split(CONSTANTS.sep_btw_study_trait)[1].split("_")[1:]).strip(".log").strip("sumstats")


        with open(line, "r") as file:
            content = file.read()
        
        # Extract information from "Heritability of phenotype 1"
        # heritability1 = re.search(r"Heritability of phenotype 1\n-+\nTotal Observed scale h2: ([\d.-]+) \([\d.-]+\)\nLambda GC: ([\d.]+)\nMean Chi\^2: ([\d.]+)\nIntercept: ([\d.]+) \([\d.]+\)\nRatio: ([\d.-]+) \([\d.-]+\)", content)
        # heritability1 = re.search(r"Heritability of phenotype 1\n-+\nTotal Observed scale h2: (.*?) \((.*?)\)\nLambda GC: (.*?)\nMean Chi\^2: (.*?)\nIntercept: (.*?) \((.*?)\)\nRatio: (.*?) \((.*?)\)", content)
        heritability1 = re.search(r"Total Observed scale h2: (.*?) \((.*?)\)", content)
        
        heritability2 = re.search(r"Lambda GC: (.*?)\nMean Chi\^2: (.*?)\nIntercept: (.*?) \((.*?)\)\nRatio: (.*?) \((.*?)\)", content)
        
        
        
        her_1_dict={"Trait Name":trait_name}

        if heritability1:
            total_h2_1 = heritability1.group(1)
            total_h2_1_se = heritability1.group(2)
            her_1_dict["h\u00b2"]=total_h2_1
            her_1_dict["h\u00b2 SE"]=total_h2_1_se
        if heritability2:
            lambda_gc_1 = heritability2.group(1)
            mean_chi2_1 = heritability2.group(2)
            intercept_1 = heritability2.group(3)
            intercept_1_se = heritability2.group(4)
            ratio_1 = heritability2.group(5)
            ratio_1_se = heritability2.group(6)
            her_1_dict["Mean Chi\u00b2"]=mean_chi2_1
            her_1_dict["LambdaGC"]=lambda_gc_1
            her_1_dict["Intercept"]=intercept_1
            her_1_dict["Intercept SE"]=intercept_1_se
            her_1_dict["Ratio"]=ratio_1

        new_row_df = pd.DataFrame([her_1_dict])

        df_final = pd.concat([df_final, new_row_df], ignore_index=True)


    df_final.to_csv(CONSTANTS.Extra_temp_files_dict["extras_LDSC_Heri__results"]+f"/{out_file_name}traits.csv",index=False)


def CellTypeLDSC(files_Munged:str,CellAnalysisName:str,LDCTSFile:str):
    
    import pandas as pd
    import re
    import subprocess
    import sys
    sys.path.append('/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/')
    import CONSTANTS

    fileNamesList=files_Munged.split(" ")
    all_files=",".join(fileNamesList)



    fileOutName=CellAnalysisName.lower()+"_"+fileNamesList[0].split("/")[-1].split(".")[0]    

    cmd_line_n="/ifs/loni/faculty/njahansh/nerds/ankush/software/Anaconda/envs/ldsc/bin/python2.7 \
    /ifs/loni/faculty/njahansh/nerds/ankush/webApplication_Genome/Git_Genome/GenomeAPI/ldsc/ldsc.py \
    --h2-cts "+all_files+" \
    --ref-ld-chr /ifs/loni/faculty/njahansh/nerds/ravi/genetics/ldsc/1000G_EUR_Phase3_baseline/baseline. \
    --out "+CONSTANTS.Extra_temp_files_dict["extra_cell_h2"]+"/"+fileOutName+" \
    --ref-ld-chr-cts "+CONSTANTS.cell_type_ldsc_ref+LDCTSFile+" \
    --w-ld-chr /ifs/loni/faculty/njahansh/nerds/ravi/genetics/ldsc/weights_hm3_no_hla/weights."


    subprocess.call(cmd_line_n, shell=True)
    return CONSTANTS.Extra_temp_files_dict["extra_cell_h2"]+fileOutName



def trait_Combinations_for_rG(trait_input_txt1:list, trait_input_txt2:list):#For creating combination between two files (replace the 2nd txt with Enigma)
    import pandas as pd
    import re
    
    import sys
    sys.path.append('/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/')
    import CONSTANTS

    # df1=pd.read_csv(trait_input_txt1,header=None)
    df1 = pd.DataFrame({'0': trait_input_txt1})

    # df2=pd.read_csv(trait_input_txt2,header=None)
    df2 = pd.DataFrame({'0': trait_input_txt2})


    df_12_joinCombn=df1.merge(df2, how='cross')
    df_12_joinCombn["combined"]=df_12_joinCombn["0_x"]+","+df_12_joinCombn["0_y"]

    return df_12_joinCombn["combined"]

    # new_comb_series=df_12_joinCombn["combined"]
    # new_comb_series.to_csv(abs_path+"Combined_Meta_ukbb_abcd_Neuro.txt",index=False)

def rG_LDSC(trait_files:str):

    import subprocess
    import sys
    sys.path.append('/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/')
    import CONSTANTS

    trait_files=trait_files.replace("\n","").replace("\r","")

    trait_file_1=trait_files.split(",")[0] #
    trait_file_2=trait_files.split(",")[1] #
    

    fileOutName=trait_file_1.split("/")[-1].split(".")[0]+"___"+trait_file_2.split("/")[-1].split(".")[0]
        
    out_file_comp=CONSTANTS.Extra_temp_files_dict["extras_LDSC_rG_results"]+"/"+ fileOutName

    cmd_line_n="/ifs/loni/faculty/njahansh/nerds/ankush/software/Anaconda/envs/ldscENV_3_14_2024/bin/python \
        /ifs/loni/faculty/njahansh/nerds/ankush/webApplication_Genome/Git_Genome/GenomeAPI_2/ldsc/ldsc.py  \
        --out "+out_file_comp +" \
        --rg "+ trait_files +"\
        --ref-ld-chr /ifs/loni/faculty/njahansh/nerds/ravi/genetics/ldsc/eur_w_ld_chr/ \
        --w-ld-chr /ifs/loni/faculty/njahansh/nerds/ravi/genetics/ldsc/eur_w_ld_chr/ "

    print(cmd_line_n)
    subprocess.call(cmd_line_n, shell=True)

    return out_file_comp


def rG_Log_Extraction(filePath:list):
    
    import pandas as pd
    import sys
    sys.path.append('/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/')
    import CONSTANTS
    
    df_final=pd.DataFrame(columns = ['p1', 'p2', 'rg','se', 'z','p'])
    

    for f_name in filePath:
        f_name=f_name.strip("\n").strip("\r")
        with open(f_name) as w:
            f_list=[]
            dict_temp={}

            bool_val_p_z=False

            for line in w:
                if line.strip(" ").startswith("p1"):
                    line2=line.replace(" ","\t").split("\t")
                    res = list(filter(lambda item: item != "", line2))
                    res[-1]=res[-1].rstrip("\n")
                    f_list.append(res)
                    bool_val_p_z=True

                elif bool_val_p_z:
                    line2=line.replace(" ","\t").split("\t")
                    res = list(filter(lambda item: item != "", line2))
                    res[-1]=res[-1].rstrip("\n")
                    f_list.append(res)
                    bool_val_p_z=False

            if len(f_list)==2:
                for u,v in zip(f_list[0],f_list[1]):
                    dict_temp[u]=v
                dict_temp["Log_File_result"]=f_name

                new_row_df = pd.DataFrame([dict_temp])

                df_final = pd.concat([df_final, new_row_df], ignore_index=True)

    # df_final[['p1', 'p2', 'rg','se', 'z','p','Log_File_result']].to_csv(CONSTANTS.Extra_temp_files_dict["rG_CSV_files"]+"/Global_rG_"+file_name_out+".csv",index=False)
    df_final[['p1', 'p2', 'rg','se', 'z','p','Log_File_result']].to_csv("/ifs/loni/faculty/njahansh/nerds/ankush/for_Ravi/Pain_Genetics/Pain_rG/Global_rG_COPC.csv",index=False)

