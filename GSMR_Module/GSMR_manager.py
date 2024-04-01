

import pandas as pd
import numpy as np
from scipy import stats
import math


def GSMRPREP(all_paths:list):
    for line in all_paths:
        print(line)
        line=line.strip("\n")
        line=line.strip("\r")
        df=pd.read_csv(line,delim_whitespace=True)
        
        # print(df)
        df=df.rename({"MarkerName":"SNP","variant_id":"SNP",
                      "Allele2":"A2","other_allele":"A2",
                      "Allele1":"A1","effect_allele":"A1",
                      "PvalueARE":"p","PVAL":"p","PvalueARE":"p","p_value":"p","P":"p",
                      "Freq1":"freq","FRQ":"freq","effect_allele_freq":"freq","effect_allele_frequency":"freq",
                      "EffectARE":"b","BETA":"b","beta":"b",
                      "TotalN":"N",
                      "StdErrARE":"se","SE":"se","standard_error":"se"
                      },axis=1)
        
        print(df.columns)
        df["A2"]=df["A2"].str.upper()
        df["A1"]=df["A1"].str.upper()

        if "b" not in df.columns and "OR" in df.columns:
            df["b"]=np.log10(df["OR"])

        if "N_CASE" in df.columns and "N_CONTROL" in df.columns:
            df["N"]=df["N_CASE"]+df["N_CONTROL"]
        
        
        print(df.columns)
        if "b" in df.columns and "p" in df.columns and "se" not in df.columns:
            #z-score=Effect/StdErr
            #dd$SE_fixed_qnorm<-abs(log(dd$OR)/qnorm(dd$P/2))
            df['se'] = np.abs(np.log10(df['OR']) / np.abs(stats.norm.ppf(df['p'] / 2)))

            # df["se"]=df["b"]/df["Z"]


        name_line=line.split("/")[-1].split(".")[0]

        coll_selector=["SNP","A1","A2","freq","b","se","p","N"]
        # if "b" in df.columns:
        #     coll_selector=["SNP","A1","A2","freq","b","se","p","N"]
        # df[coll_selector].to_csv("/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/Exta_temp_files/GSMR_Out/GSMR_temp/"+line.split("/")[-1].split(".")[0]+".txt",sep="\t",index=None)
        df[coll_selector].to_csv("/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/Exta_temp_files/GSMR_Out/GSMR_nonDuplicates_input/"+line.split("/")[-1].split(".")[0]+".txt",sep="\t",index=None)
class MyDict(dict):
    def __missing__(self, key):
        return key

class FilePrep_CSV():
    def function_Name(self,file_name):
        #df_enigma_files["Enigma_name_2"].str.split("_", expand=True)[5]+("_SA" if df_enigma_files["Enigma_name_2"].str.split("_", expand=True)[6]=="surfavg" else "_CT") + " " + df_enigma_files["Enigma_name"]
        return file_name.split("_")[5]+("_SA" if file_name.split("_")[6]=="surfavg" else "_CT")

    def __init__(self,mainFile:str):
        filePath=mainFile.strip("\n")
        filePath=filePath.strip("\r")

        #/ifs/loni/faculty/njahansh/GAMBIT/genome_WebApplication/static/MungedData/test_cc_UKBB/test_cc_ukbb_white_Witelson5_AnteriorBody_Area.sumstats.gz


        main_folder_name=filePath.split("/")[-1].split("_")[0]

        df_ukbb_name_dict=pd.read_csv("/ifs/loni/faculty/njahansh/nerds/ankush/webApplication_Genome/Git_Genome/GenomeAPI_2/File_Path_text_Files/metaAnalysisGSMRPrepFile.txt",header=None,names=["file_path_ukbb"])
        df_ukbb_name_dict["key_ukbb_file"]=df_ukbb_name_dict["file_path_ukbb"].str.split("/",expand=True)[10].str.split(".",expand=True)[0]
        dict_ukk_fname=MyDict(df_ukbb_name_dict.set_index('key_ukbb_file')['file_path_ukbb'].to_dict())

        if filePath.rstrip(".csv").split("_")[-2]=="GC":
            df_enigmaGC_name_dict=pd.read_csv("/ifs/loni/faculty/njahansh/nerds/ankush/webApplication_Genome/Git_Genome/GenomeAPI_2/File_Path_text_Files/enigmaGC_gsmr_file.txt",header=None,names=["file_path_enigma"])
            df_enigmaGC_name_dict["key_enigma_file"]=df_enigmaGC_name_dict["file_path_enigma"].str.split("/",expand=True)[10].str.split(".",expand=True)[0]
            dict_enigmaGC_fname=MyDict(df_enigmaGC_name_dict.set_index('key_enigma_file')['file_path_enigma'].to_dict())
        else:
            df_enigmaGC_name_dict=pd.read_csv("/ifs/loni/faculty/njahansh/nerds/ankush/webApplication_Genome/Git_Genome/GenomeAPI_2/File_Path_text_Files/enigmaNoGC_gsmr_file.txt",header=None,names=["file_path_enigma"])
            df_enigmaGC_name_dict["key_enigma_file"]=df_enigmaGC_name_dict["file_path_enigma"].str.split("/",expand=True)[10].str.split(".",expand=True)[0]
            dict_enigmaGC_fname=MyDict(df_enigmaGC_name_dict.set_index('key_enigma_file')['file_path_enigma'].to_dict())


        if filePath!="":
            df=pd.read_csv(filePath,usecols=["p1Full","p2Full"])
            df["p1Full"]=df["p1Full"].str.split("/", expand=True)[10].str.split(".sumstats", expand=True)[0]
            df["p2Full"]=df["p2Full"].str.split("/", expand=True)[10].str.split(".sumstats", expand=True)[0]

            df2=df.groupby('p1Full').apply(lambda x: list(x.p2Full)).reset_index(name="Enigm_files_name").to_dict(orient="records")
            
            print(df2)
            # for el in df2:
            #     # file_name=el["p1Full"].lstrip("test_cc_ukbb_white_")+ ("_Enigma_GC" if filePath.rstrip(".csv").split("_")[-1]=="GC" else "_Enigma_NoGC")
            #     #JHU3_Area_Enigma_NoGC.csv
            #     print(el["p1Full"])
            #     file_name=el["p1Full"].split("1RE")[0]+ ("_Enigma_GC" if filePath.rstrip(".csv").split("_")[-2]=="GC" else "_Enigma_NoGC")+"_"+filePath.rstrip(".csv").split("_")[-1]
            #     enigma_files=el["Enigm_files_name"]
            #     df_enigma_files=pd.DataFrame(data={"Enigma_name":enigma_files,"Enigma_name_2":enigma_files})
            #     df_enigma_files["Enigma_name"]=df_enigma_files["Enigma_name"].map(dict_enigmaGC_fname   )
            #     # df_enigma_files["Final_name_list"]=df_enigma_files["Enigma_name_2"].str.split("_", expand=True)[5]+("_SA" if df_enigma_files["Enigma_name_2"].str.split("_", expand=True)[6]=="surfavg" else "_CT") + " " + df_enigma_files["Enigma_name"]
            #     df_enigma_files["Final_name_list"]=df_enigma_files["Enigma_name_2"].apply(self.function_Name) + " " + df_enigma_files["Enigma_name"]
            #     df_enigma_files["Final_name_list"].to_csv(r'/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/Exta_temp_files/GSMR_Out/GSMR_input/'+file_name+"__brain.txt",index=None,header=None)
                
            #     # with open('/ifs/loni/faculty/njahansh/GAMBIT/genome_WebApplication/static/GSMR_Preprocess/MetaAnalysis_GSMR_InputFiles/'+main_folder_name+'/'+file_name+"__brain.txt", 'w') as fw:
            #     #     fw.write(df_enigma_files["Final_name_list"])

                
            #     with open('/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/Exta_temp_files/GSMR_Out/GSMR_input/'+file_name+'__meta.txt', 'w') as f:
            #         f.write(el["p1Full"].split("1RE")[0]+" "+dict_ukk_fname[el["p1Full"]])

class Pain_FilePrep_CSV():
    def function_Name(self,file_name):
        #df_enigma_files["Enigma_name_2"].str.split("_", expand=True)[5]+("_SA" if df_enigma_files["Enigma_name_2"].str.split("_", expand=True)[6]=="surfavg" else "_CT") + " " + df_enigma_files["Enigma_name"]
        return file_name.split("_")[5]+("_SA" if file_name.split("_")[6]=="surfavg" else "_CT")

    def __init__(self,mainFile:str):
        filePath=mainFile.strip("\n")
        filePath=filePath.strip("\r")

        #/ifs/loni/faculty/njahansh/GAMBIT/genome_WebApplication/static/MungedData/test_cc_UKBB/test_cc_ukbb_white_Witelson5_AnteriorBody_Area.sumstats.gz


        main_folder_name="COPC"

        
        df_pain_name_dict=pd.read_csv("/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/TEMP_inputs/pain_gsmr_inp.txt",header=None,names=["key_enigma_file","file_path_enigma"])
        # df_enigmaGC_name_dict["key_enigma_file"]=df_enigmaGC_name_dict["file_path_enigma"].str.split("/",expand=True)[10].str.split(".",expand=True)[0]
        dict_pain_fname=MyDict(df_pain_name_dict.set_index('key_enigma_file')['file_path_enigma'].to_dict())

        dict_p2_handler={"Suicide Attempt":"Suicide_Attempt",
                         "Neuroticism":"Neuroticism",
                         "Insomnia":"Insomnia",
                         "ADHD":"ADHD",
                         "Intelligence":"Intelligence",
                         "ASD":"ASD",
                         "IBS":"IBS",
                         "Anxiety":"Anxiety",
                         "COPC":"COPC"}

        if filePath!="":
            df=pd.read_csv(filePath,usecols=["p1","p2","p"])
            df=df[df["p"]<0.00625]
            df["p1Full"]=df["p1"]
            df["p1Full"]= df["p1Full"].str.split("/").str[-1].str.split("_").str[0] 
            df["p2Full"]=df["p2"].map(dict_p2_handler)

            df2=df.groupby('p1Full').apply(lambda x: list(x.p2Full)).reset_index(name="Enigm_files_name").to_dict(orient="records")
            
            print(df2)
            for el in df2:
                file_name=el["p1Full"]
                enigma_files=el["Enigm_files_name"]
                df_enigma_files=pd.DataFrame(data={"Enigma_name":enigma_files,"Enigma_name_2":enigma_files})
                df_enigma_files["Enigma_name"]=df_enigma_files["Enigma_name"].map(dict_pain_fname)
                # df_enigma_files["Final_name_list"]=df_enigma_files["Enigma_name_2"].str.split("_", expand=True)[5]+("_SA" if df_enigma_files["Enigma_name_2"].str.split("_", expand=True)[6]=="surfavg" else "_CT") + " " + df_enigma_files["Enigma_name"]
                df_enigma_files["Final_name_list"]=df_enigma_files["Enigma_name_2"] + " " + df_enigma_files["Enigma_name"]
                df_enigma_files["Final_name_list"].to_csv(r'/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/Exta_temp_files/GSMR_Out/GSMR_input/'+file_name+"__brain.txt",index=None,header=None)
                
                #     fw.write(df_enigma_files["Final_name_list"])

                
                with open('/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/Exta_temp_files/GSMR_Out/GSMR_input/'+file_name+'__meta.txt', 'w') as f:
                    f.write(el["p1Full"]+" "+dict_pain_fname[el["p1Full"]])



def func_fixer(dict_path:str):
    #rs2189970
    df_root=pd.read_csv(dict_path,delim_whitespace=True)#usecols=["CHR","BP","SNP","A1","A2","N","P","new_freq","OR"]
    df_uk=pd.read_csv("/ifs/loni/faculty/njahansh/datasets/UKBB/dataset/processed/genetics/imputed_plink_files/plink.frq",delim_whitespace=True,usecols=["CHR","SNP","MAF"])
    glob_df=pd.read_csv("/ifs/loni/faculty/njahansh/nerds/ravi/genetics/LAVA-main/support_data/eur/plink.frq",delim_whitespace=True,usecols=["CHR","SNP","MAF"])
    glob_df=glob_df.rename({"MAF":"freq_1k"},axis=1)
    # # # print(glob_df.columns)
    df_root=df_root.rename({"MarkerName":"SNP","variant_id":"SNP",
                      "Allele2":"A2","other_allele":"A2",
                      "Allele1":"A1","effect_allele":"A1",
                      "PvalueARE":"p","PVAL":"p","PvalueARE":"p","p_value":"p","P":"p",
                      "Freq1":"freq","FRQ":"freq","effect_allele_freq":"freq","effect_allele_frequency":"freq","new_freq":"freq",
                      "EffectARE":"b","BETA":"b","beta":"b",
                      "TotalN":"N",
                      "StdErrARE":"se","SE":"se","standard_error":"se"
                      },axis=1)
    

    # print(df_root[df_root['SNP']=="rs4553305"])
    if "b" not in df_root.columns and "OR" in df_root.columns:
        df_root["b"]=np.log10(df_root["OR"])

    if "b" in df_root.columns and "p" in df_root.columns and "se" not in df_root.columns:
    #     #z-score=Effect/StdErr
    #     #dd$SE_fixed_qnorm<-abs(log(dd$OR)/qnorm(dd$P/2))
        print("entered")
        df_root['se'] = np.abs(df_root["b"] / np.abs(stats.norm.ppf(df_root['p'] / 2)))
    coll_selector=["SNP","A1","A2","freq","b","se","p","N"]
    
    # # df_root[coll_selector].to_csv("/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/Exta_temp_files/GSMR_Out/GSMR_temp/"+dict_path.split("/")[-1].split("-")[0]+".txt",sep="\t",index=None)

    df_merged=pd.merge(df_root,df_uk,how="left",on=["CHR","SNP"])
    df_merged=pd.merge(df_merged,glob_df,how="left",on=["CHR","SNP"])


    # # # # #rs12380191
    def funct_MAF_noNa(row):
        # if pd.isna(row['MAF']):
        #     return row["freq_1k"]
        # return row['MAF']
    
    
        if pd.isna(row['freq_1k']):
            return row["MAF"]
        return row['freq_1k']
    
    df_merged["freq"]=df_merged.apply(funct_MAF_noNa, axis=1)
    # # # print(df_merged.dtypes)
    # df_merged=df_merged[~df_merged["freq"].isna()]
    # df_merged[coll_selector].to_csv("/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/Exta_temp_files/GSMR_Out/GSMR_nonDuplicates_input/Anxiety_03_21_2024.txt",index=False)#, sep="\t"
    
    df_merged_non_Dups=df_merged[df_merged.duplicated(subset=['SNP'], keep=False)]
    # df_merged=pd.read_csv("/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/Exta_temp_files/GSMR_Out/GSMR_nonDuplicates_input/Anxiety_noNA.txt",sep="\t")
    df_merged_non_Dups=df_merged_non_Dups[~df_merged_non_Dups["se"].isna()]
    # print(df2[df2["se"].isna()])
    df_merged_non_Dups.to_csv("/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/Exta_temp_files/GSMR_Out/GSMR_nonDuplicates_input/Anxiety_noNA_se_freq.txt",sep="\t",index=False)
        


#  
    # df[df.duplicated(subset=['SNP'], keep=False)].sort_values(by="p")


def rem_duplicates(): 
    # Neuroticism,/ifs/loni/faculty/njahansh/nerds/ravi/genetics/neurotic_gwas/neuroticism.ma
    # Insomnia,/ifs/loni/faculty/njahansh/nerds/ravi/genetics/insomnia/insomnia.ma
    # ADHD,/ifs/loni/faculty/njahansh/nerds/ravi/genetics/adhd/adhd.ma
    # Suicide_Attempt,/ifs/loni/faculty/njahansh/nerds/ravi/genetics/suicide_gwas/suicide.ma
    # Intelligence,/ifs/loni/faculty/njahansh/nerds/ravi/genetics/intelligence/intelligence.ma
    # ASD,/ifs/loni/faculty/njahansh/nerds/ravi/genetics/asd/asd.ma
    # IBS,/ifs/loni/faculty/njahansh/nerds/ravi/pain_gwases/ibs/GCST90016564_buildGRCh37.txt
    # Anxiety,/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/Exta_temp_files/GSMR_Out/GSMR_temp/Anxiety.txt
    # COPC,/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/Exta_temp_files/GSMR_Out/GSMR_temp/COPC_Z_QC_noMHC.txt

    list_path=[
        # "/ifs/loni/faculty/njahansh/nerds/ravi/genetics/neurotic_gwas/neuroticism.ma",
        #        "/ifs/loni/faculty/njahansh/nerds/ravi/genetics/insomnia/insomnia.ma",
        #        "/ifs/loni/faculty/njahansh/nerds/ravi/genetics/adhd/adhd.ma",
        #        "/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/Exta_temp_files/GSMR_Out/GSMR_temp/suicide.ma",
        #        "/ifs/loni/faculty/njahansh/nerds/ravi/genetics/intelligence/intelligence.ma",
        #        "/ifs/loni/faculty/njahansh/nerds/ravi/genetics/asd/asd.ma",
        #        "/ifs/loni/faculty/njahansh/nerds/ravi/pain_gwases/ibs/GCST90016564_buildGRCh37.txt",
        #        "/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/Exta_temp_files/GSMR_Out/GSMR_temp/Anxiety.txt",
               "/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/Exta_temp_files/GSMR_Out/GSMR_temp/COPC_Z_QC_noMHC.txt"
               ]
    
    for list_i in list_path:
        df=pd.read_csv(list_i,sep="\t")
        df=df[~df.duplicated(subset=['SNP'], keep=False)]#.sort_values(by="p")
        
        df.to_csv("/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/Exta_temp_files/GSMR_Out/GSMR_nonDuplicates_input/"+list_i.split("/")[-1], sep="\t", index=False)

# rem_duplicates()

class GSMR_MendelianRandomisation():
    def __init__(self,brainFile:str):
        import subprocess
        #/ifs/loni/faculty/njahansh/GAMBIT/genome_WebApplication/static/GSMR_Preprocess/z_temp_file/JHU3/JHU3_Body_Area_Enigma_GC__brain.txt
         
        fileUkbbPath=brainFile.split("__")[0]+"__meta.txt"

        out_file_name=brainFile.split("/")[-1].split("__")[0]

        folder_format_name=brainFile.split("/")[-1].split("_")[0]
        
        # print()
        # if brainFile.split("/")[-3]=="MetaAnalysis_Neuro_GSMR_InputFiles":
        #     out_main_folder_static="MetaAnalysis_Neuro_GSMR"
        # elif brainFile.split("/")[-3]=="MetaAnalysis_GSMR_InputFiles":  
        #     out_main_folder_static= "MetaAnalysis_Enigma_GSMR"

        # cmd_line="/ifs/loni/faculty/njahansh/nerds/iyad/software/gcta/gcta-1.94.1 \
        #     --bfile /ifs/loni/faculty/njahansh/nerds/shruti/for_neda/grant_things/cc/genetics/gwas/ukb.qc5 \
        #     --gsmr-file "+brainFile+" "+fileUkbbPath+" \
        #     --gsmr-direction 2 \
        #     --out /ifs/loni/faculty/njahansh/nerds/ankush/for_Ravi/Pain_Genetics/Pain_GSMR/"+out_file_name+" \
        #     --effect-plot \
        #     --gsmr2-beta \
        #     --gwas-thresh 5e-8\
        #     --threads 5"
        
        cmd_line="/ifs/loni/faculty/njahansh/nerds/iyad/software/gcta/gcta-1.94.1 \
            --bfile /ifs/loni/faculty/njahansh/nerds/ravi/genetics/LAVA-main/support_data/eur/g1000_eur \
            --gsmr-file "+brainFile+" "+fileUkbbPath+" \
            --gsmr-direction 2 \
            --out /ifs/loni/faculty/njahansh/nerds/ankush/for_Ravi/Pain_Genetics/Pain_GSMR/"+out_file_name+" \
            --effect-plot \
            --gsmr2-beta \
            --gwas-thresh 5e-8\
            --diff-freq 1\
            --threads 5"
        
        subprocess.call(cmd_line, shell=True)

        print(cmd_line)

# /ifs/loni/faculty/njahansh/nerds/iyad/software/gcta/gcta-1.94.1 \
#           --bfile /ifs/loni/faculty/njahansh/nerds/ravi/genetics/LAVA-main/support_data/eur/g1000_eur \
#                 --ld /ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/Exta_temp_files/GSMR_Out/GSMR_nonDuplicates_input/neuroticism_only_sigs.ma  \
#                     --ld-wind 5000  \
#                         --ld-sig 0.05  \
#                             --out /ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/Exta_temp_files/GSMR_Out/GSMR_LD_SNPS/neuroticism_ld_SNPS

GSMR_MendelianRandomisation('/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/Exta_temp_files/GSMR_Out/GSMR_input/COPC__brain.txt')

# Anxiety /ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/Exta_temp_files/GSMR_Out/GSMR_nonDuplicates_input/Anxiety_noNA_se_freq.txt

# GSMRPREP(["/ifs/loni/faculty/njahansh/nerds/ravi/genetics/li_anxiety_2023/Anxiety-Disorders-plink.meta.P_Sorted.2021.3.31_Final.txt"])
# FilePrep_CSV("/ifs/loni/faculty/njahansh/GAMBIT/genome_WebApplication/static/GSMR_Preprocess/rG_MetaAnalysis_Enigma_SE_Significant/Hits_JHU_Total_Wit5/Wit5_MeanThickness_Enigma_GC_CT.csv")
# Pain_FilePrep_CSV("/ifs/loni/faculty/njahansh/nerds/ankush/for_Ravi/Pain_Genetics/Pain_rG/Global_rG_COPC.csv")

# func_fixer("/ifs/loni/faculty/njahansh/nerds/ravi/genetics/li_anxiety_2023/Anxiety-Disorders-plink.meta.P_Sorted.2021.3.31_Final.txt")
# func_fixer("/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/Exta_temp_files/GSMR_Out/GSMR_temp/Anxiety-Disorders-plink.txt")#"/ifs/loni/faculty/njahansh/nerds/ravi/genetics/li_anxiety_2023/Anxiety-Disorders-plink.meta.P_Sorted.2021.3.31_Final.txt")

# trait_file_1="/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/TEMP_inputs/Pain_traits_txt1.txt"
# trait_file_2="/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/TEMP_inputs/Pain_traits_txt2.txt"
# with open(trait_file_1, 'r') as file1:
#     trait_file_1_files = file1.read().splitlines()

# with open(trait_file_2, 'r') as file2:
#     trait_file_2_files = file2.read().splitlines()

# iterable_new_list=trait_file_1_files+trait_file_2_files

# GSMRPREP(iterable_new_list)