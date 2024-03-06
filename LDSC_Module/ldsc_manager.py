import subprocess


def General_Munge(filePathInp:str,kwargs):
    import subprocess
    extras_LDSC_files="/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/Exta_temp_files/LDSC_Out/Munged_results/"
    # kwargs=dict(N_col="TotalN",snp="MarkerName",frq="Freq1",signed_sumstats="EffectARE",out=extras_LDSC_files,a1="Allele1",a2="Allele2",p="PvalueARE")
        
    fileName=""
    filePathInp=filePathInp.replace("\n","")
    filePathInp=filePathInp.replace("\r","")
    
    fileName=filePathInp.split("/")[-1]#.split(".")[0]
    
    output_file_name=fileName.replace(".tbl","").replace(".csv","").replace(".txt","").replace(".gz","").replace(".tsv","").replace(".vcf","")
    
    new_el=""
    out_file_name=kwargs.pop('out',"/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/Exta_temp_files/LDSC_Out/Munged_results/")+output_file_name
    new_el+=f"""--sumstats {filePathInp} \
        --out {out_file_name} \
        --a1 {kwargs.pop('a1',"A1")} \
        --a2 {kwargs.pop('a2',"A2")} \
            """

    if "signed_sumstats" in kwargs:
        new_el+=f"""--signed-sumstats {kwargs.pop("signed_sumstats","OR")},{kwargs.pop("signed_sumstats_int","0")} \
            """

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

    kwargs_Munged=dict(out="/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/Exta_temp_files/LDSC_Out/Heritability_Results/",
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
    out_heritable_file=kwargs.pop('out',"/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/Exta_temp_files/LDSC_Out/Heritability_Results/")+fileOutName

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
    # input_file="/ifs/loni/faculty/njahansh/nerds/ankush/webApplication_Genome/Git_Genome/GenomeAPI_2/File_Path_text_Files/metaAnalysisHeritabilityLog.txt"
    df_final=pd.DataFrame(columns = ["Trait Name","h\u00b2", "h\u00b2 SE", "Mean Chi\u00b2","LambdaGC", "Intercept","Intercept SE","Ratio"])

    # with open(input_file,"r") as lines:
    for line in input_file_arr:
        line=line.strip("\r")
        line=line.strip("\n")
        trait_name="_".join(line.split("/")[-1].split("1RE")[0].split("_")[1:])

        if not (trait_name.split("_")[0] in ["Total","Witelson5"] and trait_name.split("_")[-1] in ["Area","MeanThickness"]):
            continue

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

        print("Heritability of phenotype 1:")
        print("Total Observed scale h2:", total_h2_1)
        print("Total Observed scale h2 in bracket:", total_h2_1_se)
        print("Lambda GC:", lambda_gc_1)
        print("Mean Chi^2:", mean_chi2_1)
        print("Intercept:", intercept_1)
        print("Intercept in bracket:", intercept_1_se)
        print("Ratio:", ratio_1)
        print("Ratio in bracket:", ratio_1_se)

        df_final=df_final.append(her_1_dict,ignore_index = True)

    df_final.to_csv("/ifs/loni/faculty/njahansh/GAMBIT/genome_WebApplication/static/Tables_LDSC/MetaAnalysis_Heritability_Log_12Traits.csv",index=False)

