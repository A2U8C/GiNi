# GiNi-GWAS post-GWAS process: 
## Python-based, post-GWAS processing and analysis pipeline 

GiNi is a automated pipeline built for post-GWAS analysis. With GWAS summary statistics as input, the pipeline processes and perform a set of comprehensive analyses. Uses brain gene-expression and splicing data resources for analysis. GiNi includes the following modules:
* [Global Genetic Correlation](https://github.com/bulik/ldsc)
* [Local Genetic Correlation](https://github.com/josefin-werme/LAVA)
* [Mendelian Randomization](https://github.com/jean997/cause)
* [Heritability](https://github.com/bulik/ldsc)
* Transcriptome-Wide Associations



## Overview of the pipeline
<img width="800" alt="workflow" src="https://github.com/A2U8C/GiNi_post_GWAS_processing/assets/69053051/735647e4-02b6-44e6-ace3-e0e395de959a">
</p>

## How to use the package:
Clone the github directory using:
```bash
https://github.com/A2U8C/GiNi_post_GWAS_processing.git
```
 
## Virtual environment:
Navigate to the "smacc" folder and then create a virtual environment using the requirements.txt file:
```bash
conda create -n GiNi_Env python== 3.8 -y
conda activate GiNi_Env
pip install -r requirements.txt
pip install .
```

## Prepare underlying tools:
```bash
## Run the following python script, to create environments and other tools which will be required by GiNi
python setup_tools_Gini.py
wget -O output_filename.ext https://www.dropbox.com/s/GiNi_underlying_data?dl=1
```

## Input file fomat:
GiNi expect `input_txt` argument which is path to the input file which contains all GWAS files path. Example of `input_txt` input
```
/path/to/GWAS/stats/studyName__TraitName
```

```
/path/to/GWAS/stats/ukbb_regenie/ukbb__Total_MeanThickness.regenie
/path/to/GWAS/stats/ukbb_regenie/ukbb__Witelson5_Genu_Area.regenie
/path/to/GWAS/stats/ukbb_regenie/ukbb__Witelson5_Isthmus_MeanThickness.regenie
/path/to/GWAS/stats/abcd_regenie/abcd__Total_MeanThickness.regenie
/path/to/GWAS/stats/abcd_regenie/abcd__Witelson5_Genu_Area.regenie
/path/to/GWAS/stats/abcd_regenie/abcd__Witelson5_Isthmus_MeanThickness.regenie
```

## Test the tool:
```bash
python gini_main.py input-module-wrapper \
	--input_txt /path/to/GWAS/stats/ABCD_UKBB_input.txt \
	--n_studies 1 \
	--ethnicity European \
	--analysis_list Heritability \
	--tissues_cells Astrocytes \
	--meta random \
	--gwas_format regenie \
	--lava_control_cases /path/to/Case_Controls/COPC_Case_control.txt
```

--input_txt : To the input file

--n_studies : Number of studies, to verify if Meta-Analysis is required

--ethnicity : Ethnicity for selecting the appropriate reference 

--analysis_list : Comma seperated analysis list or 'ALL' for all analysis

--tissues_cells : Brain tissue and cells to performe TWAS and heritability

--meta : Analysis type for meta-analysis

--gwas_format : Format of the input GWAS

--control_cases : File containing count of control/cases for each trait

--computing_options : If HPC is available, else tests will run sequentially

## If you use this code, please cite the following:
* Willer, C. J., et al. METAL: fast and efficient meta-analysis of genomewide association scans. Bioinformatics 26, 2190–2191 (2010).
* Bulik-Sullivan, B. et al. LD score regression distinguishes confounding from polygenicity in genome-wide association studies. Nat. Genet. 47, 291–295 (2015).
* Werme, J., et al.  An integrated framework for local genetic correlation analysis. Nat. Genet. 54, 274–282 (2022).
* Morrison J, et al. Mendelian randomization accounting for correlated and uncorrelated pleiotropic effects using genome-wide summary statistics. Nat Genet. 2020.
* Finucane, H. K. et al. Heritability enrichment of specifically expressed genes identifies disease-relevant tissues and cell types. Nat. Genet. 50, 621–629 (2018).
* de Leeuw, C., et al. On the interpretation of transcriptome-wide association studies. PLoS Genet. 19, e1010921 (2023).

