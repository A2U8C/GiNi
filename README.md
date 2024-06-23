# GiNi-GWAS post-GWAS process: 
## Segment, Measure and AutoQC the midsagittal Corpus Callosum

GiNi is a automated pipeline built for post-GWAS analysis. With GWAS summary statistics as input, the pipeline processes and perform a set of comprehensive analyses. GiNi includes the following modules:
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

## To install underlying tools:
```bash
## Run the following python script, to create environments and other tools which will be required by GiNi
python setup_tools_Gini.py
```

## Input Preprocessing:
All the MR images should be registred to MNI 1mm template(182 X 218 X 182) with 6dof. You can use the template provided in the "model" folder on github. You can use the FSL's flirt command for linear registration:
```bash
flirt -in ${inpdir}/${subj}.nii.gz \
	-ref ${MNI_1mm_template} \
  	-out ${outdir}/${subj} \
 	-dof 6 \
  	-cost mutualinfo \
  	-omat ${outdir}/matrices/${subj}_MNI_6p.xfm
```

## Test the tool:
```bash
smacc -f ./subject_list.txt -o ./smacc_output -m t1
```
-f : Text file with a list of absolute paths to the niftis to be processed and names to save the outputs for each subject. Check example text file "subject_list.txt" provided. <br />
-o : Absolute path of output folder <br />
-m : Modality of the images to be processed (t1/t2/flair) <br />
-q : Optional flag to perform Automated QC on the segmentations <br />
The final output is a csv which will contain all the extracted shape metrics and a column "QC label" indicating whether the segmentations were accurate(0)/fail(1) if the QC flag is provided.


## If you use this code, please cite the following:
* Willer, C. J., et al. METAL: fast and efficient meta-analysis of genomewide association scans. Bioinformatics 26, 2190–2191 (2010).
* Bulik-Sullivan, B. et al. LD score regression distinguishes confounding from polygenicity in genome-wide association studies. Nat. Genet. 47, 291–295 (2015).
* Werme, J., et al.  An integrated framework for local genetic correlation analysis. Nat. Genet. 54, 274–282 (2022).
* Morrison J, et al. Mendelian randomization accounting for correlated and uncorrelated pleiotropic effects using genome-wide summary statistics. Nat Genet. 2020.
* Finucane, H. K. et al. Heritability enrichment of specifically expressed genes identifies disease-relevant tissues and cell types. Nat. Genet. 50, 621–629 (2018).
* de Leeuw, C., et al. On the interpretation of transcriptome-wide association studies. PLoS Genet. 19, e1010921 (2023).

