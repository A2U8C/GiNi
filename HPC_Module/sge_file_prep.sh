#!/bin/sh
#$ -S /bin/bash
#$ -o /ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/Exta_temp_files/Logs/ -j y 
#$ -wd /ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/



export PATH=/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/GiNi_Env/bin:$PATH

main_package_path="/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/"





if [ "$1" == "Meta_Analysis" ]; then
    dirS=${2}
    SUBJECTS=($(cat "${dirS}"))
    subject=${SUBJECTS[${SGE_TASK_ID}-1]}

    python3 "${main_package_path}HPC_Module/hpc_manager.py" --trait_Name_loc "${subject}" --analysis_name "${1}" --meta_type "${3}"

fi