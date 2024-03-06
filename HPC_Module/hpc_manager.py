
import sys
import os

sys.path.append('/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/')
# print(os.getcwd())
# os. chdir('/ifs/loni/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/')
print(os.getcwd())


from gini_main import * ##file_checker
from METAL_Module.metal_script import * #metal_execution_function




argp=argparse.ArgumentParser()
argp.add_argument("--trait_Name_loc")
argp.add_argument("--analysis_name")
argp.add_argument("--meta_type")
args=argp.parse_args()

if args.analysis_name=="Meta_Analysis":
    trait_path_str=args.trait_Name_loc
    trait_i=trait_path_str.split(":")[0]
    trait_i_studies_path=trait_path_str.split(":")[1].split(",")
    metal_execution_function(trait_i,trait_i_studies_path,args.meta_type)