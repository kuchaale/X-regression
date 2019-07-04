#!/bin/bash

output_folder=/mnt/2data/solarisheppa/outputs_new/
reg_folder=/mnt/2data/solarisheppa/regressors/bernd/
sim=$4 #refC2

s_year=$1 #1960
e_year=$2 #2010
variable=$3 #zmo3,zmnoy

model=SOCOL3
ens=r1i1p1
config_name=config_${model}${ens}_${variable}_${s_year}-${e_year}.yaml
variable=${variable} sim=${sim} model=${model} ens=${ens} s_year=${s_year} e_year=${e_year} output_folder=${output_folder} reg_folder=${reg_folder} ./templater.sh config_template.yaml > ${config_name}
echo "running "${config_name}
python3 solaris_heppa_analysis-script.py ${config_name}

model=CESM1-WACCM
ens=r1i1p1                                                                                                                                               
config_name=config_${model}${ens}_${variable}_${s_year}-${e_year}.yaml                                                                              
variable=${variable} sim=${sim} model=${model} ens=${ens} s_year=${s_year} e_year=${e_year} output_folder=${output_folder} reg_folder=${reg_folder} ./templater.sh config_template.yaml > ${config_name}                                                                                                                                   
echo "running "${config_name}
python3 solaris_heppa_analysis-script.py ${config_name}

ens=r2i1p1                                                                                                                                               
config_name=config_${model}${ens}_${variable}_${s_year}-${e_year}.yaml
variable=${variable} sim=${sim} model=${model} ens=${ens} s_year=${s_year} e_year=${e_year} output_folder=${output_folder} reg_folder=${reg_folder} ./templater.sh config_template.yaml > ${config_name}                                                                                                                         
echo "running "${config_name}
python3 solaris_heppa_analysis-script.py ${config_name}

ens=r3i1p1                                                                                                                                               
config_name=config_${model}${ens}_${variable}_${s_year}-${e_year}.yaml
variable=${variable} sim=${sim} model=${model} ens=${ens} s_year=${s_year} e_year=${e_year} output_folder=${output_folder} reg_folder=${reg_folder} ./templater.sh config_template.yaml > ${config_name}                                                                                                                                   
echo "running "${config_name}
python3 solaris_heppa_analysis-script.py ${config_name}
