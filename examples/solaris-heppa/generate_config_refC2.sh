#!/bin/bash

output_folder=/mnt/2data/solarisheppa/outputs_new/
reg_folder=/mnt/2data/solarisheppa/regressors/bernd/
sim=$4 #refC2

s_year=$1 #1960
e_year=$2 #2010
variable=$3 #zmo3,zmnoy

if [ ${variable} = "zmta,zmua" ] || [ ${variable} = "zmua,zmta" ]  
then
	epp_reg=epp_for_t_and_u
	spe_reg=spe_for_t_and_u
elif [ ${variable} = "zmo3,zmnoy" ] || [ ${variable} = "zmnoy,zmo3" ]
then
	epp_reg=eep_for_o3_and_noy                                                                                                                      
	spe_reg=spe_for_o3_and_noy
else
	echo "none of these options: "${variable}" is valid"
fi


model=SOCOL3
ens=r1i1p1
config_name=config_${model}_${sim}_${ens}_${variable}_${s_year}-${e_year}.yaml
variable=${variable} sim=${sim} model=${model} ens=${ens} s_year=${s_year} e_year=${e_year} output_folder=${output_folder} reg_folder=${reg_folder} epp_reg=${epp_reg} spe_reg=${spe_reg} ./templater.sh config_template.yaml > ${config_name}
echo "running "${config_name}
python3 solaris_heppa_analysis-script.py ${config_name}

model=CESM1-WACCM
ens=r1i1p1                                                                                                                                               
config_name=config_${model}_${sim}_${ens}_${variable}_${s_year}-${e_year}.yaml                                                                              
variable=${variable} sim=${sim} model=${model} ens=${ens} s_year=${s_year} e_year=${e_year} output_folder=${output_folder} reg_folder=${reg_folder} epp_reg=${epp_reg} spe_reg=${spe_reg} ./templater.sh config_template.yaml > ${config_name}                                                                                                                                   
echo "running "${config_name}
python3 solaris_heppa_analysis-script.py ${config_name}

ens=r2i1p1                                                                                                                                               
config_name=config_${model}_${sim}_${ens}_${variable}_${s_year}-${e_year}.yaml
variable=${variable} sim=${sim} model=${model} ens=${ens} s_year=${s_year} e_year=${e_year} output_folder=${output_folder} reg_folder=${reg_folder} epp_reg=${epp_reg} spe_reg=${spe_reg} ./templater.sh config_template.yaml > ${config_name}                                                                                                                         
echo "running "${config_name}
python3 solaris_heppa_analysis-script.py ${config_name}

ens=r3i1p1                                                                                                                                               
config_name=config_${model}_${sim}${ens}_${variable}_${s_year}-${e_year}.yaml
variable=${variable} sim=${sim} model=${model} ens=${ens} s_year=${s_year} e_year=${e_year} output_folder=${output_folder} reg_folder=${reg_folder} epp_reg=${epp_reg} spe_reg=${spe_reg} ./templater.sh config_template.yaml > ${config_name}                                                                                                                                   
echo "running "${config_name}
python3 solaris_heppa_analysis-script.py ${config_name}
