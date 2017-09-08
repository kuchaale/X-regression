#!/bin/bash

export n_samples=5000 #
export reg_dir='../regressors/'
export nc_gen=True
export pdf_gen=True
s_year=( 1980 )
i_year=1980
e_year=2010
config=all_trend
vars=( zmta zmua ) 
#r1i1p1 models
models=( CMAM IPSL MRI-ESM1r1 CCSRNIES CHASER-MIROC-ESM EMAC-L90MA EMAC-L47MA CESM1-WACCMSD )
#models=( CCSRNIES )
for i in ${models[@]}; do
	export in_dir='/mnt/4data/CCMI/'${i}'/'
	export out_dir='../regressors/'
	time python eofs_univ_ccmi_xarray.py -v ${i} zmua ${i_year} ${i_year} ${e_year} zmua_monthly_${i}_refC1SD_r1i1p1_${i_year}01-${e_year}12.nc
	export out_dir='/var/www/ccmi-sd/'
	for var in ${vars[@]}; do
		time python lin_reg_univ_ccmi_xarray.py -v ${i} ${var} ${i_year} ${s_year[@]} ${e_year}  ${var}_monthly_${i}_refC1SD_r1i1p1_${i_year}01-${e_year}12.nc all_trend
	done
done

#r1i1p2 models
models=( CNRM-CM5-3 EMAC-L90MA EMAC-L47MA )
for i in ${models[@]}; do
	export in_dir='/mnt/4data/CCMI/'${i}'/'
	export out_dir='../regressors/p2/'
	time python eofs_univ_ccmi_xarray.py -v ${i} zmua ${i_year} ${i_year} ${e_year} zmua_monthly_${i}_refC1SD_r1i1p2_${i_year}01-${e_year}12.nc
	export out_dir='/var/www/ccmi-sd/p2/'
	export reg_dir='../regressors/p2/'
	for var in ${vars[@]}; do
		time python lin_reg_univ_ccmi_xarray.py -v ${i} ${var} ${i_year} ${s_year[@]} ${e_year}  ${var}_monthly_${i}_refC1SD_r1i1p2_${i_year}01-${e_year}12.nc all_trend
	done
done
