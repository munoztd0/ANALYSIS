#!/bin/bash

codeDir="/home/OBIWAN/ANALYSIS/spm_scripts/GLM/hedonicreactivity"
matlab_script="GLM_02_stLevel"
matlabSubmit="/home/OBIWAN/ANALYSIS/spm_scripts/matlab_oneSubj.sh"


# Loop over subjects

for subj in obese224 obese225 obese226 obese227
do
	# prep for each session's data
		qsub -o /home/OBIWAN/ClusterOutput -j oe -l walltime=4:00:00,pmem=4GB -M eva.pool@unige.ch -m e -l nodes=1  -q queue1 -N GLM-3_sub-${subj} -F "${subj} ${codeDir} ${matlab_script}" ${matlabSubmit}

done
