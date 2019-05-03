#!/bin/bash

codeDir="/home/REWOD/ANALYSIS/spm_scripts/copy_spm_scripts/GLM/PIT"
matlab_script="GLM_01_stLevel"
matlabSubmit="/home/REWOD/ANALYSIS/spm_scripts/copy_spm_scripts/matlab_oneSubj.sh"


# Loop over subjects

for subj in 02
do
	# prep for each session's data
		qsub -o /home/REWOD/ClusterOutput -j oe -l walltime=4:00:00,pmem=4GB -M david.munoz@etu.unige.ch -m e -l nodes=1  -q queue1 -N GLM-3_sub-${subj} -F "${subj} ${codeDir} ${matlab_script}" ${matlabSubmit}

done
