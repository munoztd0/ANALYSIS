#!/bin/bash

codeDir="/home/REWOD/ANALYSIS/spm_scripts/GLM/hedonicreactivity"
matlab_script="GLM_01_stLevel"
matlabSubmit="/home/REWOD/ANALYSIS/spm_scripts/matlab_oneSubj.sh"


# Loop over subjects

for subID in 03 #01 02 03 04 05 06 07 09 10 11
do
	# prep for each session's data
		qsub -o /home/REWOD/ClusterOutput -j oe -l walltime=4:00:00,pmem=4GB -M david.munoz@etu.unige.ch -m e -l nodes=1  -q queue1 -N GLM-1_sub-${sbj} -F "${sbj} ${codeDir} ${matlab_script}" ${matlabSubmit}

done
