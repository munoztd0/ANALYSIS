#!/bin/bash

codeDir="/home/REWOD/ANALYSIS/spm_scripts/GLM/PIT"
matlab_script="GLM_02_stLevel"
matlabSubmit="/home/REWOD/ANALYSIS/spm_scripts/matlab_oneSubj.sh"


# Loop over subjects

for subj in 01 #02 03 04 05 06 07 09 10 11 12 13 14 15 16 17 18 20 21 22 23 24 25 26
do
	# prep for each session's data
		qsub -o /home/REWOD/ClusterOutput -j oe -l walltime=4:00:00,pmem=8GB -M david.munoz@etu.unige.ch -m e -l nodes=1  -q queue1 -N PIT_GLM-02_sub-${subj} -F "${subj} ${codeDir} ${matlab_script}" ${matlabSubmit}

done
