#!/bin/bash

codeDir="/home/REWOD/ANALYSIS/spm_scripts/GLM/hedonic"
matlab_script="GLM_07_stLevel"
matlabSubmit="/home/REWOD/ANALYSIS/spm_scripts/matlab_oneSubj.sh"


# Loop over subjects
# no variance 03 - 04 - 13 - 20 - 23
for subj in 01 02 03 04 05 06 07 09 10 11 12 13 14 15 16 17 18 20 21 22 23 24 25 26
do
	# prep for each session's data
		qsub -o /home/REWOD/ClusterOutput -j oe -l walltime=0:40:00,pmem=4GB -M david.munoz@etu.unige.ch -m e -l nodes=1  -q queue1 -N HED_GLM-07_s-${subj} -F "${subj} ${codeDir} ${matlab_script}" ${matlabSubmit}

done
