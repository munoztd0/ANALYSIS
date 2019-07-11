#!/bin/bash

codeDir="/home/REWOD/ANALYSIS/spm_scripts/GLM/hedonic"
matlab_script="GLM_03b_ndLevel"
matlabSubmit="/home/REWOD/ANALYSIS/spm_scripts/matlab_oneScript.sh"

qsub -o /home/REWOD/ClusterOutput -j oe -l walltime=0:20:00,pmem=2GB -M david.munoz@etu.unige.ch -m e -q queue1 -N HED_GLM-03b_ttests -F " ${codeDir} ${matlab_script}" ${matlabSubmit}
#qsub -o ~/ClusterOutput -j oe -l walltime=2:00:00,pmem=4GB -M david.munoz@etu.unige.ch -m e -q queue1 -N GLM-03i_sub-${subj} -F "${subj} ${codeDir} ${matlab_script}" ${matlabSubmit}
