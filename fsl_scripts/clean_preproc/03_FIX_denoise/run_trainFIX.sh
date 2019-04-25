#!/bin/bash

# script to run
trainscript=home/REWOD/ANALYSIS/fsl_scripts/fsl_ANTS/clean_preproc/03_FIX_denoise/trainClassifier.sh


#submit to cluster (estimation will need to be adapted)
qsub -o /home/REWOD/ClusterOutput -j oe -l walltime=06:00:00,pmem=8GB -M david.munoz@etu.unige.ch -m e -l nodes=1 -q queue1 -N training_classifier ${trainscript}
