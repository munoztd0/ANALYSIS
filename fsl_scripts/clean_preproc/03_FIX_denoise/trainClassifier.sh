#!/bin/bash

#name the classfier
classifierName=FIX_REWOD

# REWOD datadir // After have moved all the hand_labels_noise into the REWODdataf folder
#OR REWODdata=/home/REWOD/ANALYSIS/fsl_scripts/clean_preproc/03_FIX_denoise/
REWODdata=$(find /home/shared/fix/data/REWOD/ -type d  -name "melodic_run.ica")
#mkdir -p /home/shared/fix/data/REWOD/

echo "data used to train the classifier: ${REWODdata}"
# train the classifier: will generate a an Rdata file that is the classifier
echo "training classifier started at $(date +"%T")"
/usr/local/fix1.066/fix -t ${classifierName} -l ${REWODdata}
echo "training classifier finished at $(date +"%T")"
