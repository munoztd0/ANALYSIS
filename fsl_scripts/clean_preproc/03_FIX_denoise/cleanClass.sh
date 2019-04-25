#!/bin/bash

# pull in the subject we should be working on
subjectID=$1

#choose task OR runID=$2
taskID=$2

echo "Preparing subject ${subjectID} task ${taskID}"


# Directory containing functionals, high-res reference scans after ICA
funcDir=/home/REWOD/DATA/STUDY/DERIVED/ICA_ANTS/sub-${subjectID}/ses-second/func/task-${taskID}.ica

# the component classification rejection threshold // from the REWDOdata
threshold=20


###################
# classify components as artifact/signal, remove artifacts
# generates filtered_func_data_clean in the ICA directory
echo "started classification at $(date +"%T")"

# the classifier path
classifierPath=/home/shared/fix/classifiers/FIX_REWOD.RData

# classify components (approx 30 min)
# will generate a label file at the specified threshold called fix4melview_FIX_REWOD_thr<threshold>.txt
/usr/local/fix1.066/fix -c ${funcDir} ${classifierPath} ${threshold}

echo "Classification done at $(date +"%T")"

# remove bad ones (using the manually corrected classification) and filters out movement (flag -m)
/usr/local/fix1.066/fix -a ${funcDir}/fix4melview_FIX_REWOD_thr${threshold}.txt -m


echo "finished cleanup at $(date +"%T")"
