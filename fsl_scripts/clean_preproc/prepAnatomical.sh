#!/bin/bash

# ID for the subject we're working
#subjectID=$1
subjectID=$1
echo "Preparing subject ${subjectID} for FEAT"

# Directory containing un-processed nifti data
dataDir=/home/REWOD/DATA/STUDY/RAW/sub-${subjectID}/ses-second/
# Output directory for preprocessed files
outDir=/home/REWOD/DATA/STUDY/DERIVED/ICA_ANTS/sub-${subjectID}/ses-second/anat/


# make the subject level directory
mkdir -p ${outDir}


# ###################
# T1 SCAN: Reorient and extract brain
# Expects a file called T1 in the source directory
echo "Started working on T1 scan at $(date +"%T")"
# Reorient T1 scan to standard, and extract brain
fslreorient2std ${dataDir}/anat/sub-${subjectID}_ses-second_run-01_T1w.nii.gz ${outDir}sub-${subjectID}_ses-second_run-01_T1w_reoriented
bet ${outDir}sub-${subjectID}_ses-second_run-01_T1w_reoriented ${outDir}sub-${subjectID}_ses-second_run-01_T1w_reoriented_brain -f 0.2 -B
echo "Done working on T1 scan at $(date +"%T")"

# ###################
# T2 SCAN: Reorient and extract brain
# Expects a file called T2 in the source directory
#echo "Started working on T2 scan at $(date +"%T")"
# Reorient T2 scan to standard, and extract brain
#fslreorient2std ${dataDir}ses-first/anat/*T2.nii.gz ${outDir}sub-${subjectID}_ses-first_run-01_T2_reoriented
#bet ${outDir}sub-${subjectID}_ses-first_run-01_T2_reoriented ${outDir}sub-${subjectID}_ses-first_run-01_T2_reoriented_brain -f 0.2 -B
#echo "Done working on T2 scan at $(date +"%T")"
