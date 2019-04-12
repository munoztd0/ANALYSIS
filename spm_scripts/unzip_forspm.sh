#!/bin/bash


dataDir="/home/REWOD/DATA/STUDY/CLEAN/"
codeDir="/home/REWOD/ANALYSIS/spm_scripts/"
#sessionID="second"

# Loop over subjects

for subjID in 01 02 03 04 05 06 07 09 10 11 12 13 14 15 16 17 18 20 21 22 23 24 25 26
do
	# unzip functional
	#for runID in PIT pavlovianlearning hedonicreactivity
	#do
		# unzip for use in SPM
		 #echo "Expanding functionals Subject ${subjID} Run ${runID} at $(date +"%T")"
		 #gunzip -f ${dataDir}sub-${subjID}/ses-${sessionID}/func/sub-${subjID}_ses-${sessionID}_task-${runID}_run-01_bold.nii.gz
		 #echo "Done expanding functional Subject ${subjID} Session ${runID} at $(date +"%T")"
	#done

  # unzip anatomical
	  echo "Expanding anatomicals Subject ${subjID} at $(date +"%T")"
		gunzip -f ${dataDir}sub-${subjID}/anat/sub-${subjID}_ses-second_T1w.nii.gz
		echo "Done expanding anatomical Subject ${subjID} at $(date +"%T")"
done
