#!/bin/bash

codeDir=/home/OBIWAN/ANALYSIS/fsl_scripts/preproc/fsl_ANTS/
##codeDir=/home/evapool/PAVMOD/ANALYSIS/fsl_script/FSL_ANTS_tolman/
anatomicalScript=${codeDir}prepAnatomical.sh
functionalScript=${codeDir}prepFunctional.sh

# Loop over subjects

# for subj in control105
 for subj in control126

do
	# work on each subject's anatomical scans
	qsub -o /home/OBIWAN/ClusterOutput -j oe -l walltime=2:00:00,pmem=6GB -M eva.pool@unige.ch -m e -l nodes=1 -q queue1 -N prepFEAT_Subject_sub-${subj} -F "${subj}" ${anatomicalScript}
##	qsub -o ~/ClusterOutput -j oe -l walltime=1:00:00 -M evapool@caltech.edu -m e -l nodes=1 -q batch -N prepFEAT_Subject_${subj} -F "${subj}" ${anatomicalScript}

	# prep for each session's data
  for taskID in pavlovianlearning PIT hedonicreactivity
#	for runID in 01
	do
			qsub -o /home/OBIWAN/ClusterOutput -j oe -l walltime=2:00:00,pmem=6GB -M eva.pool@unige.ch -m e -l nodes=1 -q queue1 -N prepFEAT_Subject_sub-${subj}_task-${taskID} -F "${subj} ${taskID}" ${functionalScript}
#		qsub -o ~/ClusterOutput -j oe -l walltime=0:10:00 -M evapool@caltech.edu -m e -l nodes=1 -q batch -N prepFEAT_Subject_${subj}_${runID} -F "${subj} ${runID}" ${functionalScript}

	done
done
