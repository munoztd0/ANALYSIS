#!/bin/bash
# session level script
melodicScript=/home/REWOD/ANALYSIS/fsl_scripts/clean_preproc/02_MELODIC_ICA/melodicICA.sh
## sessionScript=/home/evapool/PAVMOD/ANALYSIS/fsl_script/FSL_ANTS_tolman/melodicICA.sh

# Loop over
for subjectID in 24 #06 07 09 10 11 15 16 17 18 20 # 01 02 03 04 05 12 13 14#10 11 12 13 14  21 22 23 24 25 26 14 15 16 17 18 19 20 21 22 23 24 25 26 #
do
	# Loop over runs
for taskID in hedonic PIT
do

	#for session in second
		#do
			# EXTRA LONG #spawn session jobs to the cluster after the subject level work is complete
      qsub -o /home/REWOD/ClusterOutput -j oe -l walltime=40:00:00,pmem=10GB -M david.munoz@etu.unige.ch -m e -l nodes=1 -q queue1 -N ICA_${subjectID}_${taskID} -F "${subjectID} ${taskID} " ${melodicScript}
			## qsub -o ~/ClusterOutput -j oe -l nodes=1,walltime=10:00:00,pmem=5GB -M evapool@caltech.edu -m e -q batch -N ICA_${subj}_${session} -F "${subj} ${session}" ${sessionScript}
	done
done
#done
