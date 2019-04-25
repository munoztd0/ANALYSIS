#!/bin/bash


# set this to the directory containing antsRegistration
#ANTSPATH=/usr/local/ANTs/build/bin/
ANTSPATH=/usr/local/ants/bin/

# ITK thread count
ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=1
export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS

# Check args
if [ $# -lt 1 ]; then
echo "USAGE: $0 <subjID>"
  exit
fi

subjID=$1
echo subj: $subjID

# path to afine transform tool
c3d_affine_tool=/usr/local/c3d-1.1.0-Linux-gcc64/bin/c3d_affine_tool

# path to the warp tool
warpTool=/usr/local/ants/bin/WarpImageMultiTransform

# paths to the T1 structurals,
subAnatDir=/home/REWOD/DATA/STUDY/DERIVED/ICA_ANTS/sub-${subjID}/ses-second/anat/

# paths to the standard anatomical images
standardAnatDir=/home/REWOD/DATA/CANONICALS/

fixedT1=CIT168_T1w_MNI

# align the T1 anatomicals and downsample to functional resolution
echo "Running Flirt to downsample T1  $(date +"%T")"

flirt -ref ${standardAnatDir}${fixedT1} -in ${standardAnatDir}${fixedT1} -out ${standardAnatDir}${fixedT1}_lowres -applyisoxfm 1.8 -omat ${standardAnatDir}${fixedT1}_lowres.mat
echo "Done Flirt to downsample T1  $(date +"%T")"


#############
# Define target space (fixed) T1 & lowres T1 mask & normal masks file names
fixed_T1=${standardAnatDir}CIT168_T1w_MNI.nii.gz
fixed_T1lowres=${standardAnatDir}CIT168_T1w_MNI_lowres.nii.gz
fixed_mask=${standardAnatDir}`basename ${fixed_T1} .nii.gz`_mask.nii.gz

# Define subject images (moving) T1/T2 & masks file names
moving_T1=${subAnatDir}sub-${subjID}_ses-second_run-01_T1w_reoriented_brain.nii.gz
moving_mask=${subAnatDir}`basename ${moving_T1} .nii.gz`_mask.nii.gz

# Prefix for output transform files
outPrefix=${moving_T1%%.nii.gz}

#############
echo "anat dir: ${subAnatDir}"
echo "fixed T1: ${fixed_T1}"
echo "fixed mask: ${fixed_mask}"
echo "moving T1: ${moving_T1}"
echo "moving mask: ${moving_mask}"
echo "out prefix: ${outPrefix}"


# apply warp to the T1 anatomical (also check coreg quality)
echo "Apply ANTs wrap to T1 at $(date +"%T")"

${warpTool} 3 ${moving_T1} ${subAnatDir}sub-${subjID}_ses-second_run-01_T1w_reoriented_brain_ANTsCoreg.nii.gz \
        -R ${fixed_T1lowres} \
        ${outPrefix}_xfm1Warp.nii.gz ${outPrefix}_xfm0GenericAffine.mat

echo "Done ANTs warp to T1 at $(date +"%T")"


# resample the original so it can be compared (just ot check everything is alright)
echo "resample the original so it can be compared $(date +"%T")"

flirt -in ${moving_T1} -ref ${fixed_T1lowres} -dof 6 -out ${subAnatDir}sub-${subjID}_ses-second_run-01_T1w_reoriented_brain_standardAlign.nii.gz
echo "done resampling original $(date +"%T")"
