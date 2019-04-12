#!/usr/bin/env python

import math
import os, sys
from nibabel import load
import nipype.interfaces.io as nio           # Data i/o
import nipype.interfaces.fsl as fsl          # fsl
import nipype.interfaces.utility as util     # utility
import nipype.pipeline.engine as pe          # pypeline engine
import nipype.algorithms.modelgen as model   # model generation
import nipype.algorithms.rapidart as ra      # artifact detection
import nibabel as nib
import nipype.interfaces.c3 as c3
import nipype.interfaces.ants as ants
from nipype.interfaces.fsl import fix        # fix automated denoising
from nipype.interfaces.utility import Function

# location (str) where the data are stored (should have sub-directories mri, openfmri)
data_dir = None

if not data_dir:
    print("please look at this file (sys.argv[0]), you still need to provide some more info")






# ---------------------------- Aux functions 

def Get_nVols(filename):
    """ Get the number of volumes in a 4D nifti file
    
    Arguments:
        filename {str} -- filename
    
    Returns:
        integer -- number of volumes in 4D nifti file
    """
    from nipype import logging
    import nibabel as nib
    v = nib.load(filename)
    data = v.get_data()
    npts = data.shape[3]
    iflogger = logging.getLogger('interface')
    iflogger.info("npts: %s" % npts)
    return npts

def Get_TotalVoxels(filename):
    """ Get the total number of voxels in a 4D file, i.e. x*y*z*t
    
    Arguments:
        filename {filename} -- filename
    
    Returns:
        integer -- the total number of voxels in 4d nifti file
    """

    from nipype import logging
    import nibabel as nib
    import numpy as np
    v = nib.load(filename)
    data = v.get_data()
    total_voxels = np.prod(data.shape)
    iflogger = logging.getLogger('interface')
    iflogger.info("Total Voxels: %s" % total_voxels)
    return total_voxels


def get_voxel_size(filename):
    """ determine resolution of nifti file
    
    Arguments:
        filename {str} -- filename
    
    Returns:
        float -- voxel size in mm
    """
    img_header = load(filename).get_header().get_zooms()
    return float(img_header[0])


def Prepare_Design_FSF(feat_files, initial_highres_files, highres_files, npts=99, total_voxels=191):
    """ This functions loads a template melodic.fsf file, and customizes it for the data of one participant.
    
        Sample call could be something like this
        feat_files = "/home/pauli/Development/core_shell/data/mri/JOD-WP-CS2-005/f1_short"
        initial_highres_files = "/home/pauli/Development/core_shell/data/mri/JOD-WP-CS2-005/f1_ref"
        highres_files = "/home/pauli/Development/core_shell/data/mri/JOD-WP-CS2-005/t2"
        total_voxels = 35000000
        npts = 100
        design_fsf = prepare_design_fsf(feat_files, initial_highres_files, highres_files, npts, total_voxels)

    Arguments:
        feat_files {str} -- filename of functional file
        initial_highres_files {str} -- filename of highres T1 anatomical file
        highres_files {str} -- filename of highres T2 file
    
    Keyword Arguments:
        npts {int} -- number of volumes in functional (default: {99})
        total_voxels {int} -- total number of voxels, i.e. x*y*x*t (default: {191})
    
    Returns:
        str -- filename of created design.fsf file
    """

    from nipype import logging
    import os
    # open design.fsf template 
    design_fsf = open('./templates/melodic.fsf','r')
    
    # assemble dicts from inputs
    str_rep = dict()
    str_rep['feat_files(1)'] = feat_files.strip('.nii.gz')
    str_rep['initial_highres_files(1)'] = initial_highres_files.strip('.nii.gz')
    str_rep['highres_files(1)'] = highres_files.strip('.nii.gz')
    int_rep = dict()
    int_rep['fmri(totalVoxels)'] = total_voxels
    int_rep['fmri(npts)'] = npts

    out_lines = []
    lines = design_fsf.readlines()
    for line in lines:
        items = line.split(' ')
        if len(items) == 3:
            if str_rep.has_key(items[1]): # the second item in the line matches one of the dict keys with a string value
                # print line
                out_line = "%s %s \"%s\"\n" % (items[0], items[1], str_rep[items[1]])
                # print out_line
            elif int_rep.has_key(items[1]): # the second item in the line matches one of the dict keys with a int value
                # print line
                out_line = "%s %s %s\n" % (items[0], items[1], int_rep[items[1]])
                # print out_line
            else:
                out_line = line
            out_lines.append(out_line)
        else:
            out_lines.append(line)
    iflogger = logging.getLogger('interface')
    out_file = os.path.abspath('design.fsf')
    iflogger.info(out_file)
    with open(out_file, 'w') as fp:
       fp.writelines(out_lines)
    return out_file



def unlist(mylist):
    """ returns first element of a list of filename. Hard to explain why we need this, but you'll see below how it is used.
    
    Arguments:
        mylist {list} -- a list of filenames
    
    Returns:
        str -- a filename
    """
    r = mylist[0]
    return r


def reverse(mylist):
    """ reverses a list, no magic here
    
    Arguments:
        mylist {list} -- a list
    
    Returns:
        list -- the input list, reversed, amazing
    """
    mylist.reverse()
    return mylist



def flatten(l):
    """ turn 2D list into 1D
    
    Arguments:
        l {list} -- a list
    """
    l = sum(l, [])
    return(l)

def to_opt_string(out_stat):
    """ this function is used to create the input arguments for a fslmath command, using the minimum act found somewhere in the func
    
    Arguments:
        out_stat {float/integer} -- minimum act in the fun
    
    Returns:
        str -- option string for fslmask command
    """

    r = []
    for i in xrange(len(out_stat)):
        r.append("-thr %s -bin -fillh -dilF" % (float(out_stat[i][0]) + .5))
    return r

fsl.FSLCommand.set_default_output_type('NIFTI_GZ')







# ------------------------------ data definition and grabbing


# create preproc workflow
preproc = pe.Workflow(name='preproc')
preproc.base_dir = os.path.abspath('./nipype') #  results of preprocessing will be stored here
ds_ventricle_mask = os.path.join(data_dir, 'openfmri/group/ventricle_mask.nii.gz')

# List of subjects with usable data
subjects = ['JOD-WP-CS2-005', 'JOD-WP-CS2-006', 'JOD-WP-CS2-007','JOD-WP-CS2-009', 'JOD-WP-CS2-010', 'JOD-WP-CS2-011', 'JOD-WP-CS2-012', 'JOD-WP-CS2-013', 'JOD-WP-CS2-014', 'JOD-WP-CS2-015', 'JOD-WP-CS2-016', 'JOD-WP-CS2-017', 'JOD-WP-CS2-018' , 'JOD-WP-CS2-019', 'JOD-WP-CS2-020', 'JOD-WP-CS2-021', 'JOD-WP-CS2-022', 'JOD-WP-CS2-023', 'JOD-WP-CS2-024', 'JOD-WP-CS2-025', 'JOD-WP-CS2-026', 'JOD-WP-CS2-027', 'JOD-WP-CS2-028', 'JOD-WP-CS2-029' , 'JOD-WP-CS2-030' ]

# dictionary with info for what files are associate with each participant. (func(tional), func_ref(whole brain functional for referene), t1, t2, fm_pos/neg is for fielmap generation)
info = dict(func=[['subject_id', ['f1', 'f2', 'f3', 'f4']]],
            func_ref=[['subject_id', ['f1_ref','f2_ref','f3_ref','f4_ref']]],
            t1=[['subject_id', 't1']],
            t2=[['subject_id', 't2']],
            fm_pos=[['subject_id', ['f1_pos','f2_pos','f3_pos','f4_pos']]],
            fm_neg=[['subject_id', ['f1_neg','f2_neg','f3_neg','f4_neg']]])

datasource = pe.Node(interface=nio.DataGrabber(infields=['subject_id'], outfields=['func', 'func_ref', 't1', 't2', 'fm_pos', 'fm_neg']), name='datasource')
datasource.inputs.base_directory = os.path.abspath(os.path.join(data_dir, 'mri/'))
datasource.inputs.template = '%s/%s.nii.gz'
datasource.inputs.template_args = info
datasource.iterables = ('subject_id', subjects)
datasource.inputs.sort_filelist = True


inputnode = pe.Node(interface=util.IdentityInterface(fields=['func', 'func_ref', 't1', 't2', 'fm_pos', 'fm_neg']), name='inputnode')

preproc.connect(datasource, 'func', inputnode, 'func')
preproc.connect(datasource, 'func_ref', inputnode, 'func_ref')
preproc.connect(datasource, 't1', inputnode, 't1')
preproc.connect(datasource, 't2', inputnode, 't2')
preproc.connect(datasource, 'fm_pos', inputnode, 'fm_pos')
preproc.connect(datasource, 'fm_neg', inputnode, 'fm_neg')







# ------------- prepare anatomical data -------------

prep_anatomicals = pe.Workflow(name='prep_anatomicals')
prep_anatomicals.base_dir = os.path.abspath('./nipype')

# reorient t1 to standard
T1_to_standard = pe.Node(interface=fsl.Reorient2Std(output_type = "NIFTI_GZ"), name='T1_to_standard') 
preproc.connect(inputnode, 't1', T1_to_standard, 'in_file')


# reorient t2 to standard
T2_to_standard = pe.Node(interface=fsl.Reorient2Std(output_type = "NIFTI_GZ"), name='T2_to_standard') 
preproc.connect(inputnode, 't2', T2_to_standard, 'in_file')


# coreg t2 to t1
t2tot1 = pe.Node(interface=fsl.FLIRT(dof=6, output_type = "NIFTI_GZ"), name='t2tot1')
preproc.connect(T2_to_standard, 'out_file', t2tot1, 'in_file')
preproc.connect(T1_to_standard, 'out_file', t2tot1, 'reference')
    

# downsample caltech atlas t1 to t1/t2 resolution
downsample_atlas_t1 = pe.Node(interface=fsl.FLIRT(dof=6, output_type = "NIFTI_GZ", in_file=os.path.join(data_dir, 'CIT_Midbrain_Atlas/CIT168_T1w_700um_MNI.nii.gz', reference='/home/pauli/Development/core_shell/data/CIT_Midbrain_Atlas/CIT168_T1w_700um_MNI.nii.gz', in_matrix_file='/usr/share/fsl/5.0/etc/flirtsch/ident.mat'), name='downsample_atlas_t1')
preproc.connect(inputnode, ('t1', get_voxel_size), downsample_atlas_t1, 'apply_isoxfm')

downsample_atlas_t2 = pe.Node(interface=fsl.FLIRT(dof=6, output_type = "NIFTI_GZ", in_file=os.path.join(data_dir, 'CIT_Midbrain_Atlas/CIT168_T2w_700um_MNI.nii.gz', reference='/home/pauli/Development/core_shell/data/CIT_Midbrain_Atlas/CIT168_T2w_700um_MNI.nii.gz', in_matrix_file='/usr/share/fsl/5.0/etc/flirtsch/ident.mat'), name='downsample_atlas_t2')
preproc.connect(inputnode, ('t2', get_voxel_size), downsample_atlas_t2, 'apply_isoxfm')

# # create mask file for atlas
mask_atlas = pe.Node(interface=fsl.ImageMaths(op_string = '-thr .001 -bin -dilF', suffix = '_mask'), name='mask_atlas')
preproc.connect(downsample_atlas_t1, 'out_file', mask_atlas, 'in_file')

# skull strip t1
skullstrip = pe.Node(interface=fsl.BET(mask = True, frac=0.2, reduce_bias = True), name = 'stripstruct')
preproc.connect(T1_to_standard, 'out_file', skullstrip, 'in_file')

# we inflate our mask a bit, to make sure we don't remove any tissue by accident
inflate_mask = pe.Node(interface=fsl.ImageMaths(op_string = '-dilF', suffix = '_dilF'), name='inflate_mask')
preproc.connect(skullstrip, 'mask_file', inflate_mask, 'in_file')

# apply inflated mask to t2
maskT2 = pe.Node(interface=fsl.ImageMaths(op_string = '-mas', suffix = '_bet'), name='maskT2')
preproc.connect(t2tot1, 'out_file', maskT2, 'in_file')
preproc.connect(inflate_mask, 'out_file', maskT2, 'in_file2')

# apply inflated mask to t1
maskT1 = pe.Node(interface=fsl.ImageMaths(op_string = '-mas', suffix = '_bet'), name='maskT1')
preproc.connect(T1_to_standard, 'out_file', maskT1, 'in_file')
preproc.connect(inflate_mask, 'out_file', maskT1, 'in_file2')

# affine co-reg (flirt) of T1 to caltech atlas
init_t1_to_atlas_coreg = pe.Node(interface=fsl.FLIRT(dof=6, output_type = "NIFTI_GZ"), name='init_t1_to_atlas_coreg')
preproc.connect(maskT1, 'out_file', init_t1_to_atlas_coreg, 'in_file')
preproc.connect(downsample_atlas_t1, 'out_file', init_t1_to_atlas_coreg, 'reference')

# convert transform matrix to ants format (ITK)
fsl2ras = pe.Node(interface=c3.C3dAffineTool(fsl2ras = True, itk_transform=True), name='fsl2ras')
preproc.connect(downsample_atlas_t1, 'out_file', fsl2ras, 'reference_file')
preproc.connect(maskT1, 'out_file', fsl2ras, 'source_file')
preproc.connect(init_t1_to_atlas_coreg, 'out_matrix_file', fsl2ras, 'transform_file')

# merge fixed and moving images into list
merge_fixed = pe.Node(interface=util.Merge(2, axis='hstack'), name='merge_fixed')
preproc.connect(downsample_atlas_t1, 'out_file', merge_fixed, 'in1')
preproc.connect(downsample_atlas_t2, 'out_file', merge_fixed, 'in2')
merge_moving = pe.Node(interface=util.Merge(2, axis='hstack'), name='merge_moving')
preproc.connect(maskT1, 'out_file', merge_moving, 'in1')
preproc.connect(maskT2, 'out_file', merge_moving, 'in2')

# assemble the inputs for coregistration of t2/t1 to caltech atlas
structs_to_atlas_coreg_input = pe.Node(interface=util.IdentityInterface(fields=['fixed_image', 'moving_image', 'fixed_image_mask', 'moving_image_mask', 'initial_moving_transform']), name='structs_to_atlas_coreg_input')
preproc.connect(merge_fixed, ('out', unlist), structs_to_atlas_coreg_input, 'fixed_image')
preproc.connect(merge_moving, ('out', unlist), structs_to_atlas_coreg_input, 'moving_image')
preproc.connect(mask_atlas, 'out_file', structs_to_atlas_coreg_input, 'fixed_image_mask')
preproc.connect(inflate_mask, 'out_file', structs_to_atlas_coreg_input, 'moving_image_mask')
preproc.connect(fsl2ras, 'itk_transform', structs_to_atlas_coreg_input, 'initial_moving_transform')

# perform diffeomorphic mapping of t2/t1 to caltech atlas
structs_to_atlas_coreg = pe.Node(interface=ants.Registration(dimension=3, transforms=['SyN'], metric=[['CC'] * 2], radius_or_number_of_bins=[[4] * 2], metric_weight = [[.5] * 2], transform_parameters=[[.1,3.0,0.0]], number_of_iterations=[[100,100,70,50,20]], convergence_threshold=[1.e-6], convergence_window_size=[10], smoothing_sigmas=[[5.0,3.0,2.0,1.0,0.0]], shrink_factors=[[10,6,4,2,1]], use_histogram_matching=True, interpolation='Linear', invert_initial_moving_transform=False, sampling_strategy = [['Random']*2], sampling_percentage = [[0.05]*2], verbose=True), name='structs_to_atlas_coreg')
preproc.connect(structs_to_atlas_coreg_input, 'fixed_image', structs_to_atlas_coreg, 'fixed_image')
preproc.connect(structs_to_atlas_coreg_input, 'moving_image', structs_to_atlas_coreg, 'moving_image')
preproc.connect(structs_to_atlas_coreg_input, 'fixed_image_mask', structs_to_atlas_coreg, 'fixed_image_mask')
preproc.connect(structs_to_atlas_coreg_input, 'moving_image_mask', structs_to_atlas_coreg, 'moving_image_mask')
preproc.connect(structs_to_atlas_coreg_input, 'initial_moving_transform', structs_to_atlas_coreg, 'initial_moving_transform')

# configure the output of above diff mapp process
structs_to_atlas_coreg_output = pe.Node(interface=util.IdentityInterface(fields=['forward_transforms', 'warped_image', 'reverse_transforms', 'inverse_warped_image']), name='structs_to_atlas_coreg_output')
preproc.connect(structs_to_atlas_coreg, 'forward_transforms', structs_to_atlas_coreg_output, 'forward_transforms')
preproc.connect(structs_to_atlas_coreg, 'warped_image', structs_to_atlas_coreg_output, 'warped_image')
preproc.connect(structs_to_atlas_coreg, 'reverse_transforms', structs_to_atlas_coreg_output, 'reverse_transforms')
preproc.connect(structs_to_atlas_coreg, 'inverse_warped_image', structs_to_atlas_coreg_output, 'inverse_warped_image')










# --------------------- prepare functionals ---------------

# n4 bias correction of mean functional
n4biasCorrMbRef = pe.MapNode(interface=ants.N4BiasFieldCorrection(dimension = 3, bspline_fitting_distance = 300, shrink_factor = 3, n_iterations = [50,50,30,20], save_bias = True), name='n4biasCorrMbRef', iterfield=['input_image'])
preproc.connect(inputnode, 'func_ref', n4biasCorrMbRef, 'input_image')

# get the total unmber of voxels
get_totalVoxels = pe.MapNode(interface=Function(input_names=['in_file'], output_names=['output'], function=Get_TotalVoxels), name='get_totalVoxels', iterfield=['in_file'])
preproc.connect(inputnode, 'func', get_totalVoxels, 'in_file')

# get the number of volumes in a scan
get_nVols = pe.MapNode(interface=Function(input_names=['in_file'], output_names=['output'], function=Get_nVols), name='get_nVols', iterfield=['in_file'])
preproc.connect(inputnode, 'func', get_nVols, 'in_file')

# customize template melodic.fsf file for data of one participant
prepare_design_fsf = pe.MapNode(interface=Function(input_names=['feat_files','initial_highres_files','highres_files','npts','total_voxels'], output_names=['out_file'], function=Prepare_Design_FSF), name='prepare_design_fsf', iterfield=['feat_files','initial_highres_files', 'npts','total_voxels'])
preproc.connect(inputnode, 'func', prepare_design_fsf, 'feat_files')
preproc.connect(n4biasCorrMbRef, 'output_image', prepare_design_fsf, 'initial_highres_files')
preproc.connect(maskT2, 'out_file', prepare_design_fsf, 'highres_files')
preproc.connect(get_totalVoxels, 'output', prepare_design_fsf, 'total_voxels')
preproc.connect(get_nVols, 'output', prepare_design_fsf, 'npts')

# define the feat process for preprocessing the data
feat = pe.MapNode(interface=fsl.FEAT(), name='feat', iterfield=['fsf_file']) 
preproc.connect(prepare_design_fsf, 'out_file', feat, 'fsf_file')







# # ----------------- denoise ---------

# call the ica-fix feature extraction, mel_ica is the output directory of a melodic run
extract_features = pe.MapNode(interface=fix.FeatureExtractor(), name='extract_features', iterfield=['mel_ica'])
preproc.connect(feat, 'feat_dir', extract_features, 'mel_ica')

# assemble inputs for training the ica-fix classifier
training_input = pe.JoinNode(interface=util.IdentityInterface(fields=['mel_ica']), joinfield=['mel_ica'], joinsource='datasource', name='training_input')
preproc.connect(extract_features, 'mel_ica', training_input, 'mel_ica')

# create a training set
create_training_set = pe.Node(interface=fix.TrainingSetCreator(), name='trainingset_creator')
preproc.connect(training_input, ('mel_ica', flatten), create_training_set, 'mel_icas_in')

# do the training
train_node = pe.Node(interface=fix.Training(trained_wts_filestem='core_shell_py'), name='train_node')
preproc.connect(create_training_set, 'mel_icas_out', train_node, 'mel_icas')

# classify ica components as noise/signal
classify_node = pe.MapNode(interface=fix.Classifier(thresh=5), name='classify', iterfield=['mel_ica'])
preproc.connect(train_node, 'trained_wts_file', classify_node, 'trained_wts_file')
preproc.connect(feat, 'feat_dir', classify_node, 'mel_ica')

# remove noise components
cleaner_node = pe.MapNode(interface=fix.Cleaner(cleanup_motion=True,), name='cleaner', iterfield=['artifacts_list_file'])
preproc.connect(classify_node, 'artifacts_list_file', cleaner_node, 'artifacts_list_file')

# extract mean func
meanfunc = pe.MapNode(interface=fsl.ImageMaths(op_string = '-Tmean', suffix='_mean'), name='meanfunc', iterfield = ['in_file'])
preproc.connect(cleaner_node, 'cleaned_functional_file', meanfunc, 'in_file')

# get intensity range of mean func
func_intensity_range = pe.MapNode(interface=fsl.ImageStats(op_string= '-r'), name='func_intensity_range', iterfield=['in_file'])
preproc.connect(meanfunc, 'out_file', func_intensity_range, 'in_file')

# create func mask
func_mask = pe.MapNode(interface=fsl.ImageMaths(suffix='_mask'), name='func_mask', iterfield=['in_file', 'op_string'])
preproc.connect(meanfunc, 'out_file', func_mask, 'in_file')
preproc.connect(func_intensity_range, ('out_stat', to_opt_string), func_mask, 'op_string')

# n4 bias correction of mean
n4biasCorrMFunc = pe.MapNode(interface=ants.N4BiasFieldCorrection(dimension = 3, bspline_fitting_distance = 300, shrink_factor = 3, n_iterations = [50,50,30,20], save_bias = True), name='n4biasCorrMFunc', iterfield=['input_image', 'mask_image'])
preproc.connect(meanfunc, 'out_file', n4biasCorrMFunc, 'input_image')
preproc.connect(func_mask, 'out_file', n4biasCorrMFunc, 'mask_image')

# co-register mbRef to merged fm epis
mFunc_to_mbref = pe.MapNode(interface=fsl.FLIRT(dof=6, output_type = "NIFTI_GZ"), name='mFunc_to_mbref', iterfield=['in_file', 'reference'])
preproc.connect(n4biasCorrMFunc, 'output_image', mFunc_to_mbref, 'in_file')
preproc.connect(n4biasCorrMbRef, 'output_image', mFunc_to_mbref, 'reference')








# # # --------------------------  prepare fieldmaps -----------------------

# merge fm_pos and fm_neg into one list
mergenode = pe.Node(interface=util.Merge(2, axis='hstack'), name='merge', iterfield=['in1','in2'])
preproc.connect(inputnode, 'fm_pos', mergenode, 'in1')
preproc.connect(inputnode, 'fm_neg', mergenode, 'in2')

# concatenate pos and neg image into one for topup
mergeposneg = pe.MapNode(interface=fsl.Merge(dimension='t'), name='mergeposneg', iterfield=['in_files'])
preproc.connect(mergenode, 'out',  mergeposneg, 'in_files')

# create fieldmap and field_cof
# subtle change to configuration file, to use subsampling of 5, instead of 2, to get started. This is necessary because of the un-even number of slices (35)
topup = pe.MapNode(interface=fsl.TOPUP(config='/home/pauli/Development/core_shell/analysis/mri/b02b0_wolf.cnf', encoding_file = "/home/pauli/Development/core_shell/analysis/mri/topup_encoding.txt", output_type = "NIFTI_GZ"), name='topup', iterfield=['in_file'])
preproc.connect(mergeposneg, 'merged_file', topup, 'in_file')

# convert fielmaps to rad/s
fm_to_fmRads = pe.MapNode(interface=fsl.ImageMaths(op_string = ('-mul %s' % (math.pi * 2)), suffix = '_mask'), name='fm_fmRads', iterfield=['in_file'])
preproc.connect(topup, 'out_field', fm_to_fmRads, 'in_file')

# co-register mbRef to merged fm epis
mbRef_to_corrected_fm_epi = pe.MapNode(interface=fsl.FLIRT(dof=6, output_type = "NIFTI_GZ"), name='mbRef_to_corrected_fm_epi', iterfield=['in_file', 'reference'])
preproc.connect(n4biasCorrMbRef, 'output_image', mbRef_to_corrected_fm_epi, 'in_file')
preproc.connect(topup, 'out_corrected', mbRef_to_corrected_fm_epi, 'reference')

# unwarp co-registered mb reference image
unwarp = pe.MapNode(interface=fsl.FUGUE(save_shift = True, unwarp_direction = 'y', output_type = "NIFTI_GZ", dwell_time=.00072), name='unwarp', iterfield=['fmap_in_file', 'in_file'])
preproc.connect(fm_to_fmRads, 'out_file', unwarp, 'fmap_in_file')
preproc.connect(mbRef_to_corrected_fm_epi, 'out_file', unwarp, 'in_file')

# invert xfm 
invert_tfm = pe.MapNode(interface=fsl.ConvertXFM(output_type = "NIFTI_GZ", invert_xfm = True), name='invert_tfm', iterfield=['in_file'])
preproc.connect(mbRef_to_corrected_fm_epi, 'out_matrix_file', invert_tfm, 'in_file')

# rev apply mbRef to topup trans mat to out_shift_file
shift_file2mbRef = pe.MapNode(interface=fsl.FLIRT(apply_xfm = True, output_type = "NIFTI_GZ"), name='shift_file2mbRef', iterfield=['in_file','reference','in_matrix_file'])
preproc.connect(unwarp, 'shift_out_file', shift_file2mbRef, 'in_file')
preproc.connect(inputnode, 'func_ref', shift_file2mbRef, 'reference')
preproc.connect(invert_tfm, 'out_file', shift_file2mbRef, 'in_matrix_file')

# concat transformation file for mFunc to mbRef and mbRef to FM
concat_tfm_mFunc_mbRef_FM = pe.MapNode(interface=fsl.ConvertXFM(output_type = "NIFTI_GZ", concat_xfm = True), name='concat_tfm_mFunc_mbRef_FM', iterfield=['in_file', 'in_file2'])
preproc.connect(mbRef_to_corrected_fm_epi, 'out_matrix_file', concat_tfm_mFunc_mbRef_FM, 'in_file')
preproc.connect(mFunc_to_mbref, 'out_matrix_file', concat_tfm_mFunc_mbRef_FM, 'in_file2')

# invert concatenated transformation matrix
invert_tfm_mFunc_mbRef_FM = pe.MapNode(interface=fsl.ConvertXFM(output_type = "NIFTI_GZ", invert_xfm = True), name='invert_tfm_mFunc_mbRef_FM', iterfield=['in_file'])
preproc.connect(concat_tfm_mFunc_mbRef_FM, 'out_file', invert_tfm_mFunc_mbRef_FM, 'in_file')

# transform fm from FM to mFunc
shift_file2mFunc = pe.MapNode(interface=fsl.FLIRT(apply_xfm = True, output_type = "NIFTI_GZ"), name='shift_file2mFunc', iterfield=['in_file','reference','in_matrix_file'])
preproc.connect(unwarp, 'shift_out_file', shift_file2mFunc, 'in_file')
preproc.connect(cleaner_node, 'cleaned_functional_file', shift_file2mFunc, 'reference')
preproc.connect(invert_tfm_mFunc_mbRef_FM, 'out_file', shift_file2mFunc, 'in_matrix_file')







# # # ---------- unwarp functionals ---------------------------------

# unwarp mbRef in place
unwarp_mbRef = pe.MapNode(interface=fsl.FUGUE(save_shift = True, unwarp_direction = 'y', output_type = "NIFTI_GZ", dwell_time=.00072), name='unwarp_mbRef', iterfield=['shift_in_file','in_file'])
preproc.connect(shift_file2mbRef, 'out_file', unwarp_mbRef, 'shift_in_file')
preproc.connect(n4biasCorrMbRef, 'output_image', unwarp_mbRef, 'in_file')

# unwarp func in place
unwarp_func = pe.MapNode(interface=fsl.FUGUE(save_shift = True, unwarp_direction = 'y', output_type = "NIFTI_GZ", dwell_time=.00072), name='unwarp_func', iterfield=['shift_in_file','in_file'])
preproc.connect(shift_file2mFunc, 'out_file', unwarp_func, 'shift_in_file')
preproc.connect(cleaner_node, 'cleaned_functional_file', unwarp_func, 'in_file')








# # ----------- coregistation to standard

# resample caltech atlas T1 to resolution of functional (from size of structural)
downsample_atlas_t1_to_func = pe.Node(interface=fsl.FLIRT(dof=6, output_type = "NIFTI_GZ", reference='/home/pauli/Development/core_shell/data/CIT_Midbrain_Atlas/CIT168_T1w_700um_MNI.nii.gz', in_matrix_file='/usr/share/fsl/5.0/etc/flirtsch/ident.mat'), name='downsample_atlas_t1_to_func')
preproc.connect(inputnode, ('func', get_voxel_size), downsample_atlas_t1_to_func, 'apply_isoxfm')
preproc.connect(downsample_atlas_t1, 'out_file', downsample_atlas_t1_to_func, 'in_file')

# same for t2
downsample_atlas_t2_to_func = pe.Node(interface=fsl.FLIRT(dof=6, output_type = "NIFTI_GZ", reference='/home/pauli/Development/core_shell/data/CIT_Midbrain_Atlas/CIT168_T2w_700um_MNI.nii.gz', in_matrix_file='/usr/share/fsl/5.0/etc/flirtsch/ident.mat'), name='downsample_atlas_t2_to_func')
preproc.connect(inputnode, ('func', get_voxel_size), downsample_atlas_t2_to_func, 'apply_isoxfm')
preproc.connect(downsample_atlas_t2, 'out_file', downsample_atlas_t2_to_func, 'in_file')

# invert xfm 
downsample_atlas_t2_to_func_inv = pe.MapNode(interface=fsl.ConvertXFM(output_type = "NIFTI_GZ", invert_xfm = True), name='downsample_atlas_t2_to_func_inv', iterfield=['in_file'])
preproc.connect(downsample_atlas_t2_to_func, 'out_matrix_file', downsample_atlas_t2_to_func_inv, 'in_file')

# convert to itk format
fsl2ras_downsample_atlas_from_struct_to_func = pe.Node(interface=c3.C3dAffineTool(fsl2ras = True, itk_transform=True, reference_file='/home/pauli/Development/core_shell/data/CIT_Midbrain_Atlas/CIT168_T1w_700um_MNI.nii.gz'), name='fsl2ras_downsample_atlas_from_struct_to_func')
preproc.connect(downsample_atlas_t2_to_func, 'out_matrix_file', fsl2ras_downsample_atlas_from_struct_to_func, 'transform_file')
preproc.connect(downsample_atlas_t2, 'out_file', fsl2ras_downsample_atlas_from_struct_to_func, 'source_file')

# convert to itk format
fsl2ras_downsample_atlas_from_struct_to_func_inv = pe.Node(interface=c3.C3dAffineTool(fsl2ras = True, itk_transform=True, source_file='/home/pauli/Development/core_shell/data/CIT_Midbrain_Atlas/CIT168_T1w_700um_MNI.nii.gz'), name='fsl2ras_downsample_atlas_from_struct_to_func_inv')
preproc.connect(downsample_atlas_t2_to_func_inv, ('out_file', unlist), fsl2ras_downsample_atlas_from_struct_to_func_inv, 'transform_file')
preproc.connect(downsample_atlas_t2, 'out_file', fsl2ras_downsample_atlas_from_struct_to_func_inv, 'reference_file')

# co-register mbRef to anatomical
mbRef_to_t2 = pe.MapNode(interface=fsl.FLIRT(dof=6, output_type = "NIFTI_GZ"), name='mbRef_to_t2', iterfield=['in_file'])
preproc.connect(unwarp_mbRef, 'unwarped_file', mbRef_to_t2, 'in_file')
preproc.connect(maskT2, 'out_file', mbRef_to_t2, 'reference')

# invert xfm 
invert_t2_to_mbRef_tfm = pe.MapNode(interface=fsl.ConvertXFM(output_type = "NIFTI_GZ", invert_xfm = True), name='invert_t2_to_mbRef_tfm', iterfield=['in_file'])
preproc.connect(mbRef_to_t2, 'out_matrix_file', invert_t2_to_mbRef_tfm, 'in_file')

# transform t2 mask to func space
t2Mask_to_mbRef = pe.MapNode(interface=fsl.FLIRT(apply_xfm = True, output_type = "NIFTI_GZ", interp='nearestneighbour'), name='t2Mask_to_mbRef', iterfield=['reference','in_matrix_file'])
preproc.connect(inflate_mask, 'out_file', t2Mask_to_mbRef, 'in_file')
preproc.connect(unwarp_mbRef, 'unwarped_file', t2Mask_to_mbRef, 'reference')
preproc.connect(invert_t2_to_mbRef_tfm, 'out_file', t2Mask_to_mbRef, 'in_matrix_file')

# apply t2 mask to mbRef
mask_mbRef = pe.MapNode(interface=fsl.ImageMaths(op_string = '-mas', suffix = '_bet'), name='mask_mbRef', iterfield=['in_file','in_file2'])
preproc.connect(unwarp_mbRef, 'unwarped_file', mask_mbRef, 'in_file')
preproc.connect(t2Mask_to_mbRef, 'out_file', mask_mbRef, 'in_file2')

# convert to itk format
fsl2ras_mbRef_to_t2 = pe.MapNode(interface=c3.C3dAffineTool(fsl2ras = True, itk_transform=True), name='fsl2ras_mbRef_to_t2', iterfield=['source_file','transform_file'])
preproc.connect(maskT2, 'out_file', fsl2ras_mbRef_to_t2, 'reference_file')
preproc.connect(unwarp_mbRef, 'unwarped_file', fsl2ras_mbRef_to_t2, 'source_file')
preproc.connect(mbRef_to_t2, 'out_matrix_file', fsl2ras_mbRef_to_t2, 'transform_file')

# assemble inputs for mbref to t2 coregistration
ants_mbRef_to_t2_input = pe.MapNode(interface=util.IdentityInterface(fields=['fixed_image', 'moving_image', 'fixed_image_mask', 'moving_image_mask', 'initial_moving_transform']), name='ants_mbRef_to_t2_input', iterfield=['moving_image','moving_image_mask','initial_moving_transform'])
preproc.connect(mask_mbRef, 'out_file', ants_mbRef_to_t2_input, 'moving_image')
preproc.connect(maskT2, 'out_file', ants_mbRef_to_t2_input, 'fixed_image')
preproc.connect(t2Mask_to_mbRef, 'out_file', ants_mbRef_to_t2_input, 'moving_image_mask')
preproc.connect(inflate_mask, 'out_file', ants_mbRef_to_t2_input, 'fixed_image_mask')
preproc.connect(fsl2ras_mbRef_to_t2, 'itk_transform', ants_mbRef_to_t2_input, 'initial_moving_transform')

# perform mbRef to t2 coregistration
ants_mbRef_to_t2 = pe.MapNode(interface=ants.Registration(dimension=3, transforms=['Affine'], metric=[['Mattes']], radius_or_number_of_bins=[[32]], metric_weight = [[1]], transform_parameters=[[.01]], number_of_iterations=[[100]], convergence_threshold=[1.e-6], convergence_window_size=[10], smoothing_sigmas=[[2.0]], shrink_factors=[[2]], use_histogram_matching=True, interpolation='Linear', invert_initial_moving_transform=False, sampling_strategy = [['Random']], sampling_percentage = [[0.05]], verbose=False, output_warped_image = 'output_warped_image.nii.gz'), name='ants_mbRef_to_t2', iterfield=['moving_image','moving_image_mask','initial_moving_transform'])
preproc.connect(ants_mbRef_to_t2_input, 'fixed_image', ants_mbRef_to_t2, 'fixed_image')
preproc.connect(ants_mbRef_to_t2_input, 'moving_image', ants_mbRef_to_t2, 'moving_image')
preproc.connect(ants_mbRef_to_t2_input, 'fixed_image_mask', ants_mbRef_to_t2, 'fixed_image_mask')
preproc.connect(ants_mbRef_to_t2_input, 'moving_image_mask', ants_mbRef_to_t2, 'moving_image_mask')
preproc.connect(ants_mbRef_to_t2_input, 'initial_moving_transform', ants_mbRef_to_t2, 'initial_moving_transform')

# assemble output of this
ants_mbRef_to_t2_output = pe.MapNode(interface=util.IdentityInterface(fields=['forward_transforms', 'warped_image', 'reverse_transforms', 'inverse_warped_image']), name='ants_mbRef_to_t2_output', iterfield=['forward_transforms','warped_image','reverse_transforms','inverse_warped_image'])
preproc.connect(ants_mbRef_to_t2, 'forward_transforms', ants_mbRef_to_t2_output, 'forward_transforms')
preproc.connect(ants_mbRef_to_t2, 'warped_image', ants_mbRef_to_t2_output, 'warped_image')
preproc.connect(ants_mbRef_to_t2, 'reverse_transforms', ants_mbRef_to_t2_output, 'reverse_transforms')
preproc.connect(ants_mbRef_to_t2, 'inverse_warped_image', ants_mbRef_to_t2_output, 'inverse_warped_image')


# remove square brackets (['filename'] -> 'filename')
def get_first(mylist):
    """ get first element in a list
    
    Arguments:
        mylist {list} -- a list
    
    Returns:
        object -- first element in the list
    """
    r = mylist[0]
    return r

def get_sec(mylist):
        """ get second element in a list
    
    Arguments:
        mylist {list} -- a list
    
    Returns:
        object -- second element in the list
    """
    r = mylist[1]
    return r


def get_sec_tfm(in_list):
    """ get second transform from a list of transforms
    
    Arguments:
        in_list {list} -- list of transforms
    
    Returns:
        filename -- filename of transform
    """
    out_list = []
    for i in xrange(len(in_list)):
        out_list.append([in_list[i][0]])
    return out_list


sep_tfms = pe.Node(interface=Function(input_names=['in_list'], output_names=['out_list'], function=get_sec_tfm),name='sep_tfms')
preproc.connect(ants_mbRef_to_t2_output, 'forward_transforms', sep_tfms, 'in_list')

# combine transformations for (func to atlas)
cmb_tfm_func_to_ds_atlas = pe.MapNode(interface=util.Merge(4, axis='vstack'), name='cmb_tfm_func_to_ds_atlas', iterfield=['in4'])
# from atlas in struct resolution to func resolution
preproc.connect(fsl2ras_downsample_atlas_from_struct_to_func, 'itk_transform', cmb_tfm_func_to_ds_atlas, 'in1')
# from structural to atlas (with struct resolution)
preproc.connect(structs_to_atlas_coreg_output, ('forward_transforms', get_sec), cmb_tfm_func_to_ds_atlas, 'in2')
preproc.connect(structs_to_atlas_coreg_output, ('forward_transforms', get_first), cmb_tfm_func_to_ds_atlas, 'in3')
# transform from functional to structural
#preproc.connect(ants_mbRef_to_t2_output, ('forward_transforms', get_sec), cmb_tfm_func_to_ds_atlas, 'in4')
preproc.connect(sep_tfms, 'out_list', cmb_tfm_func_to_ds_atlas, 'in4')

# apply transformation from mbRef all the way to caltech atlas (in func resolution)
mbRef_to_ds_atlas = pe.MapNode(interface=ants.WarpImageMultiTransform(use_nearest = True), name = 'mbRef_to_ds_atlas', iterfield=['input_image', 'transformation_series'])
preproc.connect(mask_mbRef, 'out_file', mbRef_to_ds_atlas, 'input_image')
preproc.connect(downsample_atlas_t2_to_func, 'out_file', mbRef_to_ds_atlas, 'reference_image')
preproc.connect(cmb_tfm_func_to_ds_atlas, 'out', mbRef_to_ds_atlas, 'transformation_series')

# apply transformation to 4D func
func_to_ds_atlas = pe.MapNode(interface=ants.WarpTimeSeriesImageMultiTransform(use_nearest = True), name = 'func_to_ds_atlas', iterfield=['transformation_series','input_image'])
preproc.connect(unwarp_func, 'unwarped_file', func_to_ds_atlas, 'input_image')
preproc.connect(downsample_atlas_t2_to_func, 'out_file', func_to_ds_atlas, 'reference_image')
preproc.connect(cmb_tfm_func_to_ds_atlas, 'out', func_to_ds_atlas, 'transformation_series')

# combine transformations (t2 to atlas)
cmb_tfm_t2_to_ds_atlas = pe.Node(interface=util.Merge(3, axis='hstack'), name='cmb_tfm_t2_to_ds_atlas')
# from structural to atlas (with struct resolution)
preproc.connect(structs_to_atlas_coreg_output, ('forward_transforms', get_sec), cmb_tfm_t2_to_ds_atlas, 'in1')
preproc.connect(structs_to_atlas_coreg_output, ('forward_transforms', get_first), cmb_tfm_t2_to_ds_atlas, 'in2')
# from atlas in struct resolution to t2 resolution
preproc.connect(fsl2ras_downsample_atlas_from_struct_to_func, 'itk_transform', cmb_tfm_t2_to_ds_atlas, 'in3')

# apply transformation from t2 all the way to caltech atlas (in func resolution)
t2_to_ds_atlas = pe.Node(interface=ants.WarpImageMultiTransform(use_nearest = True), name = 't2_to_ds_atlas')
preproc.connect(maskT2, 'out_file', t2_to_ds_atlas, 'input_image')
preproc.connect(downsample_atlas_t2_to_func, 'out_file', t2_to_ds_atlas, 'reference_image')
preproc.connect(cmb_tfm_t2_to_ds_atlas, ('out', unlist), t2_to_ds_atlas, 'transformation_series')

# apply transformation from t1 all the way to caltech atlas (in func resolution)
t1_to_ds_atlas = pe.Node(interface=ants.WarpImageMultiTransform(use_nearest = True), name = 't1_to_ds_atlas')
preproc.connect(maskT1, 'out_file', t1_to_ds_atlas, 'input_image')
preproc.connect(downsample_atlas_t1_to_func, 'out_file', t1_to_ds_atlas, 'reference_image')
preproc.connect(cmb_tfm_t2_to_ds_atlas, ('out', unlist), t1_to_ds_atlas, 'transformation_series')

# mask downsampled multiband reference 
mask_ds_mbRef = pe.MapNode(interface=fsl.ImageMaths(op_string = '-thr 1 -bin'), name = 'mask_ds_mbRef', iterfield=['in_file'])
preproc.connect(mbRef_to_ds_atlas, 'output_image',  mask_ds_mbRef, 'in_file')

# assemble the normalized mbRef of all participants
join_mbRef = pe.JoinNode(interface=util.IdentityInterface(fields=['files']), joinfield=['files'], joinsource='datasource', name='join_mbRef')
preproc.connect(mask_ds_mbRef, 'out_file', join_mbRef, 'files')

# join normalized mbref images 
merge_mbRef = pe.Node(interface=fsl.Merge(dimension = 't'), name = 'merge_mbRef') 
preproc.connect(join_mbRef, ('files',  flatten), merge_mbRef, 'in_files')

# calculate how well normalized mbref overlap. 
calc_overlap_mbRef = pe.Node(interface=fsl.ImageMaths(op_string = '-Tmean'), name = 'calc_overlap_mbRef') 
preproc.connect(merge_mbRef, 'merged_file',  calc_overlap_mbRef, 'in_file')

# create a mask based on the mbRef 
mask_overlap_mbRef = pe.Node(interface=fsl.ImageMaths(op_string = '-thr .95 -bin'), name = 'mask_overlap_mbRef')
preproc.connect(calc_overlap_mbRef, 'out_file', mask_overlap_mbRef, 'in_file')

# determine ROI
get_roi = pe.Node(interface=fsl.ImageStats(op_string='-w'), name = 'get_roi') 
preproc.connect(mask_overlap_mbRef, 'out_file',  get_roi, 'in_file')

split_roi_coords = pe.Node(interface=util.Split(splits=[1,1,1,1,1,1,1,1]), name='split_roi_coords') 
preproc.connect(get_roi, 'out_stat',  split_roi_coords, 'inlist')


def unlist_long(mylist):
    """ returns first element of a list as long
    
    Arguments:
        mylist {list} -- list of numbers
    
    Returns:
        long -- first element of a list as long
    """

    r = mylist[0]
    return long(r)

# crop  mbRef
crop_mbRef = pe.MapNode(interface=fsl.ExtractROI(), name = 'crop_mbRef', iterfield=['in_file'])#,'x_min','x_size','y_min','y_size','z_min','z_size'])
preproc.connect(mbRef_to_ds_atlas, 'output_image',  crop_mbRef, 'in_file')
preproc.connect(split_roi_coords, ('out1', unlist_long), crop_mbRef, 'x_min')
preproc.connect(split_roi_coords, ('out2', unlist_long), crop_mbRef, 'x_size')
preproc.connect(split_roi_coords, ('out3', unlist_long), crop_mbRef, 'y_min')
preproc.connect(split_roi_coords, ('out4', unlist_long), crop_mbRef, 'y_size')
preproc.connect(split_roi_coords, ('out5', unlist_long), crop_mbRef, 'z_min')
preproc.connect(split_roi_coords, ('out6', unlist_long), crop_mbRef, 'z_size')

# crop functional
crop_func = pe.MapNode(interface=fsl.ExtractROI(), name = 'crop_func', iterfield=['in_file'])#,'x_min','x_size','y_min','y_size','z_min','z_size'])
preproc.connect(func_to_ds_atlas, 'output_image',  crop_func, 'in_file')
preproc.connect(split_roi_coords, ('out1', unlist_long), crop_func, 'x_min')
preproc.connect(split_roi_coords, ('out2', unlist_long), crop_func, 'x_size')
preproc.connect(split_roi_coords, ('out3', unlist_long), crop_func, 'y_min')
preproc.connect(split_roi_coords, ('out4', unlist_long), crop_func, 'y_size')
preproc.connect(split_roi_coords, ('out5', unlist_long), crop_func, 'z_min')
preproc.connect(split_roi_coords, ('out6', unlist_long), crop_func, 'z_size')

# crop functional and mbRef
crop_overlap_mask = pe.MapNode(interface=fsl.ExtractROI(), name = 'crop_overlap_mask', iterfield=['in_file'])#,'x_min','x_size','y_min','y_size','z_min','z_size'])
preproc.connect(mask_overlap_mbRef, 'out_file',  crop_overlap_mask, 'in_file')
preproc.connect(split_roi_coords, ('out1', unlist_long), crop_overlap_mask, 'x_min')
preproc.connect(split_roi_coords, ('out2', unlist_long), crop_overlap_mask, 'x_size')
preproc.connect(split_roi_coords, ('out3', unlist_long), crop_overlap_mask, 'y_min')
preproc.connect(split_roi_coords, ('out4', unlist_long), crop_overlap_mask, 'y_size')
preproc.connect(split_roi_coords, ('out5', unlist_long), crop_overlap_mask, 'z_min')
preproc.connect(split_roi_coords, ('out6', unlist_long), crop_overlap_mask, 'z_size')



# prepare mask

downsample_striatum_to_func = pe.Node(interface=fsl.FLIRT(output_type = "NIFTI_GZ", in_file=os.path.join(data_dir, 'openfmri/group/striatum_one_700mu_mni_bin.nii.gz'), reference=os.path.join(data_dir, 'openfmri/group/CIT168_T1w_700um_MNI.nii.gz'), in_matrix_file=os.path.join(fsl_dir, '/etc/flirtsch/ident.mat'), interp='nearestneighbour'), name='downsample_striatum_to_func')
preproc.connect(inputnode, ('func', get_voxel_size), downsample_striatum_to_func, 'apply_isoxfm')

# crop functional and mbRef
crop_striatal_mask = pe.Node(interface=fsl.ExtractROI(), name = 'crop_striatal_mask')
preproc.connect(downsample_striatum_to_func, 'out_file',  crop_striatal_mask, 'in_file')
preproc.connect(split_roi_coords, ('out1', unlist_long), crop_striatal_mask, 'x_min')
preproc.connect(split_roi_coords, ('out2', unlist_long), crop_striatal_mask, 'x_size')
preproc.connect(split_roi_coords, ('out3', unlist_long), crop_striatal_mask, 'y_min')
preproc.connect(split_roi_coords, ('out4', unlist_long), crop_striatal_mask, 'y_size')
preproc.connect(split_roi_coords, ('out5', unlist_long), crop_striatal_mask, 'z_min')
preproc.connect(split_roi_coords, ('out6', unlist_long), crop_striatal_mask, 'z_size')

# crop functional and mbRef
crop_t2 = pe.Node(interface=fsl.ExtractROI(), name = 'crop_t2')
preproc.connect(downsample_atlas_t2_to_func, 'out_file',  crop_t2, 'in_file')
preproc.connect(split_roi_coords, ('out1', unlist_long), crop_t2, 'x_min')
preproc.connect(split_roi_coords, ('out2', unlist_long), crop_t2, 'x_size')
preproc.connect(split_roi_coords, ('out3', unlist_long), crop_t2, 'y_min')
preproc.connect(split_roi_coords, ('out4', unlist_long), crop_t2, 'y_size')
preproc.connect(split_roi_coords, ('out5', unlist_long), crop_t2, 'z_min')
preproc.connect(split_roi_coords, ('out6', unlist_long), crop_t2, 'z_size')

# crop functional and mbRef
crop_t1 = pe.Node(interface=fsl.ExtractROI(), name = 'crop_t1')
preproc.connect(downsample_atlas_t1_to_func, 'out_file',  crop_t1, 'in_file')
preproc.connect(split_roi_coords, ('out1', unlist_long), crop_t1, 'x_min')
preproc.connect(split_roi_coords, ('out2', unlist_long), crop_t1, 'x_size')
preproc.connect(split_roi_coords, ('out3', unlist_long), crop_t1, 'y_min')
preproc.connect(split_roi_coords, ('out4', unlist_long), crop_t1, 'y_size')
preproc.connect(split_roi_coords, ('out5', unlist_long), crop_t1, 'z_min')
preproc.connect(split_roi_coords, ('out6', unlist_long), crop_t1, 'z_size')


# project downsampled ventricle mask into subject functional space
sep_inv_tfms = pe.Node(interface=Function(input_names=['in_list'], output_names=['out_list'], function=get_sec_tfm),name='sep_inv_tfms')
preproc.connect(ants_mbRef_to_t2_output, 'reverse_transforms', sep_inv_tfms, 'in_list')

# combine transformations for (func to atlas)
cmb_tfm_ds_atlas_to_func = pe.MapNode(interface=util.Merge(4, axis='vstack'), name='cmb_tfm_ds_atlas_to_func', iterfield=['in1'])
# from atlas in struct resolution to func resolution
preproc.connect(fsl2ras_downsample_atlas_from_struct_to_func_inv, 'itk_transform', cmb_tfm_ds_atlas_to_func, 'in4')
# from structural to atlas (with struct resolution)
preproc.connect(structs_to_atlas_coreg_output, ('reverse_transforms', get_sec), cmb_tfm_ds_atlas_to_func, 'in3')
preproc.connect(structs_to_atlas_coreg_output, ('reverse_transforms', get_first), cmb_tfm_ds_atlas_to_func, 'in2')
# transform from functional to structural
#preproc.connect(ants_mbRef_to_t2_output, ('forward_transforms', get_sec), cmb_tfm_ds_atlas_to_func, 'in4')
preproc.connect(sep_inv_tfms, 'out_list', cmb_tfm_ds_atlas_to_func, 'in1')




# apply transformation to 4D func
ds_ventricle_mask_to_func = pe.MapNode(interface=ants.WarpImageMultiTransform(use_nearest = True, input_image=ds_ventricle_mask), name = 'ds_ventricle_mask_to_func', iterfield=['transformation_series', 'reference_image'])
preproc.connect(meanfunc, 'out_file', ds_ventricle_mask_to_func, 'reference_image')
preproc.connect(cmb_tfm_ds_atlas_to_func, 'out', ds_ventricle_mask_to_func, 'transformation_series')






# ---------------- RUN! (good luck)---------------------

preproc.write_graph()
outgraph = preproc.run(plugin='MultiProc', plugin_args={'n_procs' : 5})

