# Copyright 2023 Population Health Sciences, German Center for Neurodegenerative Diseases (DZNE)
#
#    Licensed under the Apache License, Version 2.0 (the "License");
#    you may not use this file except in compliance with the License.
#    You may obtain a copy of the License at
#
#        http://www.apache.org/licenses/LICENSE-2.0
#
#    Unless required by applicable law or agreed to in writing, software
#    distributed under the License is distributed on an "AS IS" BASIS,
#    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#    See the License for the specific language governing permissions and
#    limitations under the License.


def file_selector(scans_dir,subject_id):
    import os
    import glob
    from nipype.interfaces.fsl import ExtractROI

    restfile=None
    strucfile=None
    magfile=None
    phasefile=None

    if os.path.exists(os.path.join(scans_dir,subject_id,'RestingState.nii.gz')):
        restfile=os.path.abspath(os.path.join(scans_dir,subject_id,'RestingState.nii.gz'))
    else:
        raise IOError("No RestingState.nii.gz file found.")

    if os.path.exists(os.path.join(scans_dir,subject_id,'B0.nii.gz')):
        magfile=os.path.abspath(os.path.join(scans_dir,subject_id,'B0.nii.gz'))
        exroi=ExtractROI(in_file=magfile, roi_file='B0.nii.gz', t_min=0,t_size=1)
        res=exroi.run()
        magfile=res.outputs.get()['roi_file'] 
    elif os.path.exists(os.path.join(scans_dir,subject_id,'Scout.nii.gz')):
        magfile=os.path.abspath(os.path.join(scans_dir,subject_id,'Scout.nii.gz'))
    else:
        raise IOError("No Magnitude B0 or alternative Scout file found.")

    strucfiles= [t for t in glob.glob(os.path.join(scans_dir,subject_id,'T1*.gz')) if os.path.isfile(t)]

    if strucfiles:
        strucfile=strucfiles[0]
    elif os.path.exists(os.path.join(scans_dir,subject_id,'Scout.nii.gz')):
        strucfile=os.path.abspath(os.path.join(scans_dir,subject_id,'Scout.nii.gz'))
    else:
        raise IOError("No structural file T1*.nii.gz or Scout.nii.gz found.")
    
    if os.path.exists(os.path.join(scans_dir,subject_id,'B0_Phase.nii.gz')):
        phasefile=os.path.abspath(os.path.join(scans_dir,subject_id,'B0_Phase.nii.gz'))     
    else:
        raise IOError("No Phase image found.")

    return restfile,strucfile,magfile,phasefile


def prepare_fieldmap(name='prepare_fieldmap'):
    """Return a field map creation workflow

    Based on FSL's bet,flirt and fsl_prepare_fieldmap.
    """
    import nipype.pipeline.engine as pe
    import nipype.interfaces.fsl as fsl
    from nipype.interfaces.fsl.maths import ErodeImage,BinaryMaths
    from nipype import IdentityInterface

    fieldmapwf = pe.Workflow(name=name)

    inputnode = pe.Node(interface=IdentityInterface(fields=['magfile', 'phasefile']), name='inputspec')

    flirter = pe.Node(interface=fsl.FLIRT(),name='flirt')
    flirter.inputs.apply_xfm=True
    flirter.inputs.uses_qform=True


    better = pe.Node(interface=fsl.BET(), name='bet')
    better.inputs.mask = False
    better.inputs.frac=0.65
    better.inputs.robust=True

    eroder=pe.Node(interface=ErodeImage(),name='erode')

    fieldmap = pe.Node(interface=fsl.PrepareFieldmap(),name='fieldmap')
    fieldmap.inputs.scanner='SIEMENS'
    fieldmap.inputs.delta_TE=2.46

    fieldmap2hz = pe.Node(interface=BinaryMaths(),name='fieldmap2hz')
    fieldmap2hz.inputs.operand_value=0.159155
    fieldmap2hz.inputs.operation='mul'
    fieldmap2hz.inputs.out_file='fmap_hz.nii.gz'


    outputnode = pe.Node(interface=IdentityInterface(fields=['fieldmapfile']), name='outputspec')

    fieldmapwf.connect(inputnode   , 'magfile'     , flirter     , 'in_file')
    fieldmapwf.connect(inputnode   , 'phasefile'   , flirter     , 'reference')
    fieldmapwf.connect(flirter     , 'out_file'    , better      , 'in_file')
    fieldmapwf.connect(better      , 'out_file'    , eroder      , 'in_file')
    fieldmapwf.connect(eroder      , 'out_file'    , fieldmap    , 'in_magnitude')
    fieldmapwf.connect(inputnode   , 'phasefile'   , fieldmap    , 'in_phase')
    fieldmapwf.connect(fieldmap    , 'out_fieldmap',fieldmap2hz  , 'in_file')
    fieldmapwf.connect(fieldmap2hz , 'out_file'    , outputnode  , 'fieldmapfile')

    return fieldmapwf

def _mel_ica_dir(in_ica_dir,in_filtered_func,in_hires,mean_func,in_brainmask,
                in_motion_file, roi_file, out_matrix_file):

    import os
    from nipype.utils.filemanip import copyfile
    from nipype.interfaces.io import copytree

    out_ica_dir=os.path.abspath(os.getcwd())
    
    regdir=os.path.join(out_ica_dir,'reg')
    mcdir =os.path.join(out_ica_dir,'mc')

    os.makedirs(regdir,exist_ok=True)
    os.makedirs(mcdir,exist_ok=True)

    copytree(in_ica_dir, os.path.join(os.getcwd(),os.path.basename(in_ica_dir)), use_hardlink=True)

    copyfile( in_motion_file, os.path.join(mcdir,'prefiltered_func_data_mcf.par'),use_hardlink=True)
    copyfile( in_hires, os.path.join(regdir, 'highres.nii.gz'), use_hardlink=True)
    copyfile( roi_file, os.path.join(regdir, 'example_func.nii.gz'), use_hardlink=True)
    copyfile( out_matrix_file, os.path.join(regdir, 'highres2example_func.mat'), use_hardlink=True)
    copyfile( in_brainmask, os.path.join(os.getcwd(),'brain_mask.nii.gz'), use_hardlink=True)
    copyfile( mean_func, os.path.join(os.getcwd(),'mean_func.nii.gz'), use_hardlink=True)
    copyfile( in_filtered_func, os.path.join(os.getcwd(),'filtered_func_data.nii.gz'), use_hardlink=True)
    copyfile( os.path.join(in_ica_dir,'mask.nii.gz'), os.path.join(os.getcwd(),'mask.nii.gz'),use_hardlink=True)

    return out_ica_dir


def prepare_fix_wf(name='fix_workflow'):

    import nipype.pipeline.engine as pe
    import nipype.interfaces.fsl as fsl
    from nipype import IdentityInterface
    from nipype.interfaces.utility import Function

    fixwf = pe.Workflow(name=name)

    inputnode = pe.Node(interface=IdentityInterface(fields=['in_filtered_func', 'in_hires']),
                                                            name='inputspec')

    #fslroi _hp.ica/filtered_func_data _hp.ica/reg/example_func 526 1
    fslroi = pe.Node(interface=fsl.ExtractROI(),name='extract_roi')
    fslroi.inputs.roi_file='example_func.nii.gz'
    fslroi.inputs.t_min=526
    fslroi.inputs.t_size=1

    hires2func = pe.Node(interface=fsl.FLIRT(),name='hires2func')
    hires2func.inputs.out_matrix_file='highres2example_func.mat'

    outputnode = pe.Node(interface=IdentityInterface(fields=['roi_file','out_matrix_file']), name='outputspec')

    fixwf.connect(inputnode      , 'in_filtered_func',     fslroi             , 'in_file')
    fixwf.connect(inputnode      , 'in_hires',             hires2func         , 'in_file')
    fixwf.connect(fslroi         , 'roi_file',             hires2func         , 'reference')
    fixwf.connect(fslroi         , 'roi_file',             outputnode         , 'roi_file')
    fixwf.connect(hires2func     , 'out_matrix_file',      outputnode         , 'out_matrix_file')

    return fixwf
    

