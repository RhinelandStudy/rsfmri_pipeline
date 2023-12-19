# -*- coding: utf-8 -*-

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


"""A nipype pre-processing pipeline for timeseries MR data (fMRI/EPI).

Created on Wed Mar 17 2021

@authors: hweeling-lee, shahidm

Copyright 2023 Population Health Sciences, German Center for Neurodegenerative Diseases (DZNE)
Licensed under the Apache License, Version 2.0 (the "License"); 
you may not use this file except in compliance with the License. 
You may obtain a copy of the License at  http://www.apache.org/licenses/LICENSE-2.0 
Unless required by applicable law or agreed to in writing, software distributed under the 
License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and limitations under the License.

"""

#%%
# Imports
from __future__ import division
import nipype.pipeline.engine as pe
import argparse
import sys
import os
import glob
from itertools import chain
from nipype import config, logging


from .spm_interface import SPMPreprocess, SPMGLMProcess, SPMFCMeans, SPMGLMProcessST, SPMGLMProcessVS, SPMGLMProcessGordon
from .utils import file_selector,prepare_fieldmap, prepare_fix_wf,_mel_ica_dir
from .custom_fix_interface import FixCleaner
from .fmri_qc import fmri_qc_workflow

import nipype.interfaces.fsl as fsl
from nipype.interfaces.fsl.model import MELODIC

from nipype.interfaces.fsl.maths import MeanImage,MultiImageMaths,ChangeDataType

from nipype import IdentityInterface, DataSink
from nipype.interfaces.utility import Function, Merge



def create_preprocess_wf(scans_dir, work_dir, subj_ids, traindata,fixth, mthreads,spm_exe_script, mcr_path, name='preprocess_workflow'):

    ppwf = pe.Workflow(name=name)
    ppwf.base_dir=work_dir


    inputnode = pe.Node(interface=IdentityInterface(fields=['subj_ids', 'outputdir']), name='inputnode')
    inputnode.iterables = [('subj_ids', subj_ids)]

    fileselector = pe.Node(interface=Function(input_names=['scans_dir','subject_id'],
                                     output_names=['restfile','strucfile','magfile','phasefile'],
                                     function=file_selector),name='file_selector')
    fileselector.inputs.scans_dir=scans_dir

    # prepare fieldmap  
    prepfieldmap = prepare_fieldmap(name='prepare_fieldmap')


    #SPM preproc 
    spmpreproc=pe.Node(interface=SPMPreprocess(), name='spmpreproc')
    spmpreproc.n_procs=mthreads
    spmpreproc.inputs.TR='0.57'
    spmpreproc.inputs.TE='0.3'
    spmpreproc.inputs.ndummy='17'
    spmpreproc.inputs.spm_exe_script = spm_exe_script
    spmpreproc.inputs.mcr_path = mcr_path
    spmpreproc.terminal_output='file'

    #%###################
    #preFIX ica
    #####################

    #tmean
    mean_func = pe.Node(interface=MeanImage(), name='mean_func')
    mean_func.inputs.dimension='T'
    mean_func.inputs.out_file='mean_func.nii.gz'

    #get bet brain mask of tmean
    bet_mean_func=pe.Node(interface=fsl.BET(), name='bet_mean_func')
    bet_mean_func.inputs.mask=True
    bet_mean_func.inputs.frac=0.4
    bet_mean_func.inputs.vertical_gradient=0
    bet_mean_func.inputs.threshold=True
    bet_mean_func.inputs.no_output=True

    #get list of input files for filtering using fslmaths
    merger=pe.Node(interface=Merge(2), name='list_mean_func')

    #fslmaths ${fname} -sub mean_func -bptf 187.1 -1 -add mean_func filtered_func_data
    filter_func = pe.Node(interface=MultiImageMaths(), name='filter_func')
    filter_func.inputs.op_string = "-sub %s -bptf 187.1 -1 -add %s"
    filter_func.inputs.out_file='filtered_func_data.nii.gz'

    #melodic -i _hp.ica/filtered_func_data -o _hp.ica/filtered_func_data.ica -m _hp.ica/brain_mask.nii.gz --tr=0.57 -a concat -d -250 --nobet --Oall
    melodic = pe.Node(interface=MELODIC(), name='melodic')
    melodic.inputs.no_bet=True
    melodic.inputs.out_all=True
    melodic.inputs.approach='concat'
    melodic.inputs.dim=250
    melodic.inputs.tr_sec=0.57
    melodic.inputs.out_dir='filtered_func_data.ica'

    #reg hires to avg filtered func
    prefixwf = prepare_fix_wf(name='pre_fix_workflow')

    mel_ica_dir = pe.Node(interface=Function(input_names=['in_ica_dir','in_filtered_func','in_hires','mean_func',
                                                          'in_brainmask','in_motion_file',
                                                          'roi_file', 'out_matrix_file'],
                                             output_names=['out_ica_dir'],
                                             function=_mel_ica_dir),name='mel_ica_dir')

    #%# FIX cleaner
    cleaner = pe.Node(interface=FixCleaner(), name='fix_cleaner')
    cleaner.inputs.thresh=fixth
    cleaner.inputs.trained_wts_file=traindata

    #change dtype from FLOAT to short as SPM has problems converting 4d to 3d
    chdt = pe.Node(interface=ChangeDataType(), name='chdt')
    chdt.inputs.output_datatype='short'

    #inject here fmri qc wf
    #updated 27.06.2021 qc part in matlab scripts
    fmri_qc_wf=fmri_qc_workflow(spm_exe_script,mcr_path)
    
    #SPM 17-7 Networks
    spmglmproc = pe.Node(interface=SPMGLMProcess(), name='spmglmproc')
    spmglmproc.n_procs=mthreads
    spmglmproc.inputs.spm_exe_script = spm_exe_script
    spmglmproc.inputs.mcr_path = mcr_path
    spmglmproc.terminal_output='file'

    #get FC means for FC networks
    spmfcmeans = pe.Node(interface=SPMFCMeans(), name='spmfcmeans')
    spmfcmeans.n_procs=mthreads
    spmfcmeans.inputs.spm_exe_script = spm_exe_script
    spmfcmeans.inputs.mcr_path = mcr_path
    spmfcmeans.terminal_output='file'



    #SPM GLM proc 
    #NN# spmglmproc = pe.Node(interface=SPMGLMProcess(), name='spmglmproc')
    #NN# spmglmproc.n_procs=mthreads
    #NN# spmglmproc.inputs.spm_exe_script = spm_exe_script
    #NN# spmglmproc.inputs.mcr_path = mcr_path
    #NN# spmglmproc.terminal_output='file'

    #SPM GLM proc 200_100
    #NN# spmglmproc_st = pe.Node(interface=SPMGLMProcessST(), name='spmglmproc_st')
    #NN# spmglmproc_st.n_procs=mthreads
    #NN# spmglmproc_st.inputs.spm_exe_script = spm_exe_script
    #NN# spmglmproc_st.inputs.mcr_path = mcr_path
    #NN# spmglmproc_st.terminal_output='file'


    #SPM GLM proc 
    #NN# spmglmproc_vs = pe.Node(interface=SPMGLMProcessVS(), name='spmglmproc_vs')
    #NN# spmglmproc_vs.n_procs=mthreads
    #NN# spmglmproc_vs.inputs.spm_exe_script = spm_exe_script
    #NN# spmglmproc_vs.inputs.mcr_path = mcr_path
    #NN# spmglmproc_vs.terminal_output='file'

    #SPM GLM proc Gordon atlas
    spmglmproc_gordon = pe.Node(interface=SPMGLMProcessGordon(), name='spmglmproc_gordon')
    spmglmproc_gordon.n_procs=mthreads
    spmglmproc_gordon.inputs.spm_exe_script = spm_exe_script
    spmglmproc_gordon.inputs.mcr_path = mcr_path
    spmglmproc_gordon.terminal_output='file'

 
    #%# datasinker
    datasink = pe.Node(interface=DataSink(), name='datasinker')
    datasink.inputs.parameterization=False
    datasink.inputs.regexp_substitutions = [(r'/[0-9]+_(.*).nii.gz', r'/\1.nii.gz'),
                                            (r'/[0-9]+_(.*).nii', r'/\1.nii'),
                                            (r'/[0-9]+_(.*).mat', r'/\1.mat'),
                                            (r'/[0-9]+_(.*).csv', r'/\1.csv')]

    #%# workflow connections
    ppwf.connect(inputnode       , 'subj_ids'     , fileselector  , 'subject_id')
    ppwf.connect(fileselector    , 'magfile'      , prepfieldmap  , 'inputspec.magfile')
    ppwf.connect(fileselector    , 'phasefile'    , prepfieldmap  , 'inputspec.phasefile')

    ppwf.connect(prepfieldmap    , 'outputspec.fieldmapfile' , spmpreproc   , 'fieldmapfile')
    ppwf.connect(fileselector    , 'magfile'      , spmpreproc   , 'magfile')
    ppwf.connect(fileselector    , 'phasefile'    , spmpreproc   , 'phasefile')
    ppwf.connect(fileselector    , 'strucfile'    , spmpreproc   , 'strucfile')
    ppwf.connect(fileselector    , 'restfile'     , spmpreproc   , 'restfile')

    ppwf.connect(spmpreproc     , 'swurestfile'  , mean_func     , 'in_file')
    ppwf.connect(mean_func      , 'out_file'     , bet_mean_func , 'in_file')

    ppwf.connect(mean_func      , 'out_file'     , merger        , 'in1')
    ppwf.connect(mean_func      , 'out_file'     , merger        , 'in2')

    ppwf.connect(spmpreproc     , 'swurestfile'  , filter_func   , 'in_file')
    ppwf.connect(merger         , 'out'          , filter_func   , 'operand_files')

    ppwf.connect(bet_mean_func  , 'mask_file'    , melodic       , 'mask')
    ppwf.connect(filter_func    , 'out_file'     , melodic       , 'in_files')

    ppwf.connect(filter_func    , 'out_file'     , prefixwf      , 'inputspec.in_filtered_func')
    ppwf.connect(spmpreproc     , 'wt1hiresfile' , prefixwf      , 'inputspec.in_hires')

    ppwf.connect(melodic        , 'out_dir'                     , mel_ica_dir       , 'in_ica_dir')
    ppwf.connect(filter_func    , 'out_file'                    , mel_ica_dir       , 'in_filtered_func')
    ppwf.connect(spmpreproc     , 'wt1hiresfile'                , mel_ica_dir       , 'in_hires')
    ppwf.connect(mean_func      , 'out_file'                    , mel_ica_dir       , 'mean_func')
    ppwf.connect(bet_mean_func  , 'mask_file'                   , mel_ica_dir       , 'in_brainmask')
    ppwf.connect(spmpreproc     , 'fslmotionfile'               , mel_ica_dir       , 'in_motion_file')
    ppwf.connect(prefixwf       , 'outputspec.roi_file'         , mel_ica_dir       , 'roi_file')
    ppwf.connect(prefixwf       , 'outputspec.out_matrix_file'  , mel_ica_dir       , 'out_matrix_file')

    ppwf.connect(mel_ica_dir    , 'out_ica_dir'                 , cleaner       , 'mel_ica')
    ppwf.connect(cleaner        , 'cleaned_functional_file'     , chdt          , 'in_file')

    ppwf.connect(spmpreproc     , 'swurestfile'       , fmri_qc_wf    , 'inputnode.swufileb4ica')
    ppwf.connect(chdt           , 'out_file'          , fmri_qc_wf    , 'inputnode.swufileafterica')
    ppwf.connect(spmpreproc     , 'wc2t1file'         , fmri_qc_wf    , 'inputnode.whitematterfile')
    ppwf.connect(spmpreproc     , 'wc1t1file'         , fmri_qc_wf    , 'inputnode.graymatterfile')
    ppwf.connect(spmpreproc     , 'wc3t1file'         , fmri_qc_wf    , 'inputnode.csffile')
    ppwf.connect(spmpreproc     , 'fslmotionfile'     , fmri_qc_wf    , 'inputnode.rpfile')
    ppwf.connect(spmpreproc     , 't1segmat'          , fmri_qc_wf    , 'inputnode.segmentfile')


    ppwf.connect(chdt           , 'out_file'                    , spmglmproc    , 'swurestfile') 
    ppwf.connect(spmpreproc     , 'wc2t1file'                   , spmglmproc    , 'wc2t1file')
    ppwf.connect(spmpreproc     , 'wc3t1file'                   , spmglmproc    , 'wc3t1file')
    ppwf.connect(spmpreproc     , 'fslmotionfile'               , spmglmproc    , 'fslmotionfile')

    #NN
    #ppwf.connect(chdt           , 'out_file'                    , spmglmproc_st , 'swurestfile')
    #ppwf.connect(spmpreproc     , 'wc2t1file'                   , spmglmproc_st , 'wc2t1file')
    #ppwf.connect(spmpreproc     , 'wc3t1file'                   , spmglmproc_st , 'wc3t1file')
    #ppwf.connect(spmpreproc     , 'fslmotionfile'               , spmglmproc_st , 'fslmotionfile')
 
    #NN
    #ppwf.connect(chdt           , 'out_file'                    , spmglmproc_vs , 'swurestfile')
    #ppwf.connect(spmpreproc     , 'wc2t1file'                   , spmglmproc_vs , 'wc2t1file')
    #ppwf.connect(spmpreproc     , 'wc3t1file'                   , spmglmproc_vs , 'wc3t1file')
    #ppwf.connect(spmpreproc     , 'fslmotionfile'               , spmglmproc_vs , 'fslmotionfile')

    ppwf.connect(spmglmproc     , 'fc_17full'                   , spmfcmeans , 'fc_17full')
    ppwf.connect(spmglmproc     , 'fc_7full'                    , spmfcmeans , 'fc_7full')

    ppwf.connect(chdt           , 'out_file'                    , spmglmproc_gordon , 'swurestfile')
    ppwf.connect(spmpreproc     , 'wc2t1file'                   , spmglmproc_gordon , 'wc2t1file')
    ppwf.connect(spmpreproc     , 'wc3t1file'                   , spmglmproc_gordon , 'wc3t1file')
    ppwf.connect(spmpreproc     , 'fslmotionfile'               , spmglmproc_gordon , 'fslmotionfile')


    #outputs
    ppwf.connect(inputnode       , 'subj_ids'       , datasink, 'container')
    ppwf.connect(inputnode       , 'outputdir'      , datasink, 'base_directory')

    #ppwf.connect(spmpreproc      , 'swurestfile'    , datasink, 'rs_preproc.@swrestfile')
    ppwf.connect(spmpreproc      , 'motionfile'     , datasink, 'rs_preproc.@motionfile')
    ppwf.connect(spmpreproc      , 'fslmotionfile'  , datasink, 'rs_preproc.@fslmotionfile')
    ppwf.connect(spmpreproc      , 'qcfile'         , datasink, 'rs_preproc.@qcfile')

    ppwf.connect(spmpreproc      , 'wt1hiresfile'   , datasink, 't1.@hires')

    ppwf.connect(spmpreproc      , 'wc1t1file'   , datasink, 't1.@wc1t1file')
    ppwf.connect(spmpreproc      , 'wc2t1file'   , datasink, 't1.@wc2t1file')
    ppwf.connect(spmpreproc      , 'wc3t1file'   , datasink, 't1.@wc3t1file')
    ppwf.connect(spmpreproc      , 't1segmat'    , datasink, 't1.@t1segmat')

    ppwf.connect(chdt            , 'out_file'    , datasink, 'ica_fix.@cleaned')
    ppwf.connect(cleaner         , 'fix4melview' , datasink, 'ica_fix.@fix4melview')

    ppwf.connect(spmglmproc      , 'fc_17raw'   , datasink, 'fc.@fc_17raw')
    ppwf.connect(spmglmproc      , 'fc_17full'  , datasink, 'fc.@fc_17full')
    ppwf.connect(spmglmproc      , 'fc_7raw'    , datasink, 'fc.@fc_7raw')
    ppwf.connect(spmglmproc      , 'fc_7full'   , datasink, 'fc.@fc_7full')

    ppwf.connect(spmfcmeans      , 'fc17meanfc' , datasink, 'fc.@fc17meanfc')
    ppwf.connect(spmfcmeans      , 'fc7meanfc'  , datasink, 'fc.@fc7meanfc')


    ##NN
    #ppwf.connect(spmglmproc_st   , 'fc_conn_file_200p'   , datasink, 'fc.@fc_conn_file_200p')
    #ppwf.connect(spmglmproc_st   , 'fc_conn_file_100p'   , datasink, 'fc.@fc_conn_file_100p')

    #NN
    #ppwf.connect(spmglmproc_vs   , 'vs_conn_file'   , datasink, 'con.vstriatum')
    #ppwf.connect(spmglmproc_vs   , 'lvs_conn_file'  , datasink, 'con.l_vstriatum')
    #ppwf.connect(spmglmproc_vs   , 'rvs_conn_file'  , datasink, 'con.r_vstriatum')

    ppwf.connect(spmglmproc_gordon, 'fc_conn_file'  , datasink, 'fc.@gordon_fc_conn_file')

    ppwf.connect(fmri_qc_wf      , 'outputnode.fdcarpb4ica'      , datasink, 'qc.@fdcarpb4icafile')
    ppwf.connect(fmri_qc_wf      , 'outputnode.fdcarpafterica'   , datasink, 'qc.@fdcarpaftericafile')
    ppwf.connect(fmri_qc_wf      , 'outputnode.allmetrics'       , datasink, 'qc.@allmetricsfile')


    return ppwf


def create_preprocess_wf_nocln(scans_dir, work_dir, subj_ids, traindata,fixth, mthreads,spm_exe_script, mcr_path, name='preprocess_workflow'):

    ppwf = pe.Workflow(name=name)
    ppwf.base_dir=work_dir

    inputnode = pe.Node(interface=IdentityInterface(fields=['subj_ids', 'outputdir']), name='inputnode')
    inputnode.iterables = [('subj_ids', subj_ids)]

    fileselector = pe.Node(interface=Function(input_names=['scans_dir','subject_id'],
                                     output_names=['restfile','strucfile','magfile','phasefile'],
                                     function=file_selector),name='file_selector')
    fileselector.inputs.scans_dir=scans_dir

    # prepare fieldmap  
    prepfieldmap = prepare_fieldmap(name='prepare_fieldmap')

    #SPM preproc 
    spmpreproc=pe.Node(interface=SPMPreprocess(), name='spmpreproc')
    spmpreproc.n_procs=mthreads
    spmpreproc.inputs.TR='0.57'
    spmpreproc.inputs.TE='0.3'
    spmpreproc.inputs.ndummy='17'
    spmpreproc.inputs.spm_exe_script = spm_exe_script
    spmpreproc.inputs.mcr_path = mcr_path
    spmpreproc.terminal_output='file'

    #tmean
    mean_func = pe.Node(interface=MeanImage(), name='mean_func')
    mean_func.inputs.dimension='T'
    mean_func.inputs.out_file='mean_func.nii.gz'

    #get list of input files for filtering using fslmaths
    merger=pe.Node(interface=Merge(2), name='list_mean_func')

    #fslmaths ${fname} -sub mean_func -bptf 187.1 -1 -add mean_func filtered_func_data
    filter_func = pe.Node(interface=MultiImageMaths(), name='filter_func')
    filter_func.inputs.op_string = "-sub %s -bptf 187.1 -1 -add %s"
    filter_func.inputs.out_file='filtered_func_data.nii.gz'

    #change dtype from FLOAT to short as SPM has problems converting 4d to 3d
    chdt = pe.Node(interface=ChangeDataType(), name='chdt')
    chdt.inputs.output_datatype='short'


    #SPM GLM proc
    spmglmproc = pe.Node(interface=SPMGLMProcess(), name='spmglmproc')
    spmglmproc.n_procs=mthreads
    spmglmproc.inputs.spm_exe_script = spm_exe_script
    spmglmproc.inputs.mcr_path = mcr_path
    spmglmproc.terminal_output='file'

    #SPM GLM proc
    spmglmproc_vs = pe.Node(interface=SPMGLMProcessVS(), name='spmglmproc_vs')
    spmglmproc_vs.n_procs=mthreads
    spmglmproc_vs.inputs.spm_exe_script = spm_exe_script
    spmglmproc_vs.inputs.mcr_path = mcr_path
    spmglmproc_vs.terminal_output='file'

    #SPM GLM proc Gordon atlas
    spmglmproc_gordon = pe.Node(interface=SPMGLMProcessGordon(), name='spmglmproc_gordon')
    spmglmproc_gordon.n_procs=mthreads
    spmglmproc_gordon.inputs.spm_exe_script = spm_exe_script
    spmglmproc_gordon.inputs.mcr_path = mcr_path
    spmglmproc_gordon.terminal_output='file'

    #%# datasinker
    datasink = pe.Node(interface=DataSink(), name='datasinker')
    datasink.inputs.parameterization=False
    datasink.inputs.regexp_substitutions = [(r'/[0-9]+_(.*).nii.gz', r'/\1.nii.gz'),
                                            (r'/[0-9]+_(.*).nii', r'/\1.nii'),
                                            (r'/[0-9]+_(.*).mat', r'/\1.mat'),
                                            (r'/[0-9]+_(.*).csv', r'/\1.csv')]

    #%# workflow connections
    ppwf.connect(inputnode       , 'subj_ids'     , fileselector  , 'subject_id')
    ppwf.connect(fileselector    , 'magfile'      , prepfieldmap  , 'inputspec.magfile')
    ppwf.connect(fileselector    , 'phasefile'    , prepfieldmap  , 'inputspec.phasefile')

    ppwf.connect(prepfieldmap    , 'outputspec.fieldmapfile' , spmpreproc   , 'fieldmapfile')
    ppwf.connect(fileselector    , 'magfile'      , spmpreproc   , 'magfile')
    ppwf.connect(fileselector    , 'phasefile'    , spmpreproc   , 'phasefile')
    ppwf.connect(fileselector    , 'strucfile'    , spmpreproc   , 'strucfile')
    ppwf.connect(fileselector    , 'restfile'     , spmpreproc   , 'restfile')

    ppwf.connect(spmpreproc     , 'swurestfile'   , mean_func     , 'in_file')

    ppwf.connect(mean_func      , 'out_file'      , merger        , 'in1')
    ppwf.connect(mean_func      , 'out_file'      , merger        , 'in2')

    ppwf.connect(spmpreproc     , 'swurestfile'   , filter_func   , 'in_file')
    ppwf.connect(merger         , 'out'           , filter_func   , 'operand_files')

    ppwf.connect(filter_func    , 'out_file'      , chdt          , 'in_file')


    ppwf.connect(chdt           , 'out_file'      , spmglmproc    , 'swurestfile')
    ppwf.connect(spmpreproc     , 'wc2t1file'     , spmglmproc    , 'wc2t1file')
    ppwf.connect(spmpreproc     , 'wc3t1file'     , spmglmproc    , 'wc3t1file')
    ppwf.connect(spmpreproc     , 'fslmotionfile' , spmglmproc    , 'fslmotionfile')

    ppwf.connect(chdt           , 'out_file'      , spmglmproc_vs , 'swurestfile')
    ppwf.connect(spmpreproc     , 'wc2t1file'     , spmglmproc_vs , 'wc2t1file')
    ppwf.connect(spmpreproc     , 'wc3t1file'     , spmglmproc_vs , 'wc3t1file')
    ppwf.connect(spmpreproc     , 'fslmotionfile' , spmglmproc_vs , 'fslmotionfile')

    ppwf.connect(chdt           , 'out_file'      , spmglmproc_gordon , 'swurestfile')
    ppwf.connect(spmpreproc     , 'wc2t1file'     , spmglmproc_gordon , 'wc2t1file')
    ppwf.connect(spmpreproc     , 'wc3t1file'     , spmglmproc_gordon , 'wc3t1file')
    ppwf.connect(spmpreproc     , 'fslmotionfile' , spmglmproc_gordon , 'fslmotionfile')

    #outputs
    ppwf.connect(inputnode       , 'subj_ids'       , datasink, 'container')
    ppwf.connect(inputnode       , 'outputdir'      , datasink, 'base_directory')

    ppwf.connect(spmglmproc      , 'fc_conn_file'   , datasink, 'fc.no_ica_fix.@fc_conn_file')

    ppwf.connect(spmglmproc_vs   , 'vs_conn_file'   , datasink, 'fc.no_ica_fix.vstriatum')
    ppwf.connect(spmglmproc_vs   , 'lvs_conn_file'  , datasink, 'fc.no_ica_fix.l_vstriatum')
    ppwf.connect(spmglmproc_vs   , 'rvs_conn_file'  , datasink, 'fc.no_ica_fix.r_vstriatum')

    ppwf.connect(spmglmproc_gordon, 'fc_conn_file'   , datasink, 'fc.@gordon_fc_conn_file')

    return ppwf

