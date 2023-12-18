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
from nipype import config, logging, SelectFiles

from .spm_interface import VBMPreprocess,VBMETIV

from nipype import IdentityInterface, DataSink


def create_preprocess_wf(scans_dir, work_dir, subj_ids, spm_exe_script, mcr_path, name='preprocess_workflow'):

    ppwf = pe.Workflow(name=name)
    ppwf.base_dir=work_dir


    inputnode = pe.Node(interface=IdentityInterface(fields=['subj_ids', 'outputdir']), name='inputnode')
    inputnode.iterables = [('subj_ids', subj_ids)]

    templates={"T1":"{subject_id}/T1*.nii.gz"}

    fileselector = pe.Node(SelectFiles(templates), name='fileselect')
    fileselector.inputs.base_directory = scans_dir


    #SPM VBM process: this is removed from here to be run as a separate pipeline.
    vbmpreprocess = pe.Node(interface=VBMPreprocess(), name='vbmpreprocess')
    vbmpreprocess.inputs.spm_exe_script = spm_exe_script
    vbmpreprocess.inputs.mcr_path = mcr_path

    calc_etiv = pe.Node(interface=VBMETIV(), name='calc_etiv')
    calc_etiv.inputs.spm_exe_script = spm_exe_script
    calc_etiv.inputs.mcr_path = mcr_path
 
    #%# datasinker
    datasink = pe.Node(interface=DataSink(), name='datasinker')
    datasink.inputs.parameterization=False
    datasink.inputs.regexp_substitutions = [(r'/[0-9]+_(.*).nii.gz', r'/\1.nii.gz'),
                                            (r'/[0-9]+_(.*).nii', r'/\1.nii'),
                                            (r'/[0-9]+_(.*).mat', r'/\1.mat'),
                                            (r'/[0-9]+_(.*).csv', r'/\1.csv')]

    #%# workflow connections
    ppwf.connect(inputnode       , 'subj_ids'     , fileselector  , 'subject_id')
    ppwf.connect(fileselector    , 'T1'           , vbmpreprocess , 't1file')
    ppwf.connect(vbmpreprocess   , 'segmat'       , calc_etiv     , 'in_segmat')

    ppwf.connect(inputnode       , 'subj_ids'       , datasink, 'container')
    ppwf.connect(inputnode       , 'outputdir'      , datasink, 'base_directory')

    ppwf.connect(vbmpreprocess   , 'wcfiles'      , datasink      , '@wcfiles')
    ppwf.connect(vbmpreprocess   , 'mwcfiles'     , datasink      , '@mwcfiles')
    ppwf.connect(vbmpreprocess   , 'segmat'       , datasink      , '@segmat')
    ppwf.connect(vbmpreprocess   , 'wt1file'      , datasink      , '@wt1file')
    ppwf.connect(vbmpreprocess   , 'qcfile'       , datasink      , '@qcfile')
    ppwf.connect(calc_etiv       , 'etiv'         , datasink      , '@etiv')

    
    return ppwf

