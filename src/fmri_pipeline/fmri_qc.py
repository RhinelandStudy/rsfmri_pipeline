#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 10 11:22:26 2021

@author: hweeling


Copyright 2023 Population Health Sciences, German Center for Neurodegenerative Diseases (DZNE)
Licensed under the Apache License, Version 2.0 (the "License"); 
you may not use this file except in compliance with the License. 
You may obtain a copy of the License at  http://www.apache.org/licenses/LICENSE-2.0 
Unless required by applicable law or agreed to in writing, software distributed under the 
License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and limitations under the License.

"""

from nipype.pipeline import engine as pe
from nipype.interfaces import io as nio
from nipype.interfaces import utility as niu
from nipype.algorithms.confounds import ComputeDVARS
import os.path as op
import numpy as np
import nibabel as nib

from .configoptions import BOLD_MASK
from .spm_interface import CalculateQC

def fmri_qc_workflow(spm_exe_script,mcr_path,name='func_qc_wf'):

    workflow = pe.Workflow(name=name)

    inputnode = pe.Node(niu.IdentityInterface(fields=['swufileb4ica', 'swufileafterica','graymatterfile','whitematterfile','csffile','rpfile','segmentfile']), name='inputnode')
    
    outputnode = pe.Node(niu.IdentityInterface(
        fields=['fdcarpb4ica','fdcarpafterica','allmetrics']), name='outputnode')

    # calculate Compute DVARS before ICA
    #dvars = pe.Node(ComputeDVARS(), name='dvarsb4ica')
    #dvars.inputs.in_mask=BOLD_MASK
    #dvars.inputs.save_plot = False
    #dvars.inputs.save_all = True

    # calculate Compute DVARS after ICA
    #dvarsa = pe.Node(ComputeDVARS(),name='dvarsafterica')
    #dvarsa.inputs.in_mask=BOLD_MASK
    #dvarsa.inputs.save_plot = False
    #dvarsa.inputs.save_all = True


    # Compute GSR
    #gsr = pe.Node(niu.Function(input_names=['swufileb4ica','swufileafterica','maskfile'],
    #                           output_names=['gsryoutput'],
    #                           function=compute_gsr),name='gsr')
    #gsr.inputs.maskfile=BOLD_MASK

    calc_qc = pe.Node(CalculateQC(), name='calc_qc')
    calc_qc.inputs.spm_exe_script=spm_exe_script
    calc_qc.inputs.mcr_path=mcr_path

    workflow.connect([
        #(inputnode,   dvars,  [('swufileb4ica', 'in_file')]),
        #(inputnode,   dvarsa, [('swufileafterica', 'in_file')]),
        #(inputnode,   gsr,    [('swufileb4ica','swufileb4ica'),
        #                       ('swufileafterica','swufileafterica')]),

        (inputnode,  calc_qc, [('swufileb4ica', 'swufileb4ica')]),
        (inputnode,  calc_qc, [('swufileafterica', 'swufileafterica')]),
        (inputnode,  calc_qc, [('graymatterfile', 'graymatterfile')]),
        (inputnode,  calc_qc, [('whitematterfile', 'whitematterfile')]),
        (inputnode,  calc_qc, [('csffile', 'csffile')]),
        (inputnode,  calc_qc, [('rpfile', 'realignmentfile')]),
        (inputnode,  calc_qc, [('segmentfile', 'segmentfile')]),

        #(dvars,      calc_qc, [('out_all', 'dvarsfileb4ica')]),
        #(dvarsa,     calc_qc, [('out_all', 'dvarsfileafterica')]),
        #(gsr,        calc_qc, [('gsryoutput','gsryoutputfile')]),

        (calc_qc,  outputnode, [('fdcarpb4ica', 'fdcarpb4ica')]),
        (calc_qc,  outputnode, [('fdcarpafterica', 'fdcarpafterica')]),
        (calc_qc,  outputnode, [('allmetrics','allmetrics')])
    ])

    return workflow
             

def compute_gsr(swufileb4ica,swufileafterica,maskfile):

    import os.path as op
    import numpy as np
    import nibabel as nib

    epi_data1 = nib.loadsave.load(swufileb4ica)
    epi_data1 = epi_data1.get_data()
    epi_data2 = nib.loadsave.load(swufileb4ica)
    epi_data2 = epi_data2.get_data()

    mask = nib.loadsave.load(maskfile)
    mask = mask.get_data()

    axis = 1 #y-direction

    # Roll data of mask through the appropriate axis
    n2_mask = np.roll(mask, mask.shape[axis] // 2, axis=axis)

    # Step 3: remove from n2_mask pixels inside the brain
    n2_mask = n2_mask * (1 - mask)

    # Step 4: non-ghost background region is labeled as 2
    n2_mask = n2_mask + 2 * (1 - n2_mask - mask)

    # Step 5: signal is the entire foreground image
    ghost1 = np.mean(epi_data1[n2_mask == 1]) - np.mean(epi_data1[n2_mask == 2])
    signal1 = np.median(epi_data1[n2_mask == 0])
    
    gsr_y1 = ghost1 / signal1

    # Step 5: signal is the entire foreground image
    ghost2 = np.mean(epi_data2[n2_mask == 1]) - np.mean(epi_data2[n2_mask == 2])
    signal2 = np.median(epi_data2[n2_mask == 0])
    
    gsr_y2 = ghost2 / signal2

    names = np.array(['gsr-y-b4ica', 'gsr-y-afterica'])
    floats = np.array([gsr_y1, gsr_y2])
    ab = np.zeros(names.size, dtype= [('var1', 'U6'), ('var2', float)])
    ab['var1'] = names
    ab['var2'] = floats
    np.savetxt('gsryoutput.txt', ab, fmt='%20s %10.5f')

    gsryoutput=op.abspath('gsryoutput.txt')                              
    
    return gsryoutput
