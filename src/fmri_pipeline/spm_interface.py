# -*- coding: utf-8 -*-
"""
This is matlab code wrapper interface for nipype 
"""

import os
from glob import glob
from nipype.interfaces.base import BaseInterface, BaseInterfaceInputSpec, File, Directory, traits, TraitedSpec, InputMultiPath
from nipype.interfaces.base import TraitedSpec,CommandLineInputSpec, CommandLine

class CalculateQCInputSpec(CommandLineInputSpec):

    swufileb4ica    = File(exists = True , desc='input swuRestingState before ICA',mandatory=True )
    swufileafterica = File(exists = True , desc='input swuRestingState after ICA',mandatory=True )

    graymatterfile  = File(exists = True , desc='input wc1T1',mandatory=True )
    whitematterfile = File(exists = True , desc='input wc2T1',mandatory=True )
    csffile         = File(exists = True , desc='input wc3T1',mandatory=True )

    realignmentfile   = File(exists = True , desc='input rp_RestingState motion file',mandatory=True )
    segmentfile       = File(exists = True , desc='seg8.mat file',mandatory=True )

    spm_exe_script = File(exists=True, desc='FileName of SPM12 exe script (run_spm12.sh)',mandatory=True)
    mcr_path = Directory(exists=True, desc='MCR path (/path/to/MCR2019b/v97', mandatory=True)


class CalculateQCOutputSpec(TraitedSpec):


    fdcarpb4ica    =  File(exists = True, desc='output FD_carpet_b4ica.pdf')
    fdcarpafterica =  File(exists = True, desc='output FD_carpet_afterica.pdf')
    allmetrics     =  File(exists = True, desc='output all_metrics.csv')


class CalculateQC(CommandLine):

    input_spec = CalculateQCInputSpec
    output_spec = CalculateQCOutputSpec
    _cmd = ' script spmscript.m'

    def _my_script(self):
        swb4=self.inputs.swufileb4ica
        swaf=self.inputs.swufileafterica
        gm=self.inputs.graymatterfile
        wm=self.inputs.whitematterfile
        csf=self.inputs.csffile
        rp=self.inputs.realignmentfile
        seg=self.inputs.segmentfile

        script = """ 
        warning('off','all');
        disp('Starting Calculate QC preprocess...');
        calc_qualitycontrol('%s','%s','%s','%s','%s','%s','%s');       
        """%(swb4,swaf,gm,wm,csf,rp,seg)
        return script

    def _run_interface(self, runtime):
        myscript = self._my_script()
        with open('spmscript.m', 'w') as outf:
            outf.write(myscript)
        self._cmd = self.inputs.spm_exe_script + ' ' + self.inputs.mcr_path + self._cmd
        runtime = super(CalculateQC, self)._run_interface(runtime)
        if runtime.stderr:
           self.raise_exception(runtime)
        return runtime

    def _list_outputs(self):
        outputs = self._outputs().get()
        outputs['fdcarpb4ica'] = os.path.abspath('FD_carpet_b4ica.pdf')
        outputs['fdcarpafterica'] = os.path.abspath('FD_carpet_afterica.pdf')
        outputs['allmetrics'] = os.path.abspath('all_metrics.csv')
        return outputs




class VBMETIVInputSpec(CommandLineInputSpec):

    in_segmat      = File(exists = True , desc='input T1 seg8.mat',mandatory=True )
    spm_exe_script = File(exists=True, desc='FileName of SPM12 exe script (run_spm12.sh)',mandatory=True)
    mcr_path = Directory(exists=True, desc='MCR path (/path/to/MCR2019b/v97', mandatory=True)

class VBMETIVOutputSpec(TraitedSpec):

    etiv      = File(exists = True, desc='output etiv.csv file')

class VBMETIV(CommandLine):
    input_spec = VBMETIVInputSpec
    output_spec = VBMETIVOutputSpec
    _cmd = ' script spmscript.m'

    def _my_script(self):
        in_segmat=self.inputs.in_segmat

        script = """ 
        warning('off','all');
        disp('Starting VBM eTIV preprocess...');
        SPM_eTIV('%s');       
        """%(in_segmat)
        return script

    def _run_interface(self, runtime):
        myscript = self._my_script()
        with open('spmscript.m', 'w') as outf:
            outf.write(myscript)
        self._cmd = self.inputs.spm_exe_script + ' ' + self.inputs.mcr_path + self._cmd
        runtime = super(VBMETIV, self)._run_interface(runtime)
        if runtime.stderr:
           self.raise_exception(runtime)
        return runtime

    def _list_outputs(self):
        outputs = self._outputs().get()
        outputs['etiv'] = os.path.abspath(os.path.join(os.path.dirname(self.inputs.in_segmat),'etiv.csv'))
        return outputs



class VBMPreprocessInputSpec(CommandLineInputSpec):

    t1file         = File(exists = True , desc='input T1 image' )
    spm_exe_script = File(exists=True, desc='FileName of SPM12 exe script (run_spm12.sh)',mandatory=True)
    mcr_path = Directory(exists=True, desc='MCR path (/path/to/MCR2019b/v97', mandatory=True)

class VBMPreprocessOutputSpec(TraitedSpec):

    wcfiles     = InputMultiPath(File(exists=True, desc='output smoothed T1 files.'))
    mwcfiles    = InputMultiPath(File(exists=True, desc='output smoothed T1 files.'))
    segmat      = File(exists = True, desc='output .mat file')
    qcfile      = File(exists = True, desc='output norm qc pdf.')
    wt1file     = File(exists = True, desc='output norm T1 file.')

class VBMPreprocess(CommandLine):
    input_spec = VBMPreprocessInputSpec
    output_spec = VBMPreprocessOutputSpec
    _cmd = ' script spmscript.m'

    def _my_script(self):
        t1file=self.inputs.t1file

        script = """ 
        warning('off','all');
        disp('Starting VBM preprocess...');
        VBM_preprocess_pipeline('%s');       
        """%(t1file)
        return script


    def _run_interface(self, runtime):
        myscript = self._my_script()
        with open('spmscript.m', 'w') as outf:
            outf.write(myscript)
        self._cmd = self.inputs.spm_exe_script + ' ' + self.inputs.mcr_path + self._cmd
        runtime = super(VBMPreprocess, self)._run_interface(runtime)
        if runtime.stderr:
           self.raise_exception(runtime)
        return runtime

    def _list_outputs(self):
        outputs = self._outputs().get()
        outputs['wcfiles']       = [os.path.abspath(f) for f in glob(os.path.join(os.getcwd(),'wc*T1*.nii'))]
        outputs['mwcfiles']      = [os.path.abspath(f) for f in glob(os.path.join(os.getcwd(),'mwc*T1*.nii'))]
        outputs['segmat']        = [os.path.abspath(f) for f in glob(os.path.join(os.getcwd(),'T1*_seg8.mat'))][0]
        outputs['wt1file']       = [os.path.abspath(f) for f in glob(os.path.join(os.getcwd(),'wT1*.nii'))][0]
        outputs['qcfile']        = os.path.abspath('normalized_VBM_QC.pdf')

        return outputs

class SPMPreprocessInputSpec(CommandLineInputSpec):

    restfile      = File(exists= True,  desc = 'Name of input restingstate nii time series.')
    strucfile     = File(exists= True,  desc = 'Name of input structural file.')
    magfile       = File(exists= True,  desc = 'Name of input magnitude file.') 
    phasefile     = File(exists= True,  desc = 'Name of input phase image file.')
    fieldmapfile  = File(exists= True,  desc = 'Name of input fieldmap file.')
    TR            = traits.String(desc='TR')
    TE            = traits.String(desc='TE')
    ndummy        = traits.String(desc='Number of dummy vols to drop')

    spm_exe_script = File(exists=True, desc='FileName of SPM12 exe script (run_spm12.sh)',mandatory=True)
    mcr_path = Directory(exists=True, desc='MCR path (/path/to/MCR2019b/v97', mandatory=True)

class SPMPreprocessOutputSpec(TraitedSpec):
    motionfile    = File(exists = True, desc='output motion parameters file.')
    fslmotionfile = File(exists = True, desc='output FSLmotion.txt file.')
    swurestfile   = File(exists = True, desc='output swuRestingState.nii file.')
    wt1hiresfile  = File(exists = True, desc='output warped T1 hires file')
    wc1t1file     = File(exists = True, desc='output wc1T1.nii file.')
    wc2t1file     = File(exists = True, desc='output wc2T1.nii file.')
    wc3t1file     = File(exists = True, desc='output wc3T1.nii file.')
    t1segmat      = File(exists = True, desc='output T1_seg8.mat file.')
    strucdir      = Directory(exists = True , desc='T1 directory created' )
    qcfile        = File(exists = True, desc='output norm qc pdf.')


class SPMPreprocess(CommandLine):
    input_spec = SPMPreprocessInputSpec
    output_spec = SPMPreprocessOutputSpec
    _cmd = ' script spmscript.m'

    def _my_script(self):
        restfile=os.path.abspath(self.inputs.restfile)
        strucfile=os.path.abspath(self.inputs.strucfile)
        magfile=os.path.abspath(self.inputs.magfile)
        phasefile=os.path.abspath(self.inputs.phasefile)
        fieldmapfile=os.path.abspath(self.inputs.fieldmapfile)
        TR=self.inputs.TR
        TE=self.inputs.TE
        ndummy=self.inputs.ndummy

        script = """ 
        warning('off','all');
        disp('Starting SPM preprocess...');
        SPM_preprocess_pipeline('%s','%s','%s','%s','%s','%s','%s','%s');	
        """%(restfile,strucfile,magfile,phasefile,fieldmapfile,TR,TE,ndummy)
        return script
    

    def _run_interface(self, runtime):
        #spm12 used is in-house-compiled version with embedded SPM_preprocess and SPM_rsfmri_glm scripts 
        #so we have to force nipype to use mcr with this version of spm
        myscript = self._my_script()
        with open('spmscript.m', 'w') as outf:
            outf.write(myscript)
        self._cmd = self.inputs.spm_exe_script + ' ' + self.inputs.mcr_path + self._cmd
        runtime = super(SPMPreprocess, self)._run_interface(runtime)
        if runtime.stderr:
           self.raise_exception(runtime)
        return runtime

    def _list_outputs(self):
        outputs = self._outputs().get()
        outputs['motionfile']      = [os.path.abspath(f) for f in glob('RestingState/rp*.txt')][0]
        outputs['fslmotionfile']   = os.path.abspath('FSL_motion.txt')
        outputs['swurestfile']     = [os.path.abspath(f) for f in glob('swuRestingState*.nii')][0] 
        outputs['wt1hiresfile']    = [os.path.abspath(f) for f in glob('T1/wT1*.nii')][0]
        outputs['wc1t1file']       = [os.path.abspath(f) for f in glob('T1/wc1T1*.nii')][0]
        outputs['wc2t1file']       = [os.path.abspath(f) for f in glob('T1/wc2T1*.nii')][0]
        outputs['wc3t1file']       = [os.path.abspath(f) for f in glob('T1/wc3T1*.nii')][0]
        outputs['t1segmat']        = [os.path.abspath(f) for f in glob('T1/T1*_seg8.mat')][0]
        outputs['strucdir']        = os.path.abspath(os.path.join(os.getcwd(), 'T1'))
        outputs['qcfile']          = os.path.abspath('normalized_RS_QC.pdf')



        return outputs


class SPMGLMProcessInputSpec(CommandLineInputSpec):

    swurestfile    = File(exists= True,  desc = 'Name of input swuRestingState nii time series.')
    wc2t1file      = File(exists= True,  desc = 'Name of input structural preprocess file wc2T1*.')
    wc3t1file      = File(exists= True,  desc = 'Name of input structural preprocess file wc3T1*.')
    fslmotionfile  = File(exists= True,  desc = 'Name of input motion params file.')
    spm_exe_script = File(exists=True, desc='FileName of SPM12 exe script (run_spm12.sh)',mandatory=True)
    mcr_path = Directory(exists=True, desc='MCR path (/path/to/MCR2019b/v97', mandatory=True)

class SPMGLMProcessOutputSpec(TraitedSpec):
    fc_17raw  = File(exists = True, desc='output functional connectivity file.')
    fc_17full = File(exists = True, desc='output functional connectivity file.')
    fc_7raw   = File(exists = True, desc='output functional connectivity file.')
    fc_7full = File(exists = True, desc='output functional connectivity file.')

class SPMGLMProcess(CommandLine):
    input_spec = SPMGLMProcessInputSpec
    output_spec = SPMGLMProcessOutputSpec
    _cmd = ' script spmscript.m'

    def _my_script(self):
        swurestfile  = os.path.abspath(self.inputs.swurestfile)
        wc2t1file    = os.path.abspath(self.inputs.wc2t1file)
        wc3t1file    = os.path.abspath(self.inputs.wc3t1file)
        fslmotionfile= os.path.abspath(self.inputs.fslmotionfile)

        script = """ 
        warning('off','all');
        disp('Starting SPM GLM process...');
        SPM_rsfmri_glm_200_17_7networks('%s','%s','%s','%s');;       
        """%(swurestfile,fslmotionfile,wc2t1file,wc3t1file)
        return script

    def _run_interface(self, runtime):
        #spm12 used is in-house-compiled version with embedded SPM_preprocess and SPM_rsfmri_glm scripts 
        #so we have to force nipype to use mcr with this version of spm
        myscript = self._my_script()
        with open('spmscript.m', 'w') as outf:
            outf.write(myscript)
        self._cmd = self.inputs.spm_exe_script + ' ' + self.inputs.mcr_path + self._cmd
        runtime = super(SPMGLMProcess, self)._run_interface(runtime)
        if runtime.stderr:
           self.raise_exception(runtime)
        return runtime

    def _list_outputs(self):
        outputs = self._outputs().get()
        outputs['fc_17raw']   = os.path.abspath('FFX2/FC_Schaefer_17Networks_200p_raw.csv')
        outputs['fc_17full']   = os.path.abspath('FFX2/FC_Schaefer_17Networks_200p_fullcorr.csv')
        outputs['fc_7raw']   = os.path.abspath('FFX2/FC_Schaefer_7Networks_200p_raw.csv')
        outputs['fc_7full']   = os.path.abspath('FFX2/FC_Schaefer_7Networks_200p_fullcorr.csv')

        return outputs

class SPMFCMeansInputSpec(CommandLineInputSpec):
    fc_17full    = File(exists= True,  desc = 'Name of input 17N FC file.')
    fc_7full     = File(exists= True,  desc = 'Name of input 7N FC file.')
    spm_exe_script = File(exists=True, desc='FileName of SPM12 exe script (run_spm12.sh)',mandatory=True)
    mcr_path = Directory(exists=True, desc='MCR path (/path/to/MCR2019b/v97', mandatory=True)


class SPMFCMeansOutputSpec(TraitedSpec):
    fc17meanfc  = File(exists = True, desc='output functional connectivity means.')
    fc7meanfc   = File(exists = True, desc='output functional connectivity means.')

class SPMFCMeans(CommandLine):
    input_spec = SPMFCMeansInputSpec
    output_spec = SPMFCMeansOutputSpec
    _cmd = ' script spmscript.m'

    def _my_script(self):
        fc_17full  = os.path.abspath(self.inputs.fc_17full)
        fc_7full   = os.path.abspath(self.inputs.fc_7full)

        script = """ 
        warning('off','all');
        disp('Starting FC means process...');
        extract_Schaefer_17networks('%s');
        extract_Schaefer_7networks('%s');       
        """%(fc_17full, fc_7full)
        return script

    def _run_interface(self, runtime):
        #spm12 used is in-house-compiled version with embedded SPM_preprocess and SPM_rsfmri_glm scripts 
        #so we have to force nipype to use mcr with this version of spm
        myscript = self._my_script()
        with open('spmscript.m', 'w') as outf:
            outf.write(myscript)
        self._cmd = self.inputs.spm_exe_script + ' ' + self.inputs.mcr_path + self._cmd
        runtime = super(SPMFCMeans, self)._run_interface(runtime)
        if runtime.stderr:
           self.raise_exception(runtime)
        return runtime

    def _list_outputs(self):
        outputs = self._outputs().get()
        outputs['fc17meanfc']  = os.path.abspath('FC_Schaefer_17Networks_200p_fullcorr_means.json')
        outputs['fc7meanfc']   = os.path.abspath('FC_Schaefer_7Networks_200p_fullcorr_means.json')

        return outputs




class SPMGLMProcessSTInputSpec(CommandLineInputSpec):

    swurestfile    = File(exists= True,  desc = 'Name of input swuRestingState nii time series.')
    wc2t1file      = File(exists= True,  desc = 'Name of input structural preprocess file wc2T1*.')
    wc3t1file      = File(exists= True,  desc = 'Name of input structural preprocess file wc3T1*.')
    fslmotionfile  = File(exists= True,  desc = 'Name of input motion params file.')
    spm_exe_script = File(exists=True, desc='FileName of SPM12 exe script (run_spm12.sh)',mandatory=True)
    mcr_path = Directory(exists=True, desc='MCR path (/path/to/MCR2019b/v97', mandatory=True)

class SPMGLMProcessSTOutputSpec(TraitedSpec):
    fc_conn_file_200p = File(exists = True, desc='output functional connectivity file 200p.')
    fc_conn_file_100p = File(exists = True, desc='output functional connectivity file 100p.')

class SPMGLMProcessST(CommandLine):
    input_spec = SPMGLMProcessSTInputSpec
    output_spec = SPMGLMProcessSTOutputSpec
    _cmd = ' script spmscript.m'

    def _my_script(self):
        swurestfile  = os.path.abspath(self.inputs.swurestfile)
        wc2t1file    = os.path.abspath(self.inputs.wc2t1file)
        wc3t1file    = os.path.abspath(self.inputs.wc3t1file)
        fslmotionfile= os.path.abspath(self.inputs.fslmotionfile)

        script = """ 
        warning('off','all');
        disp('Starting SPM GLM_200_100 process...');
        SPM_rsfmri_glm_200_100('%s','%s','%s','%s');;       
        """%(swurestfile,fslmotionfile,wc2t1file,wc3t1file)
        return script

    def _run_interface(self, runtime):
        #spm12 used is in-house-compiled version with embedded SPM_preprocess and SPM_rsfmri_glm scripts 
        #so we have to force nipype to use mcr with this version of spm
        myscript = self._my_script()
        with open('spmscript.m', 'w') as outf:
            outf.write(myscript)
        self._cmd = self.inputs.spm_exe_script + ' ' + self.inputs.mcr_path + self._cmd
        runtime = super(SPMGLMProcessST, self)._run_interface(runtime)
        if runtime.stderr:
           self.raise_exception(runtime)
        return runtime

    def _list_outputs(self):
        outputs = self._outputs().get()
        outputs['fc_conn_file_200p']   = os.path.abspath('FFX2/FC_Schaefer_17Networks_200p_simple.csv')
        outputs['fc_conn_file_100p']   = os.path.abspath('FFX2/FC_Schaefer_17Networks_100p_simple.csv')


        return outputs


class SPMGLMProcessVSInputSpec(CommandLineInputSpec):

    swurestfile    = File(exists= True,  desc = 'Name of input swuRestingState nii time series.')
    wc2t1file      = File(exists= True,  desc = 'Name of input structural preprocess file wc2T1*.')
    wc3t1file      = File(exists= True,  desc = 'Name of input structural preprocess file wc3T1*.')
    fslmotionfile  = File(exists= True,  desc = 'Name of input motion params file.')
    spm_exe_script = File(exists=True, desc='FileName of SPM12 exe script (run_spm12.sh)',mandatory=True)
    mcr_path = Directory(exists=True, desc='MCR path (/path/to/MCR2019b/v97', mandatory=True)

class SPMGLMProcessVSOutputSpec(TraitedSpec):
    vs_conn_file = File(exists = True, desc='output functional connectivity file.')
    rvs_conn_file = File(exists = True, desc='output functional connectivity file.')
    lvs_conn_file = File(exists = True, desc='output functional connectivity file.')

class SPMGLMProcessVS(CommandLine):
    input_spec  = SPMGLMProcessVSInputSpec
    output_spec = SPMGLMProcessVSOutputSpec
    _cmd = ' script spmscript.m'

    def _my_script(self):
        swurestfile  = os.path.abspath(self.inputs.swurestfile)
        wc2t1file    = os.path.abspath(self.inputs.wc2t1file)
        wc3t1file    = os.path.abspath(self.inputs.wc3t1file)
        fslmotionfile= os.path.abspath(self.inputs.fslmotionfile)

        script = """ 
        warning('off','all');
        disp('Starting SPM GLM VS process...');
        SPM_rsfmri_glm_VStriatum('%s','%s','%s','%s');;       
        """%(swurestfile,fslmotionfile,wc2t1file,wc3t1file)
        return script

    def _run_interface(self, runtime):
        #spm12 used is in-house-compiled version with embedded SPM_preprocess and SPM_rsfmri_glm_VS scripts 
        #so we have to force nipype to use mcr with this version of spm
        myscript = self._my_script()
        with open('spmscript.m', 'w') as outf:
            outf.write(myscript)
        self._cmd = self.inputs.spm_exe_script + ' ' + self.inputs.mcr_path + self._cmd
        runtime = super(SPMGLMProcessVS, self)._run_interface(runtime)
        if runtime.stderr:
           self.raise_exception(runtime)
        return runtime

    def _list_outputs(self):
        outputs = self._outputs().get()
        outputs['vs_conn_file']    = os.path.abspath('VStiatum/con_0001.nii')
        outputs['rvs_conn_file']   = os.path.abspath('R_VStriatum/con_0001.nii')
        outputs['lvs_conn_file']   = os.path.abspath('L_VStriatum/con_0001.nii')

        return outputs


class SPMGLMProcessGordonInputSpec(CommandLineInputSpec):

    swurestfile    = File(exists= True,  desc = 'Name of input swuRestingState nii time series.')
    wc2t1file      = File(exists= True,  desc = 'Name of input structural preprocess file wc2T1*.')
    wc3t1file      = File(exists= True,  desc = 'Name of input structural preprocess file wc3T1*.')
    fslmotionfile  = File(exists= True,  desc = 'Name of input motion params file.')
    spm_exe_script = File(exists=True, desc='FileName of SPM12 exe script (run_spm12.sh)',mandatory=True)
    mcr_path = Directory(exists=True, desc='MCR path (/path/to/MCR2019b/v97', mandatory=True)

class SPMGLMProcessGordonOutputSpec(TraitedSpec):
    fc_conn_file = File(exists = True, desc='output functional connectivity file.')

class SPMGLMProcessGordon(CommandLine):
    input_spec = SPMGLMProcessGordonInputSpec
    output_spec = SPMGLMProcessGordonOutputSpec
    _cmd = ' script spmscript.m'

    def _my_script(self):
        swurestfile  = os.path.abspath(self.inputs.swurestfile)
        wc2t1file    = os.path.abspath(self.inputs.wc2t1file)
        wc3t1file    = os.path.abspath(self.inputs.wc3t1file)
        fslmotionfile= os.path.abspath(self.inputs.fslmotionfile)

        script = """ 
        warning('off','all');
        disp('Starting SPM GLM process...');
        SPM_rsfmri_glm_Gordon_333('%s','%s','%s','%s');;       
        """%(swurestfile,fslmotionfile,wc2t1file,wc3t1file)
        return script

    def _run_interface(self, runtime):
        #spm12 used is in-house-compiled version with embedded SPM_preprocess and SPM_rsfmri_glm scripts 
        #so we have to force nipype to use mcr with this version of spm
        myscript = self._my_script()
        with open('spmscript.m', 'w') as outf:
            outf.write(myscript)
        self._cmd = self.inputs.spm_exe_script + ' ' + self.inputs.mcr_path + self._cmd
        runtime = super(SPMGLMProcessGordon, self)._run_interface(runtime)
        if runtime.stderr:
           self.raise_exception(runtime)
        return runtime

    def _list_outputs(self):
        outputs = self._outputs().get()
        outputs['fc_conn_file']   = os.path.abspath('FFX2/FC_Gordon_333p_simple.csv')

        return outputs

