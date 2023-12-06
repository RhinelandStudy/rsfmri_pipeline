function SPM_preprocess_pipeline(restfile,strucfile,magfile,phasefile,fmapfile,TR,TE,ndummy)

TR=str2num(TR);
TE=str2num(TE);
ndummy=str2num(ndummy);
vox_size = [2.4 2.4 2.4];

fs = filesep;
clear tmpfile;

disp('****************************************************');
disp('INPUT SETTINGS....');
disp('****************************************************');
disp(strcat('RestingState:', restfile));
disp(strcat('Structural:', strucfile));
disp(strcat('Magnitude:', magfile));
disp(strcat('Phase:',phasefile));
disp(strcat('Fmap:',fmapfile));
disp(strcat('TR:',num2str(TR)));
disp(strcat('TE:',num2str(TE)));
disp(strcat('NDummyVols:',num2str(ndummy)));
disp('***************************************************');


SubjfMRI = 'RestingState';
SubjT1 =  'T1';

if isfolder(SubjfMRI) == 0
    mkdir(SubjfMRI);
end
if isfolder(SubjT1) == 0
    mkdir(SubjT1);
end

gunzip(restfile,  pwd);
gunzip(strucfile, pwd);
gunzip(phasefile, pwd);

[fp,fn,ext] = fileparts(restfile);
restfile=fn;
clear fp fn ext;

[fp,fn,ext] = fileparts(strucfile);
strucfile=fn;
clear fp fn ext;

[fp,fn,ext] = fileparts(phasefile);
phasefile=fn;
clear fp fn ext;

if isfile(magfile)
    gunzip(magfile, pwd);
    [fp,fn,ext] = fileparts(magfile);
    magfile=fn;
    clear fp fn ext;
end

if isfile(fmapfile)
    gunzip(fmapfile, pwd);
    [fp,fn,ext] = fileparts(fmapfile);
    fmapfile=fn;
    clear fp fn ext;
end


movefile(strucfile, SubjT1);

matlabbatch{1}.spm.util.split.vol = cellstr(strcat(pwd,fs,restfile, ',1'));
matlabbatch{1}.spm.util.split.outdir = cellstr(SubjfMRI);

disp('*** converting 4D to 3D volumes *** ')
tic
spm_jobman('run', matlabbatch);
toc
clear matlabbatch;

disp('Calculate vdm - for geometric distortion correction ');
disp('this is to check if the real B0 magnitude image file is available.');

% this step is required for those who have bad quality T1 image.
if strfind(strucfile,'T1')
    struct_image = spm_select('ExtFPList',SubjT1,'^T1.*\.nii$');
else
    struct_image = spm_select('ExtFPList',SubjT1,'^Scout.*\.nii$');
end

if strfind(magfile, 'Scout') % if the real B0 magnitude image file is not available, fmap_hz.nii will be used to calculate geometric distortion map

    disp('precalculated fieldmap using Scout as magnitude image.');
    fieldmap_image = spm_select('ExtFPList',strcat(pwd, fs),'^fmap_hz*\.nii$');
    epi_image = spm_select('ExtFPList', SubjfMRI, '^RestingState_00001*\.nii$');
    
    matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.data.precalcfieldmap.precalcfieldmap = cellstr(fieldmap_image);
    matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.data.precalcfieldmap.magfieldmap = {};
    matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsfile = cellstr(strcat(spm('Dir'), fs, 'toolbox', fs, 'FieldMap', fs, 'pm_defaults.m'));
    matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.session.epi = cellstr(epi_image);
    matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.matchvdm = 1;
    matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.writeunwarped = 0;
    matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.anat = cellstr(struct_image);
    matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.matchanat = 1;
    
    disp('*** calculating voxel displacement map for epi ...')
    tic
    spm_jobman('run', matlabbatch);
    toc
    clear tmp1 matlabbatch epi_image fieldmap_image mag_image;
    
else
    
    disp('calculate fieldmap using real magnitude image.');
    
    phase_image = spm_select('ExtFPList', strcat(pwd,fs),'^B0_Phase.*\.nii$');
    mag_image = spm_select('ExtFPList', strcat(pwd,fs), '^B0.nii$',1);
    epi_image = spm_select('ExtFPList', SubjfMRI, '^RestingState_00001*\.nii$');
    
    matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.data.presubphasemag.phase = cellstr(phase_image);
    matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.data.presubphasemag.magnitude = cellstr(mag_image);
    matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.et = [4.92 7.38];
    matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.maskbrain = 0;
    matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.blipdir = -1;
    matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.tert = 50.4;
    matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.epifm = 1;
    matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.ajm = 0;
    matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.uflags.method = 'Mark3D';
    matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.uflags.fwhm = 10;
    matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.uflags.pad = 0;
    matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.uflags.ws = 1;
    matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.mflags.template = {strcat(spm('Dir'), fs, 'toolbox', fs, 'FieldMap', fs, 'T1.nii')};
    matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.mflags.fwhm = 5;
    matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.mflags.nerode = 2;
    matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.mflags.ndilate = 4;
    matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.mflags.thresh = 0.5;
    matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.mflags.reg = 0.02;
    matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.session(1).epi = cellstr(epi_image);
    matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.matchvdm = 1;
    matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.sessname = 'run';
    matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.writeunwarped = 0;
    matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.anat = cellstr(struct_image);
    matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.matchanat = 1;
    
    disp(sprintf('*** gradient field mapping ...'))
    tic
    spm_jobman('run', matlabbatch);
    toc
    clear matlabbatch;
    
end

%%%%%%%%%%%%%%%%%%%%%%
% delete unwarped if written
%%%%%%%%%%%%%%%%%%%%%%
restunwarped=strcat(SubjfMRI,fs,'uRestingState_00001.nii');
if exist(restunwarped, 'file')==2
    delete(restunwarped);
end


%%%%%%%%%%%%%%%%%%%%%%
% Realign and Unwarp %
%%%%%%%%%%%%%%%%%%%%%%
% this is motion correction (I used the mean image to correct for motion) and
% geometric distortion correction using the distortion field map
% inputs: RestingState_00018:01070.nii, vdm5_*.nii
% output: uRestingState_00018:01070.nii
% output: rp*.txt !! this should be saved for future analyses

matlabbatch{1}.spm.spatial.realignunwarp.eoptions.quality = 0.9;
matlabbatch{1}.spm.spatial.realignunwarp.eoptions.sep = 4;
matlabbatch{1}.spm.spatial.realignunwarp.eoptions.fwhm = 5;
matlabbatch{1}.spm.spatial.realignunwarp.eoptions.rtm = 1;
matlabbatch{1}.spm.spatial.realignunwarp.eoptions.einterp = 2;
matlabbatch{1}.spm.spatial.realignunwarp.eoptions.ewrap = [0 0 0];
matlabbatch{1}.spm.spatial.realignunwarp.eoptions.weight = {};
matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.basfcn = [12 12];
matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.regorder = 1;
matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.lambda = 100000;
matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.jm = 0;
matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.fot = [4 5];
matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.sot = [];
matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.uwfwhm = 4;
matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.rem = 1;
matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.noi = 5;
matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.expround = 'First';
matlabbatch{1}.spm.spatial.realignunwarp.uwroptions.uwwhich = [2 1];
matlabbatch{1}.spm.spatial.realignunwarp.uwroptions.rinterp = 4;
matlabbatch{1}.spm.spatial.realignunwarp.uwroptions.wrap = [0 0 0];
matlabbatch{1}.spm.spatial.realignunwarp.uwroptions.prefix = 'u';

% load of images
tmp1 = [] ;
tmp1 = spm_select('ExtFPList',SubjfMRI,'^R.*\.nii$');
vdm5_image = strcat(pwd, fs, '^vdm5.*\.nii$');

matlabbatch{1}.spm.spatial.realignunwarp.data(1).scans  = cellstr(tmp1(ndummy+1:end,:));
matlabbatch{1}.spm.spatial.realignunwarp.data(1).pmscan = cellstr(vdm5_image);

disp('*** realigning (coregistering) ... and unwarping (and reslicing) ...')
tic
spm_jobman('run', matlabbatch);
toc
clear tmp1 matlabbatch vdm5_image;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% creating motion parameter file for FSL melodic ICA %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% by default in SPM, motion parameters are organized in translation
% x, y, z, rotation x, y, z but for FSL melodic ICA, motion
% parameters should be organized in rotation x, y, z, translation
% x, y, z !! this should be saved for future analyses

clear d tmp1 tmp2 tmp3;
d = dir(strcat(SubjfMRI, fs, '*.txt'));
tmp1 = dlmread(char(strcat(d(1).folder, fs, d(1).name)));
tmp3 = tmp1(:,[4:6 1:3]); % for FSL Melodic ICA: rotation x, y, z translation x, y, z
wm=which('writematrix');
if ~isempty(wm)
    writematrix(tmp3, strcat(pwd, fs,'FSL_motion.txt'), 'Delimiter', ' ');
else
    dlmwrite(strcat(pwd, fs, 'FSL_motion.txt'),tmp3, 'delimiter', ' ');
end

clear tmp1 tmp2 tmp3;


%%%%%%%%%%%%%%
% Coregister %
%%%%%%%%%%%%%%
% this is to ensure that T1 image is in the same orientation as EPI
% image since T1 image will be segmented to calculate normlization
% to MNI space parameters for EPI normalization
% inputs: T1.nii (or Scout.nii), meanRestingState*.nii
% there is no output but this coregistration will change the header
% orientation parameters of the T1.nii (or Scout.nii)

matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];
matlabbatch{1}.spm.spatial.coreg.estimate.other = {''} ;


ref_image = spm_select('ExtFPList',SubjfMRI,'^mean.*\.nii$');
% this step is required for those who have bad quality T1 image.
if strfind(strucfile,'T1')
    struct_image = spm_select('ExtFPList',SubjT1,'^T1.*\.nii$');
else
    struct_image = spm_select('ExtFPList',SubjT1,'^Scout.*\.nii$');
end

matlabbatch{1}.spm.spatial.coreg.estimate.source = cellstr(struct_image);
matlabbatch{1}.spm.spatial.coreg.estimate.ref = cellstr(ref_image);

disp('*** Coregistering structural volume to mean image ...')
tic
spm_jobman('run', matlabbatch);
toc
clear matlabbatch;

%%%%%%%%%%%%%%%%%%%%%%%%
% Unified segmentation %
%%%%%%%%%%%%%%%%%%%%%%%%
% this is to segment the T1 image into different tissue classes: c1: GM; c2: WM, c3: CSF
% if T1 image is not available Scout image will be used
% inputs; T1.nii (or Scout.nii)
% outputs: c1:c6.nii, rc1:rc6.nii, iy_T1*.nii, T1*_seg8.mat !! this should be saved for future analyses

matlabbatch{1}.spm.spatial.preproc.channel.vols = cellstr(struct_image);
matlabbatch{1}.spm.spatial.preproc.channel.biasreg = 0.0001;
matlabbatch{1}.spm.spatial.preproc.channel.biasfwhm = 60;
matlabbatch{1}.spm.spatial.preproc.channel.write = [0 0];

% load TPM images
ngaus_tmp = [1 1 2 3 4 2];
native_tmp = [1 1; 1 1; 1 1; 1 1; 1 1; 1 1];
warped_tmp = [0 0; 0 0; 0 0; 0 0; 0 0; 0 0];

for k = 1:6
    matlabbatch{1}.spm.spatial.preproc.tissue(k).tpm = cellstr(spm_select('ExtFPList', strcat(spm('Dir'), fs, 'tpm'), '^TPM.*\.nii$', k));
    matlabbatch{1}.spm.spatial.preproc.tissue(k).ngaus = ngaus_tmp(k);
    matlabbatch{1}.spm.spatial.preproc.tissue(k).native = native_tmp(k,:);
    matlabbatch{1}.spm.spatial.preproc.tissue(k).warped = warped_tmp(k,:);
end

matlabbatch{1}.spm.spatial.preproc.warp.mrf = 1;
matlabbatch{1}.spm.spatial.preproc.warp.cleanup = 1;
matlabbatch{1}.spm.spatial.preproc.warp.reg = [0 0.0001 0.5 0.05 0.2];
matlabbatch{1}.spm.spatial.preproc.warp.affreg = 'mni';
matlabbatch{1}.spm.spatial.preproc.warp.fwhm = 0;
matlabbatch{1}.spm.spatial.preproc.warp.samp = 3;
matlabbatch{1}.spm.spatial.preproc.warp.write = [1 1];

disp(sprintf('*** Segmentation of the coregistered structural volume ...'))
tic
spm_jobman('run', matlabbatch);
toc
clear matlabbatch;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Run Shoot template (existing - using template from CAT12) %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this is to calculate non-linear image registration (spatial
% normalization) for better inter-subject registration
% inputs: rc1:rc2.nii
% outputs: j_rc1T1*_Template.nii, y_rc1T1*_Template.nii, v_rc1T1*_Template.nii !! this should be saved for future analyses

clear SubjGM SubjWM;

SubjGM = spm_select('ExtFPList', strcat(pwd, fs, 'T1'),'^rc1.*\.nii$');
SubjWM = spm_select('ExtFPList', strcat(pwd, fs, 'T1'),'^rc2.*\.nii$');

matlabbatch{1}.spm.tools.shoot.warp1.images{1} = cellstr(SubjGM);
matlabbatch{1}.spm.tools.shoot.warp1.images{2} = cellstr(SubjWM);
matlabbatch{1}.spm.tools.shoot.warp1.templates = cellstr(strcat(spm('Dir'), fs, 'toolbox', fs, 'Shoot', fs, 'Template_0_IXI555_MNI152_GS.nii'));

tic
spm_jobman('run', matlabbatch);
toc
clear matlabbatch;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Write normalized and smoothing %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% to normalize EPI images into MNI space
% inputs: y_rc1T1*_Template.nii, uRestingState_00018:01070.nii
% outputs: swuRestingState_00018:01070.nii

clear deformfile rest_images;
deformfile = spm_select('FPList', SubjT1, '^y_rc*.*\.nii$');
rest_images = spm_select('FPList', SubjfMRI, '^uR*.*\.nii$');

matlabbatch{1}.spm.tools.shoot.norm.data.subj.deformation = cellstr(deformfile);
matlabbatch{1}.spm.tools.shoot.norm.data.subj.images = cellstr(rest_images);
matlabbatch{1}.spm.tools.shoot.norm.template = cellstr(strcat(spm('Dir'), fs, 'toolbox', fs, 'Shoot', fs, 'Template_4_IXI555_MNI152_GS.nii'));
matlabbatch{1}.spm.tools.shoot.norm.vox = vox_size;
matlabbatch{1}.spm.tools.shoot.norm.bb = [NaN NaN NaN; NaN NaN NaN];
matlabbatch{1}.spm.tools.shoot.norm.preserve = 0;
matlabbatch{1}.spm.tools.shoot.norm.fwhm = 6;

disp('*** writing normalised functional images ...')
tic
spm_jobman('run', matlabbatch);
toc
clear tmp1 matlabbatch;

%%%%%%%%%%%%%%%%%%%%%%%
% converting 3D to 4D %
%%%%%%%%%%%%%%%%%%%%%%%
% converting images back into 4D for FSL Melodic ICA
% inputs: swuRestingState_00018:01070.nii
% output: swuRestingState.nii !! this should be saved for future analyses

matlabbatch{1}.spm.util.cat.vols = cellstr(spm_select('FPList', SubjfMRI, '^swuR.*\.nii$'));
matlabbatch{1}.spm.util.cat.name = strcat(pwd, fs, 'swuRestingState.nii');
matlabbatch{1}.spm.util.cat.dtype = 0;
matlabbatch{1}.spm.util.cat.RT = TR;

tic
spm_jobman('run', matlabbatch);
toc
clear matlabbatch;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Write normalized GM, WM, CSF + T1 images but without smoothing %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% To normalize structural images but without smoothing to use as a
% mask image in some cases. E.g. ROI signal extraction of WM and
% CSF included as regressors in GLM at a later stage.
% inputs: y_rc1T1*_Template.nii, c1:c3.nii, T1.nii (if available)
% output: wc1:c3.nii, wT1.nii (if available) !! this should be saved for future analyses

clear deformfile rest_images SubjGM SubjWM SubjCSF SubjT1file;
deformfile = spm_select('FPList', SubjT1, '^y_rc*.*\.nii$');
SubjGM = spm_select('FPList', SubjT1, '^c1T1*.*\.nii$');
SubjWM = spm_select('FPList', SubjT1, '^c2T1*.*\.nii$');
SubjCSF = spm_select('FPList', SubjT1, '^c3T1*.*\.nii$');
SubjT1file = spm_select('FPList', SubjT1, '^T1*.*\.nii$');

if isempty(SubjT1file)
    allfiles = {SubjGM; SubjWM; SubjCSF};
else
    allfiles = {SubjGM; SubjWM; SubjCSF; SubjT1file};
end

matlabbatch{1}.spm.tools.shoot.norm.data.subj.deformation = cellstr(deformfile);
matlabbatch{1}.spm.tools.shoot.norm.data.subj.images = allfiles;
matlabbatch{1}.spm.tools.shoot.norm.template = cellstr(strcat(spm('Dir'), fs, 'toolbox', fs, 'Shoot', fs, 'Template_4_IXI555_MNI152_GS.nii'));
matlabbatch{1}.spm.tools.shoot.norm.vox = vox_size;
matlabbatch{1}.spm.tools.shoot.norm.bb = [NaN NaN NaN; NaN NaN NaN];
matlabbatch{1}.spm.tools.shoot.norm.preserve = 0;
matlabbatch{1}.spm.tools.shoot.norm.fwhm = 0;

disp('*** writing normalised GM + WM + CSF + T1 images ...')
tic
spm_jobman('run', matlabbatch);
toc
clear tmp1 matlabbatch;

%%%%%%%%%%%%%%%%%%
% computing eTIV %
%%%%%%%%%%%%%%%%%%

seg8matfile = spm_select('FPList', SubjT1, 'T1.*_seg8.mat');

matlabbatch{1}.spm.util.tvol.matfiles = cellstr(seg8matfile);
matlabbatch{1}.spm.util.tvol.tmax = 3;
matlabbatch{1}.spm.util.tvol.mask = cellstr(strcat(spm('Dir'), fs, 'tpm', fs, 'mask_ICV.nii,1'));
matlabbatch{1}.spm.util.tvol.outf = strcat(pwd,fs,'etiv');

tic
spm_jobman('run', matlabbatch);
toc
clear matlabbatch;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Quality control check for normalization %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

matlabbatch{1}.spm.util.checkreg.data = {spm_select('FPListRec', SubjT1, '^wT1*.*\.nii$');...
    spm_select('FPListRec', SubjT1, '^wc1.*\.nii$');...
    spm_select('FPListRec', SubjT1, '^wc2.*\.nii$');...
    spm_select('FPListRec', SubjT1, '^wc3.*\.nii$');...
    spm_select('FPListRec', SubjfMRI, '^swuRestingState_00018.*\.nii$');...
    };
spm_jobman('run', matlabbatch);
clear matlabbatch;

printfname = strcat(pwd, fs, 'normalized_RS_QC.pdf');
print(printfname, '-dpdf');
quit;

end

