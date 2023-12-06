function VBM_preprocess_pipeline(strucfile)

fs = filesep;
clear tmpfile;

disp('****************************************************');
disp('INPUT SETTINGS....');
disp('****************************************************');
disp(strcat('Structural:', strucfile));
disp('***************************************************');

gunzip(strucfile, pwd);

[fp,fn,ext] = fileparts(strucfile);
strucfile=fn;
clear fp fn ext;

%%%%%%%%%%%%%%%%%%%%%%%%
% Unified segmentation %
%%%%%%%%%%%%%%%%%%%%%%%%
% this is to segment the T1 image into different tissue classes: c1: GM; c2: WM, c3: CSF
% if T1 image is not available Scout image will be used
% inputs; T1.nii (or Scout.nii)
% outputs: c1:c6.nii, rc1:rc6.nii, iy_T1*.nii, T1*_seg8.mat !! this should be saved for future analyses

    struct_image = spm_select('ExtFPList', pwd,'^T1.*\.nii$');

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

SubjGM = spm_select('ExtFPList', pwd,'^rc1.*\.nii$');
SubjWM = spm_select('ExtFPList', pwd,'^rc2.*\.nii$');

matlabbatch{1}.spm.tools.shoot.warp1.images{1} = cellstr(SubjGM);
matlabbatch{1}.spm.tools.shoot.warp1.images{2} = cellstr(SubjWM);
matlabbatch{1}.spm.tools.shoot.warp1.templates = cellstr(strcat(spm('Dir'), fs, 'toolbox', fs, 'Shoot', fs, 'Template_0_IXI555_MNI152_GS.nii'));

tic
spm_jobman('run', matlabbatch);
toc
clear matlabbatch;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Write modulated normalized GM, WM, CSF images and smoothing %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% To normalize structural images 
% inputs: y_rc1T1*_Template.nii, c1:c3.nii, T1.nii
% output: smwc1:c3.nii !! this should be saved for future analyses

clear deformfile SubjGM SubjWM SubjCSF SubjT1file;
deformfile = spm_select('FPList', pwd, '^y_rc*.*\.nii$');
SubjGM = spm_select('FPList', pwd, '^c1T1*.*\.nii$');
SubjWM = spm_select('FPList', pwd, '^c2T1*.*\.nii$');
SubjCSF = spm_select('FPList', pwd, '^c3T1*.*\.nii$');
SubjT1file = spm_select('FPList', pwd, '^T1*.*\.nii$');

matlabbatch{1}.spm.tools.shoot.norm.data.subj.deformation = cellstr(deformfile);
matlabbatch{1}.spm.tools.shoot.norm.data.subj.images = {SubjGM; SubjWM; SubjCSF};
matlabbatch{1}.spm.tools.shoot.norm.template = cellstr(strcat(spm('Dir'), fs, 'toolbox', fs, 'Shoot', fs, 'Template_4_IXI555_MNI152_GS.nii'));
matlabbatch{1}.spm.tools.shoot.norm.vox = [1 1 1];
matlabbatch{1}.spm.tools.shoot.norm.bb = [NaN NaN NaN; NaN NaN NaN];
matlabbatch{1}.spm.tools.shoot.norm.preserve = 1;
matlabbatch{1}.spm.tools.shoot.norm.fwhm = 0;

disp('*** writing normalised GM + WM + CSF images ...')
tic
spm_jobman('run', matlabbatch);
toc
clear tmp1 matlabbatch;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Write non-modulated normalized GM, WM, CSF images and smoothing %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

matlabbatch{1}.spm.tools.shoot.norm.data.subj.deformation = cellstr(deformfile);
matlabbatch{1}.spm.tools.shoot.norm.data.subj.images = {SubjGM; SubjWM; SubjCSF};
matlabbatch{1}.spm.tools.shoot.norm.template = cellstr(strcat(spm('Dir'), fs, 'toolbox', fs, 'Shoot', fs, 'Template_4_IXI555_MNI152_GS.nii'));
matlabbatch{1}.spm.tools.shoot.norm.vox = [1 1 1];
matlabbatch{1}.spm.tools.shoot.norm.bb = [NaN NaN NaN; NaN NaN NaN];
matlabbatch{1}.spm.tools.shoot.norm.preserve = 0;
matlabbatch{1}.spm.tools.shoot.norm.fwhm = 0;

disp('*** writing normalised GM + WM + CSF images ...')
tic
spm_jobman('run', matlabbatch);
toc
clear tmp1 matlabbatch;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Write normalized T1 image %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

matlabbatch{1}.spm.tools.shoot.norm.data.subj.deformation = cellstr(deformfile);
matlabbatch{1}.spm.tools.shoot.norm.data.subj.images = {SubjT1file};
matlabbatch{1}.spm.tools.shoot.norm.template = cellstr(strcat(spm('Dir'), fs, 'toolbox', fs, 'Shoot', fs, 'Template_4_IXI555_MNI152_GS.nii'));
matlabbatch{1}.spm.tools.shoot.norm.vox = [1 1 1];
matlabbatch{1}.spm.tools.shoot.norm.bb = [NaN NaN NaN; NaN NaN NaN];
matlabbatch{1}.spm.tools.shoot.norm.preserve = 0;
matlabbatch{1}.spm.tools.shoot.norm.fwhm = 0;

disp('*** writing normalised T1 image ...')
tic
spm_jobman('run', matlabbatch);
toc
clear tmp1 matlabbatch;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Quality control check for normalization %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

matlabbatch{1}.spm.util.checkreg.data = {spm_select('FPListRec', pwd, '^wT1*.*\.nii$');...
    spm_select('FPListRec', pwd, '^wc1.*\.nii$');...
    spm_select('FPListRec', pwd, '^wc2.*\.nii$');...
    spm_select('FPListRec', pwd, '^wc3.*\.nii$');...
    spm_select('FPListRec', pwd, '^mwc1.*\.nii$');...
    spm_select('FPListRec', pwd, '^mwc2.*\.nii$');...
    spm_select('FPListRec', pwd, '^mwc3.*\.nii$');...
    };
spm_jobman('run', matlabbatch);

printfname = strcat(pwd, fs, 'normalized_VBM_QC.pdf');
print(printfname, '-dpdf');
clear matlabbatch;
quit;

end

