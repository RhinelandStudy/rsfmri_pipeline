function SPM_rsfmri_glm_Orbitofrontal(swurestfile,rprestfile,wc2file,wc3file)

% Date created: 2021-03-12, Author: Hweeling Lee
% This script creates 2 GLM.
% 1st GLM is a simple GLM without any regressors (just global signal regressor). This is used to extract WM and CSF signal intensities to "correct" for physiological responses in the 2nd GLM.
% 2nd GLM contains discrete cosine basis function that represent low frequency fluctuations from 0.0087 - 0.1 Hz, 6 motion parameters, WM and CSF signal intensities.
% inputs: swuRestingState.nii, wc1:c3.nii, rp*.txt, rSchaefer2018_400Parcels_17Networks.nii, Schaefer2018_400Parcels_17Networks_order_FSLMNI152_2mm.Centroid_RAS.csv
% current output of the GLM is a 400 ROI functional connectivity matrix. If research question is intereested in a specific ROI, a 3rd GLM should be created with the signal intensity of that ROI.
% output to be saved: FC_Schaefer_17Networks_400p_simple.csv

fs = filesep;

dir_parcellations=strcat(spm('Dir'), fs, 'fctemplates');

n_sess = 1;
TR = 0.57;
TE = 0.3;

disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
disp('INPUTS:');
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');

clear SubjFFX1 SubjFFX2 SubjRSdir;


SubjFFX1 = strcat(pwd, fs, 'FFX1');
SubjFFX2 = strcat(pwd, fs, 'FFX2');
SubjFFX3 = strcat(pwd, fs, 'L_Orbitofrontal');
SubjFFX4 = strcat(pwd, fs, 'R_Orbitofrontal');
SubjFFX5 = strcat(pwd, fs, 'Orbitofrontal');
SubjfMRI = strcat(pwd, fs, 'swuRestingState');

if strfind(swurestfile,'gz')
    copyfile(swurestfile, 'swuRestingState.nii.gz');
    swurestfile='swuRestingState.nii.gz';
    gunzip(swurestfile,pwd);
    swurestfile='swuRestingState.nii';
else
    copyfile(swurestfile, 'swuRestingState.nii');
    swurestfile='swuRestingState.nii';
end


if isfolder(SubjfMRI) == 0
    mkdir(SubjfMRI);
end

if isfolder(SubjFFX1) == 0
    mkdir(SubjFFX1);
end
if isfolder(SubjFFX2) == 0
    mkdir(SubjFFX2);
end

if isfolder(SubjFFX3) == 0
    mkdir(SubjFFX3);
end
if isfolder(SubjFFX4) == 0
    mkdir(SubjFFX4);
end
if isfolder(SubjFFX5) == 0
    mkdir(SubjFFX5);
end


%%%%%%%%%%%%%%%%%%%%
% convert 4D to 3D %
%%%%%%%%%%%%%%%%%%%%
[fp,fn] = fileparts(swurestfile);
swurestfile=fn;

matlabbatch{1}.spm.util.split.vol = cellstr(spm_select('ExtFPList', fp, swurestfile ,1));
matlabbatch{1}.spm.util.split.outdir = cellstr(SubjfMRI);

disp('*** converting 4D to 3D volumes ...')
tic
spm_jobman('run', matlabbatch);
toc
clear matlabbatch fp fn;

%%%%%%%%%%%%%%%%%
% 1st level GLM %
%%%%%%%%%%%%%%%%%

% load files
P = spm_select('ExtFPList',SubjfMRI, '^swuRestingState.*\.nii$');
n_vols = size(P,1);

%%%%%%%%%%%%%%%%%%%%%%%
% fmri specifications %
%%%%%%%%%%%%%%%%%%%%%%%
matlabbatch{1}.spm.stats.fmri_spec.dir = cellstr(SubjFFX1);
matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'secs';
matlabbatch{1}.spm.stats.fmri_spec.timing.RT = TR;
matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = 16;
matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = 8;

matlabbatch{1}.spm.stats.fmri_spec.sess(1).scans = cellstr(P);
matlabbatch{1}.spm.stats.fmri_spec.sess(1).multi = {''};
matlabbatch{1}.spm.stats.fmri_spec.sess(1).regress = struct('name', {}, 'val', {});
matlabbatch{1}.spm.stats.fmri_spec.sess(1).multi_reg = {''};
matlabbatch{1}.spm.stats.fmri_spec.sess(1).hpf = 128;

matlabbatch{1}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [1,0];
matlabbatch{1}.spm.stats.fmri_spec.volt = 1;
matlabbatch{1}.spm.stats.fmri_spec.global = 'None';
matlabbatch{1}.spm.stats.fmri_spec.mthresh  = 0.8;
matlabbatch{1}.spm.stats.fmri_spec.mask = {''};
matlabbatch{1}.spm.stats.fmri_spec.cvi = 'AR(1)';

spm_jobman('run', matlabbatch);
clear matlabbatch;

%%%%%%%%%%%%%%
% Estimation %
%%%%%%%%%%%%%%

matlabbatch{1}.spm.stats.fmri_est.spmmat = cellstr(spm_select('FPList', SubjFFX1, '/*SPM.mat'));
matlabbatch{1}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{1}.spm.stats.fmri_est.method.Classical = 1;

spm_jobman('run', matlabbatch);
clear matlabbatch;

%%%%%%%%%%%%%%%%%%
% ROI extraction %
%%%%%%%%%%%%%%%%%%
[fp3, fn] = fileparts(wc3file);
wc3file=fn;
clear fn;

[fp2, fn] = fileparts(wc2file);
wc2file=fn;
clear fn;

matlabbatch{1}.spm.util.voi.spmmat = cellstr(spm_select('FPList', SubjFFX1, '/*SPM.mat'));
matlabbatch{1}.spm.util.voi.adjust = NaN;
matlabbatch{1}.spm.util.voi.session = 1;
matlabbatch{1}.spm.util.voi.name = 'Whole_CSF';
matlabbatch{1}.spm.util.voi.roi{1}.mask.image = cellstr(spm_select('FPList', fp3, wc3file));
matlabbatch{1}.spm.util.voi.roi{1}.mask.threshold = 0.8;
matlabbatch{1}.spm.util.voi.roi{2}.mask.image = cellstr(spm_select('FPList', SubjFFX1, '^mask.*\.nii$'));
matlabbatch{1}.spm.util.voi.expression = 'i1 & i2';

matlabbatch{2} = matlabbatch{1};
matlabbatch{2}.spm.util.voi.name = 'Whole_WM';
matlabbatch{2}.spm.util.voi.roi{1}.mask.image = cellstr(spm_select('FPList', fp2, wc2file));

spm_jobman('run',matlabbatch);
clear matlabbatch;

% resting state regressors and motion regressors

tmp1 = [];
tmp2 = [];
tmp3 = [];
tmp4 = [];

UL       = 0.1; %Hz
LL       = 1/128; %Hz

tmp1 = spm_dctmtx(n_vols, n_vols);

% Calculate lower limit components & remove
LL_com   = fix(2*(n_vols*TR)*LL+1);
tmp1(:,1:LL_com) = [];

% Calculate upper limit components & remove
UL_com   = fix(2*(n_vols*TR)/(1/UL)+1);
tmp1(:,UL_com:end) = [];

[~,n_cols] = size(tmp1);

tmp2 = dlmread(rprestfile);

t = load(strcat(SubjFFX1, fs, 'VOI_Whole_CSF_1.mat'));
tmp3 = t.Y;
clear t;

t = load(strcat(SubjFFX1, fs, 'VOI_Whole_WM_1.mat'));
tmp4 = t.Y;
clear t;
R = [tmp1, tmp2, tmp3, tmp4];

Move_fname = strcat(SubjFFX2, fs, 'RS_GLM_Para.mat');
save(Move_fname, 'R');

disp('**** running stats.fmri_spec...');

%%%%%%%%%%%%%%%%%%%%%%%
% fmri specifications %
%%%%%%%%%%%%%%%%%%%%%%%
matlabbatch{1}.spm.stats.fmri_spec.dir = cellstr(SubjFFX2);
matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'secs';
matlabbatch{1}.spm.stats.fmri_spec.timing.RT = TR;
matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = 16;
matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = 8;

matlabbatch{1}.spm.stats.fmri_spec.sess(1).scans = cellstr(P);
matlabbatch{1}.spm.stats.fmri_spec.sess(1).multi = {''};
matlabbatch{1}.spm.stats.fmri_spec.sess(1).regress = struct('name', {}, 'val', {});
matlabbatch{1}.spm.stats.fmri_spec.sess(1).multi_reg = cellstr(spm_select('FPList', SubjFFX2, '/*_Para.mat'));
matlabbatch{1}.spm.stats.fmri_spec.sess(1).hpf = 128;

matlabbatch{1}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0,0];
matlabbatch{1}.spm.stats.fmri_spec.volt = 1;
matlabbatch{1}.spm.stats.fmri_spec.global = 'None';
matlabbatch{1}.spm.stats.fmri_spec.mthresh  = 0.8;
matlabbatch{1}.spm.stats.fmri_spec.mask = {''};
matlabbatch{1}.spm.stats.fmri_spec.cvi = 'AR(1)';

spm_jobman('run', matlabbatch);
clear matlabbatch;

%%%%%%%%%%%%%%
% Estimation %
%%%%%%%%%%%%%%

matlabbatch{1}.spm.stats.fmri_est.spmmat = cellstr(spm_select('FPList', SubjFFX2, '/*SPM.mat'));
matlabbatch{1}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{1}.spm.stats.fmri_est.method.Classical = 1;

spm_jobman('run', matlabbatch);
clear matlabbatch;

%%%%%%%%%%%%%%%%%%%%%
% Contrasts manager %
%%%%%%%%%%%%%%%%%%%%%

matlabbatch{1}.spm.stats.con.spmmat = cellstr(spm_select('FPList', SubjFFX2, '/*SPM.mat'));
matlabbatch{1}.spm.stats.con.delete = 1;

% F-contrast
matlabbatch{1}.spm.stats.con.consess{1}.fcon.name = 'EOF';
matlabbatch{1}.spm.stats.con.consess{1}.fcon.weights = eye(n_cols);
matlabbatch{1}.spm.stats.con.consess{1}.fcon.sessrep = 'none';

spm_jobman('run', matlabbatch);
clear matlabbatch;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% signal extraction for Orbitofrontal %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

matlabbatch{1}.spm.util.voi.spmmat = cellstr(strcat(SubjFFX2, fs, 'SPM.mat'));
matlabbatch{1}.spm.util.voi.adjust = 1;
matlabbatch{1}.spm.util.voi.session = 1;
matlabbatch{1}.spm.util.voi.name = 'L_Orbitofrontal';
matlabbatch{1}.spm.util.voi.roi{1}.label.image = cellstr(spm_select('ExtFPList', dir_parcellations, '^rAAL3v1_2_4mm.*\.nii$'));
matlabbatch{1}.spm.util.voi.roi{1}.label.list = [27; 29; 31];
matlabbatch{1}.spm.util.voi.roi{2}.mask.image = cellstr(spm_select('FPList', SubjFFX2, '^mask.*\.nii$'));
matlabbatch{1}.spm.util.voi.expression = 'i1 & i2';

matlabbatch{2} = matlabbatch{1};
matlabbatch{2}.spm.util.voi.name = 'R_Orbitofrontal';
matlabbatch{2}.spm.util.voi.roi{1}.label.list = [28; 30; 32];

spm_jobman('run', matlabbatch);
clear matlabbatch;

clear tmp1 tmp2 tmp3 tmp4;

t = load(strcat(SubjFFX2, fs, 'VOI_L_Orbitofrontal_1.mat'));
tmp1 = t.Y;
clear t;

tmp2 = dlmread(rprestfile);

t = load(strcat(SubjFFX1, fs, 'VOI_Whole_CSF_1.mat'));
tmp3 = t.Y;
clear t;

t = load(strcat(SubjFFX1, fs, 'VOI_Whole_WM_1.mat'));
tmp4 = t.Y;
clear t;

R = [tmp1, tmp2, tmp3, tmp4];

Move_fname = strcat(SubjFFX3, fs, 'RS_GLM_Para.mat');
save(Move_fname, 'R');

clear tmp1;

t = load(strcat(SubjFFX2, fs, 'VOI_R_Orbitofrontal_1.mat'));
tmp1 = t.Y;
clear t;
R = [tmp1, tmp2, tmp3, tmp4];

Move_fname = strcat(SubjFFX4, fs, 'RS_GLM_Para.mat');
save(Move_fname, 'R');

clear tmp1;

t = load(strcat(SubjFFX2, fs, 'VOI_L_Orbitofrontal_1.mat'));
tmp1 = t.Y;
clear t;
t = load(strcat(SubjFFX2, fs, 'VOI_R_Orbitofrontal_1.mat'));
tmp1(:,2) = t.Y;
clear t;

R = [mean(tmp1,2), tmp2, tmp3, tmp4];

Move_fname = strcat(SubjFFX5, fs, 'RS_GLM_Para.mat');
save(Move_fname, 'R');

%%%%%%%%%%%%%%%%%%%%%%%
% fmri specifications %
%%%%%%%%%%%%%%%%%%%%%%%
matlabbatch{1}.spm.stats.fmri_spec.dir = cellstr(SubjFFX3);
matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'secs';
matlabbatch{1}.spm.stats.fmri_spec.timing.RT = TR;
matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = 16;
matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = 8;

matlabbatch{1}.spm.stats.fmri_spec.sess(1).scans = cellstr(P);
matlabbatch{1}.spm.stats.fmri_spec.sess(1).multi = {''};
matlabbatch{1}.spm.stats.fmri_spec.sess(1).regress = struct('name', {}, 'val', {});
matlabbatch{1}.spm.stats.fmri_spec.sess(1).multi_reg = cellstr(spm_select('FPList', SubjFFX3, '/*_Para.mat'));
matlabbatch{1}.spm.stats.fmri_spec.sess(1).hpf = 128;

matlabbatch{1}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0,0];
matlabbatch{1}.spm.stats.fmri_spec.volt = 1;
matlabbatch{1}.spm.stats.fmri_spec.global = 'None';
matlabbatch{1}.spm.stats.fmri_spec.mthresh  = 0.8;
matlabbatch{1}.spm.stats.fmri_spec.mask = {''};
matlabbatch{1}.spm.stats.fmri_spec.cvi = 'AR(1)';

matlabbatch{2} = matlabbatch{1};
matlabbatch{2}.spm.stats.fmri_spec.dir = cellstr(SubjFFX4);
matlabbatch{2}.spm.stats.fmri_spec.sess(1).multi_reg = cellstr(spm_select('FPList', SubjFFX4, '/*_Para.mat'));

matlabbatch{3} = matlabbatch{1};
matlabbatch{3}.spm.stats.fmri_spec.dir = cellstr(SubjFFX5);
matlabbatch{3}.spm.stats.fmri_spec.sess(1).multi_reg = cellstr(spm_select('FPList', SubjFFX5, '/*_Para.mat'));

spm_jobman('run', matlabbatch);
clear matlabbatch;

%%%%%%%%%%%%%%
% Estimation %
%%%%%%%%%%%%%%

matlabbatch{1}.spm.stats.fmri_est.spmmat = cellstr(spm_select('FPList', SubjFFX3, '/*SPM.mat'));
matlabbatch{1}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{1}.spm.stats.fmri_est.method.Classical = 1;

matlabbatch{2} = matlabbatch{1};
matlabbatch{2}.spm.stats.fmri_est.spmmat = cellstr(spm_select('FPList', SubjFFX4, '/*SPM.mat'));

matlabbatch{3} = matlabbatch{1};
matlabbatch{3}.spm.stats.fmri_est.spmmat = cellstr(spm_select('FPList', SubjFFX5, '/*SPM.mat'));

spm_jobman('run', matlabbatch);
clear matlabbatch;

%%%%%%%%%%%%%%%%%%%%%
% Contrasts manager %
%%%%%%%%%%%%%%%%%%%%%

matlabbatch{1}.spm.stats.con.spmmat = cellstr(spm_select('FPList', SubjFFX3, '/*SPM.mat'));
matlabbatch{1}.spm.stats.con.delete = 1;

% T-contrast
matlabbatch{1}.spm.stats.con.consess{1}.tcon.name = '+ve';
matlabbatch{1}.spm.stats.con.consess{1}.tcon.weights = [1 zeros(1,9)];
matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep = 'none';

matlabbatch{2} = matlabbatch{1};
matlabbatch{2}.spm.stats.con.spmmat = cellstr(spm_select('FPList', SubjFFX4, '/*SPM.mat'));

matlabbatch{3} = matlabbatch{1};
matlabbatch{3}.spm.stats.con.spmmat = cellstr(spm_select('FPList', SubjFFX5, '/*SPM.mat'));

spm_jobman('run', matlabbatch);
clear matlabbatch;

end
