% Copyright 2023 Population Health Sciences, German Center for Neurodegenerative Diseases (DZNE)
%
%    Licensed under the Apache License, Version 2.0 (the "License");
%    you may not use this file except in compliance with the License.
%    You may obtain a copy of the License at
%
%        http://www.apache.org/licenses/LICENSE-2.0
%
%    Unless required by applicable law or agreed to in writing, software
%    distributed under the License is distributed on an "AS IS" BASIS,
%    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
%    See the License for the specific language governing permissions and
%    limitations under the License.

function SPM_rsfmri_glm_200_100(swurestfile,rprestfile,wc2file,wc3file)

% Date created: 2021-03-12, Author: Hweeling Lee
% This script creates 2 GLM.
% 1st GLM is a simple GLM without any regressors (just global signal regressor). This is used to extract WM and CSF signal intensities to "correct" for physiological responses in the 2nd GLM.
% 2nd GLM contains discrete cosine basis function that represent low frequency fluctuations from 0.0087 - 0.1 Hz, 6 motion parameters, WM and CSF signal intensities.
% inputs: swuRestingState.nii, wc1:c3.nii, rp*.txt, rSchaefer2018_400Parcels_17Networks.nii, Schaefer2018_400Parcels_17Networks_order_FSLMNI152_2mm.Centroid_RAS.csv
% current output of the GLM is a 400 ROI functional connectivity matrix. If research question is intereested in a specific ROI, a 3rd GLM should be created with the signal intensity of that ROI.
% output to be saved: FC_Schaefer_17Networks_400p_simple.csv


fs = filesep;
n_sess = 1;
TR = 0.57;
TE = 0.3;
ndummy = 17;
roi_size = 8;

dir_parcellations=strcat(spm('Dir'), fs, 'fctemplates');

opts = detectImportOptions(strcat(dir_parcellations, fs, 'Schaefer2018_200Parcels_17Networks_order_FSLMNI152_2mm.Centroid_RAS.csv'));
t1 = readtable(strcat(dir_parcellations, fs, 'Schaefer2018_200Parcels_17Networks_order_FSLMNI152_2mm.Centroid_RAS.csv'), opts);
opts = detectImportOptions(strcat(dir_parcellations, fs, 'Schaefer2018_100Parcels_17Networks_order_FSLMNI152_2mm.Centroid_RAS.csv'));
t2 = readtable(strcat(dir_parcellations, fs, 'Schaefer2018_100Parcels_17Networks_order_FSLMNI152_2mm.Centroid_RAS.csv'), opts);

disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
display('INPUTS:');
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');

clear SubjFFX1 SubjFFX2 SubjRSdir;


SubjFFX1 = strcat(pwd, fs, 'FFX1');
SubjFFX2 = strcat(pwd, fs, 'FFX2');
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
matlabbatch{2}.spm.util.voi.roi{1}.mask.image = cellstr(spm_select('FPList', fp3, wc2file));

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

Move_fname = strcat(SubjFFX2, fs, '_RS_GLM_Para.mat');
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

%% 200 Parcels

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% signal extraction for each ROI %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for roi = 1:size(t1,1)
    
    matlabbatch{roi}.spm.util.voi.spmmat = cellstr(strcat(SubjFFX2, fs, 'SPM.mat'));
    matlabbatch{roi}.spm.util.voi.adjust = 1;
    matlabbatch{roi}.spm.util.voi.session = 1;
    matlabbatch{roi}.spm.util.voi.name = t1.ROIName{roi};
    matlabbatch{roi}.spm.util.voi.roi{1}.label.image = cellstr(spm_select('ExtFPList', dir_parcellations, '^rSchaefer2018_200Parcels_17Networks.*\.nii$'));
    matlabbatch{roi}.spm.util.voi.roi{1}.label.list = roi;
    matlabbatch{roi}.spm.util.voi.roi{2}.mask.image = cellstr(spm_select('FPList', SubjFFX2, '^mask.*\.nii$'));
    matlabbatch{roi}.spm.util.voi.expression = 'i1 & i2';
    
end

spm_jobman('run', matlabbatch);
clear matlabbatch;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% checking for any missing ROI %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear nVOI varnames FC_corr tmp1 tmp2 tmp3 tmp4;
nVOI = dir(strcat(SubjFFX2, fs, 'VOI_17Networks*.mat'));

if size(nVOI,1) < 200

    for ii = 1:size(nVOI,1)
        tmp1{ii,1} = nVOI(ii).name(5:end-6);
    end

    tmp2 = find(~ismember(t1.ROIName, tmp1));

    for jj = 1:size(tmp2,1)
        matlabbatch{jj}.spm.util.voi.spmmat = cellstr(strcat(SubjFFX2, fs, 'SPM.mat'));
        matlabbatch{jj}.spm.util.voi.adjust = 1;
        matlabbatch{jj}.spm.util.voi.session = 1;
        matlabbatch{jj}.spm.util.voi.name = t1.ROIName{tmp2(jj,1)};
        matlabbatch{jj}.spm.util.voi.roi{1}.sphere.centre     = [t1.R(tmp2(jj,1)) t1.A(tmp2(jj,1)) t1.S(tmp2(jj,1))];
        matlabbatch{jj}.spm.util.voi.roi{1}.sphere.radius     = 8;
        matlabbatch{jj}.spm.util.voi.roi{1}.sphere.move.fixed = 1;
        matlabbatch{jj}.spm.util.voi.roi{2}.mask.image = cellstr(spm_select('FPList', SubjFFX2, '^mask.*\.nii$'));
        matlabbatch{jj}.spm.util.voi.expression = 'i1 & i2';

    end
    spm_jobman('run', matlabbatch);
    clear matlabbatch tmp1 tmp2;

end

%%%%%%%%%%%%%%%%%%%%%%%%%
% FC correlation matrix %
%%%%%%%%%%%%%%%%%%%%%%%%%

clear nVOI varnames FC_corr tmp1 tmp2 tmp3 tmp4;

nVOI = dir(strcat(SubjFFX2, fs, 'VOI_17Networks*.mat'));

for ii = 1:size(nVOI,1)
    clear tmp3;
    tmp3 = regexp(nVOI(ii).name, {'_'}, 'split');
    tmp4(ii,:) = [ii tmp3{1}(:,[3:6])];
end

[G_network TID] = findgroups(tmp4(:,3));

for ii = 1:size(nVOI,1)
    allnetworks(ii,:) = [tmp4(ii,:) G_network(ii,:)];
end

count = 1;
clear tmp_raw1 varnames;
for gg = 1:size(TID,1)
    clear int_network;
    int_network = allnetworks(find(cell2mat(allnetworks(:,6))==gg),:);
    
    for ii = 1:size(int_network,1)
        varnames{1,count} = nVOI(cell2mat(int_network(ii,1))).name(16:end-6);
        clear xY;
        load(strcat(nVOI(cell2mat(int_network(ii,1))).folder, fs, nVOI(cell2mat(int_network(ii,1))).name));
        tmp_raw1(:,count) =  mean(xY.y,2);
        count = count + 1;
    end
end

FC_corr = corr(tmp_raw1, tmp_raw1);

clear t3;
t3 = array2table(FC_corr, 'VariableNames', varnames);
writetable(t3, strcat(SubjFFX2, fs, 'FC_Schaefer_17Networks_200p_simple.csv'));

%% for 100 parcels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% signal extraction for each ROI %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cd(SubjFFX2);
delete VOI_17Networks*;

for roi = 1:size(t2,1)
    
    matlabbatch{roi}.spm.util.voi.spmmat = cellstr(strcat(SubjFFX2, fs, 'SPM.mat'));
    matlabbatch{roi}.spm.util.voi.adjust = 1;
    matlabbatch{roi}.spm.util.voi.session = 1;
    matlabbatch{roi}.spm.util.voi.name = t2.ROIName{roi};
    matlabbatch{roi}.spm.util.voi.roi{1}.label.image = cellstr(spm_select('ExtFPList', dir_parcellations, '^rSchaefer2018_100Parcels_17Networks.*\.nii$'));
    matlabbatch{roi}.spm.util.voi.roi{1}.label.list = roi;
    matlabbatch{roi}.spm.util.voi.roi{2}.mask.image = cellstr(spm_select('FPList', SubjFFX2, '^mask.*\.nii$'));
    matlabbatch{roi}.spm.util.voi.expression = 'i1 & i2';
    
end

spm_jobman('run', matlabbatch);
clear matlabbatch;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% checking for any missing ROI %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear nVOI varnames FC_corr tmp1 tmp2 tmp3 tmp4 G_network TID allnetworks;
nVOI = dir(strcat(SubjFFX2, fs, 'VOI_17Networks*.mat'));

if size(nVOI,1) < 100

    for ii = 1:size(nVOI,1)
        tmp1{ii,1} = nVOI(ii).name(5:end-6);
    end

    tmp2 = find(~ismember(t2.ROIName, tmp1));

    for jj = 1:size(tmp2,1)
        matlabbatch{jj}.spm.util.voi.spmmat = cellstr(strcat(SubjFFX2, fs, 'SPM.mat'));
        matlabbatch{jj}.spm.util.voi.adjust = 1;
        matlabbatch{jj}.spm.util.voi.session = 1;
        matlabbatch{jj}.spm.util.voi.name = t2.ROIName{tmp2(jj,1)};
        matlabbatch{jj}.spm.util.voi.roi{1}.sphere.centre     = [t2.R(tmp2(jj,1)) t2.A(tmp2(jj,1)) t2.S(tmp2(jj,1))];
        matlabbatch{jj}.spm.util.voi.roi{1}.sphere.radius     = 8;
        matlabbatch{jj}.spm.util.voi.roi{1}.sphere.move.fixed = 1;
        matlabbatch{jj}.spm.util.voi.roi{2}.mask.image = cellstr(spm_select('FPList', SubjFFX2, '^mask.*\.nii$'));
        matlabbatch{jj}.spm.util.voi.expression = 'i1 & i2';

    end
    spm_jobman('run', matlabbatch);
    clear matlabbatch tmp1 tmp2;

end

%%%%%%%%%%%%%%%%%%%%%%%%%
% FC correlation matrix %
%%%%%%%%%%%%%%%%%%%%%%%%%

clear nVOI varnames FC_corr tmp1 tmp2 tmp3 tmp4 G_network TID allnetworks;

nVOI = dir(strcat(SubjFFX2, fs, 'VOI_17Networks*.mat'));

for ii = 1:size(nVOI,1)
    clear tmp3;
    tmp3 = regexp(nVOI(ii).name, {'_'}, 'split');
    tmp4(ii,:) = [ii tmp3{1}(:,[3:6])];
end

[G_network TID] = findgroups(tmp4(:,3));

for ii = 1:size(nVOI,1)
    allnetworks(ii,:) = [tmp4(ii,:) G_network(ii,:)];
end

count = 1;
clear tmp_raw1 varnames;
for gg = 1:size(TID,1)
    clear int_network;
    int_network = allnetworks(find(cell2mat(allnetworks(:,6))==gg),:);
    
    for ii = 1:size(int_network,1)
        varnames{1,count} = nVOI(cell2mat(int_network(ii,1))).name(16:end-6);
        clear xY;
        load(strcat(nVOI(cell2mat(int_network(ii,1))).folder, fs, nVOI(cell2mat(int_network(ii,1))).name));
        tmp_raw1(:,count) =  mean(xY.y,2);
        count = count + 1;
    end
end

FC_corr = corr(tmp_raw1, tmp_raw1);

clear t3;
t3 = array2table(FC_corr, 'VariableNames', varnames);
writetable(t3, strcat(SubjFFX2, fs, 'FC_Schaefer_17Networks_100p_simple.csv'));
end
