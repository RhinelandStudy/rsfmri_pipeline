function calc_qualitycontrol(swufileb4ica,swufileafterica,greymatterfile,whitematterfile,csffile,realignmentfile,segmentfile)
fs = filesep;

% inputs
TR = 0.57;
radius = 50;
vsize_cutoff = 1.5;
intensity_normalization = 1000;

%swufileb4ica = '/rs_preproc/swuRestingState.nii';
%swufileafterica = '/ica_fix/filtered_func_data_clean_chdt.nii';
%greymatterfile = '/t1/wc1T1_RMS.nii';
%whitematterfile = '/t1/wc2T1_RMS.nii';
%csffile = '/t1/wc3T1_RMS.nii';
%realignmentfile = '/rs_preproc/rp_RestingState_00018.txt';
%segmentfile = 't1/T1_RMS_seg8.mat';

maskfile = strcat(spm('Dir'), '/tpm/BOLD_mask.nii');

% output
printfile1 =  'FD_carpet_b4ica.pdf';
printfile2 =  'FD_carpet_afterica.pdf';
savedfname =  'all_metrics.csv';

if strfind(swufileafterica,'gz')
    gunzip(swufileafterica,pwd);
    [fp,fn]= fileparts(swufileafterica);
    swufileafterica=fn;
else
    [fp,fn]= fileparts(swufileafterica);
    swufileafterica=fn;
end

% load mask (without skull), grey matter, white matter and csf
mask_info = spm_vol(maskfile);
mask = round(spm_read_vols(mask_info));
greymatter = spm_read_vols(spm_vol(greymatterfile));
whitematter = spm_read_vols(spm_vol(whitematterfile));
csf = spm_read_vols(spm_vol(csffile));

%%%%%%%%
% eTIV %
%%%%%%%%

load(segmentfile, 'volumes');
saved_info = volumes.litres;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate mean framewise displacement %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

motion1 = dlmread(realignmentfile);
motion = spm_detrend(motion1,1); % de-mean and detrend

[fn, filename, ext] = fileparts(realignmentfile);
if strncmp('FSL_motion', filename, 20)
    motion(:,[1 2 3]) = motion(:,[1 2 3]).*radius; % angle in radian by average head size = displacement in mm
    motion1(:,[1 2 3]) = motion1(:,[1 2 3]).*radius;
elseif strncmp('rp_RestingState_00018', filename, 20)
    motion(:,[4 5 6]) = motion(:,[4 5 6]).*radius;
    motion1(:,[4 5 6]) = motion1(:,[4 5 6]).*radius;
end

D = diff(motion,1,1); % 1st order derivative
D = [zeros(1,6); D];  % set first row to 0
FD = sum(abs(D),2); % framewise displacement a la Powers
% RMS = sqrt(mean(detrend(D).^2,2)); % root mean square for each column a la Van Dijk
saved_info = [saved_info, mean(FD)];

%%%%%%%%%%%%%%%%%%%
% spike detection %
%%%%%%%%%%%%%%%%%%%

tmp1 = [];

for ii = 1:size(motion,1)
    if sum(abs(motion(ii,:)) > vsize_cutoff) > 0
        tmp1 = [tmp1; 1];
    else
        tmp1 = [tmp1; 0];
    end
end

saved_info = [saved_info, sum(tmp1)];

%%%%%%%%%%%%%%%%%%%
% create grayplot %
%%%%%%%%%%%%%%%%%%%
% before ica
clear Fimginfo Fimgb4ica;
Fimginfo = spm_vol(swufileb4ica);
Fimgb4ica = spm_read_vols(Fimginfo);
Nt = size(Fimginfo,1);
Ni = Fimginfo(1).dim(1);
Nj = Fimginfo(1).dim(2);
Nk = Fimginfo(1).dim(3);

% remove linear and polynomial trends from data
clear F_2D X_deisgn betas;
X_design = [(1:Nt)' ((1:Nt).^2/(Nt^2))' ((1:Nt).^3/(Nt^3))'];
X_design = X_design - mean(X_design);
X_design = [X_design ones(Nt,1)];

F_2D = reshape(Fimgb4ica, Ni*Nj*Nk, Nt);
betas = X_design\F_2D';
F_detrendedb4ica = [F_2D' - X_design(:, 1:(end-1))*betas(1:(end-1),:)]';

% mask, mean and psc
clear F_masked F_mean F_masked_psc;
F_masked = F_detrendedb4ica(find(mask),:);
F_mean = mean(F_masked,2);
F_masked_psc = 100*(F_masked./repmat(F_mean, 1, Nt)) - 100;
F_masked_psc(isnan(F_masked_psc)) = 0;
F_psc_imgb4ica = zeros(Ni, Nj, Nk, Nt);
F_2D_pscb4ica = reshape(F_psc_imgb4ica, Ni*Nj*Nk, Nt);
F_2D_pscb4ica(find(mask),:) = F_masked_psc;
F_psc_imgb4ica = reshape(F_2D_pscb4ica, Ni, Nj, Nk, Nt);

% figure
GM_b4ica = F_2D_pscb4ica(find(greymatter),:);
WM_b4ica = F_2D_pscb4ica(find(whitematter),:);
CSF_b4ica = F_2D_pscb4ica(find(csf),:);
all_b4ica = [GM_b4ica; WM_b4ica; CSF_b4ica];

line1_pos = numel(find(greymatter));
line2_pos = numel(find(greymatter)) + numel(find(whitematter));

fontsizeL = 12; fontsizeM = 12; intensity_scale = [-6 6];

figure;
ax1 = subplot(5,1,[2:5]);
imagesc(ax1, all_b4ica);
colormap(gray);
caxis(intensity_scale);
title(ax1, 'carpet plot: before ICA', 'fontsize', fontsizeL);
ylabel(ax1, 'Voxels', 'fontsize', fontsizeM);
xlabel(ax1, 'fMRI Volumes', 'fontsize', fontsizeM);
hold on;

line([1 Nt], [line1_pos line1_pos], 'Color', 'b', 'LineWidth', 2);
line([1 Nt], [line2_pos line2_pos], 'Color', 'r', 'LineWidth', 2);
hold off;

ax3 = subplot(5,1,1);
plot(ax3, motion1, 'LineWidth', 0.1); grid;
hold on;
line([1 Nt], [vsize_cutoff vsize_cutoff], 'Color', 'r', 'LineWidth', 0.1);
line([1 Nt], [-vsize_cutoff -vsize_cutoff], 'Color', 'r', 'LineWidth', 0.1);
axis tight;
set(ax3, 'XTicklabel', []);
title(ax3, 'original motion in mm', 'fontsize', fontsizeL);
ylabel(ax3, 'mm', 'fontsize', fontsizeM);

print(printfile1, '-dpdf');

% after ica
clear Fimginfo Fimgb4ica;
Fimginfo = spm_vol(swufileafterica);
Fimgafterica = spm_read_vols(Fimginfo);
Nt = size(Fimginfo,1);
Ni = Fimginfo(1).dim(1);
Nj = Fimginfo(1).dim(2);
Nk = Fimginfo(1).dim(3);

% remove linear and polynomial trends from data
clear F_2D X_deisgn betas;
X_design = [(1:Nt)' ((1:Nt).^2/(Nt^2))' ((1:Nt).^3/(Nt^3))'];
X_design = X_design - mean(X_design);
X_design = [X_design ones(Nt,1)];

F_2D = reshape(Fimgafterica, Ni*Nj*Nk, Nt);
betas = X_design\F_2D';
F_detrendedafterica = [F_2D' - X_design(:, 1:(end-1))*betas(1:(end-1),:)]';

% mask, mean and psc
clear F_masked F_mean F_masked_psc;
F_masked = F_detrendedafterica(find(mask),:);
F_mean = mean(F_masked,2);
F_masked_psc = 100*(F_masked./repmat(F_mean, 1, Nt)) - 100;
F_masked_psc(isnan(F_masked_psc)) = 0;
F_psc_imgafterica = zeros(Ni, Nj, Nk, Nt);
F_2D_pscafterica = reshape(F_psc_imgafterica, Ni*Nj*Nk, Nt);
F_2D_pscafterica(find(mask),:) = F_masked_psc;
F_psc_imgb4ica = reshape(F_2D_pscafterica, Ni, Nj, Nk, Nt);

% figure
GM_afterica = F_2D_pscafterica(find(greymatter),:);
WM_afterica = F_2D_pscafterica(find(whitematter),:);
CSF_afterica = F_2D_pscafterica(find(csf),:);
all_afterica = [GM_afterica; WM_afterica; CSF_afterica];

line1_pos = numel(find(greymatter));
line2_pos = numel(find(greymatter)) + numel(find(whitematter));

fontsizeL = 12; fontsizeM = 12; intensity_scale = [-6 6];

figure;
ax1 = subplot(5,1,[2:5]);
imagesc(ax1, all_afterica);
colormap(gray);
caxis(intensity_scale);
title(ax1, 'carpet plot: after ICA', 'fontsize', fontsizeL);
ylabel(ax1, 'Voxels', 'fontsize', fontsizeM);
xlabel(ax1, 'fMRI Volumes', 'fontsize', fontsizeM);
hold on;

line([1 Nt], [line1_pos line1_pos], 'Color', 'b', 'LineWidth', 2);
line([1 Nt], [line2_pos line2_pos], 'Color', 'r', 'LineWidth', 2);
hold off;

ax3 = subplot(5,1,1);
plot(ax3, motion, 'LineWidth', 0.1); grid;
hold on;
line([1 Nt], [vsize_cutoff vsize_cutoff], 'Color', 'r', 'LineWidth', 0.1);
line([1 Nt], [-vsize_cutoff -vsize_cutoff], 'Color', 'r', 'LineWidth', 0.1);
axis tight;
set(ax3, 'XTicklabel', []);
title(ax3, 'original motion in mm', 'fontsize', fontsizeL);
ylabel(ax3, 'mm', 'fontsize', fontsizeM);


print(printfile2, '-dpdf');

%%%%%%%%%%%%%%%%%%%
% calculate DVARS %
%%%%%%%%%%%%%%%%%%%
% before ica
clear newepi Fimginfo Nt Ni Nj Nk;

% special functions

IntnstyScl = @(Y,md,scl) (Y./md).*scl;
IQRsd = @(x) (quantile(x, 0.75) - quantile(x, 0.25))./1.349;

epi = spm_read_vols(spm_vol(swufileb4ica));
newepi = mask .* epi;
Fimginfo = spm_vol(swufileb4ica);
Nt = size(Fimginfo,1);
Ni = Fimginfo(1).dim(1);
Nj = Fimginfo(1).dim(2);
Nk = Fimginfo(1).dim(3);
newepi = reshape(newepi, Ni*Nj*Nk, Nt);

% remove voxels of zeros/NaNs
clear nan_idx zero_idx idx n_voxels;
nan_idx = find(isnan(sum(newepi,2)));
zero_idx = find(sum(newepi,2)==0);
idx = 1:(Ni*Nj*Nk);
idx([nan_idx;zero_idx]) = [];
newepi([nan_idx;zero_idx],:) = [];
n_voxels = size(newepi,1);

% intensity normalization scale by 1000
clear md;
md = median(mean(newepi,2));
scl = 1000;
newepi = IntnstyScl(newepi,md,scl);
% center data
clear mvY dmeaner;
mvY = mean(newepi,2);
dmeaner = repmat(mvY, [1 Nt]);
newepi = newepi - dmeaner;
clear dmeaner;

clear DY DVARS stdDVARS Rob_S AC;
DY = diff(newepi,1,2);
DVARS = sqrt(sum(DY.^2)./size(newepi,1));

Rob_S = IQRsd(newepi');
AC    = zeros(1, n_voxels);
for iv=1:n_voxels
    AC(iv) = madicc(newepi(iv,1:end-1),newepi(iv,2:end));
end

ACf_idx = isnan(AC);
AC(ACf_idx) = [];
Rob_S(ACf_idx) = [];

stdDVARS = DVARS./(sqrt((sum(2*(1-AC).*(Rob_S.^2)))./n_voxels));

saved_info = [saved_info, mean(DVARS), mean(stdDVARS)];
% after ica
clear newepi Fimginfo Nt Ni Nj Nk;

epi = spm_read_vols(spm_vol(swufileafterica));
newepi = mask .* epi;
Fimginfo = spm_vol(swufileafterica);
Nt = size(Fimginfo,1);
Ni = Fimginfo(1).dim(1);
Nj = Fimginfo(1).dim(2);
Nk = Fimginfo(1).dim(3);
newepi = reshape(newepi, Ni*Nj*Nk, Nt);

% remove voxels of zeros/NaNs
clear nan_idx zero_idx idx n_voxels;
nan_idx = find(isnan(sum(newepi,2)));
zero_idx = find(sum(newepi,2)==0);
idx = 1:(Ni*Nj*Nk);
idx([nan_idx;zero_idx]) = [];
newepi([nan_idx;zero_idx],:) = [];
n_voxels = size(newepi,1);

% intensity normalization scale by 1000
clear md;
md = median(mean(newepi,2));
scl = 1000;
newepi = IntnstyScl(newepi,md,scl);
% center data
clear mvY dmeaner;
mvY = mean(newepi,2);
dmeaner = repmat(mvY, [1 Nt]);
newepi = newepi - dmeaner;
clear dmeaner;

clear DY DVARS stdDVARS Rob_S AC;
DY = diff(newepi,1,2);
DVARS = sqrt(sum(DY.^2)./size(newepi,1));

Rob_S = IQRsd(newepi');
AC    = zeros(1, n_voxels);
for iv=1:n_voxels
    AC(iv) = madicc(newepi(iv,1:end-1),newepi(iv,2:end));
end

ACf_idx = isnan(AC);
AC(ACf_idx) = [];
Rob_S(ACf_idx) = [];

stdDVARS = DVARS./(sqrt((sum(2*(1-AC).*(Rob_S.^2)))./n_voxels));

saved_info = [saved_info, mean(DVARS), mean(stdDVARS)];

%%%%%%%%%%%%%%%%%%
% calculate tSNR %
%%%%%%%%%%%%%%%%%%
% before ica
clear avg sd_tmp epi;
avg = zeros(size(mask));
sd_tmp = zeros(size(mask));

% read input BOLD data
epi = spm_read_vols(spm_vol(swufileb4ica));

for i = 1:size(epi,4)
    clear new_epi;
    new_epi = mask.*epi(:,:,:,i);
    avg = avg + (new_epi/size(epi,4));
end

for i = 1:size(epi,4)
    clear new_epi;
    new_epi = mask.*epi(:,:,:,i);
    sd_tmp = sd_tmp + (new_epi - avg).^2;
end
sd = sqrt(sd_tmp/size(epi,4));
snr = avg./sd;
% general tSNR
tSNR = nanmedian(snr,'all');
saved_info = [saved_info, tSNR];

% calculate tSNR after ica
clear avg sd_tmp epi sd snr tSNR;
avg = zeros(size(mask));
sd_tmp = zeros(size(mask));

% read input BOLD data
epi = spm_read_vols(spm_vol(swufileafterica));

for i = 1:size(epi,4)
    clear new_epi;
    new_epi = mask.*epi(:,:,:,i);
    avg = avg + (new_epi/size(epi,4));
end

for i = 1:size(epi,4)
    clear new_epi;
    new_epi = mask.*epi(:,:,:,i);
    sd_tmp = sd_tmp + (new_epi - avg).^2;
end
sd = sqrt(sd_tmp/size(epi,4));
snr = avg./sd;
% general tSNR
tSNR = nanmedian(snr,'all');
saved_info = [saved_info,tSNR];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate GSR-y before ica %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear n2_mask ghost signal gsr_y;
mask_info = spm_vol(maskfile);
mask = round(spm_read_vols(mask_info));
epi = spm_read_vols(spm_vol(swufileb4ica));
n2_mask = circshift(mask, size(mask,2)./2, 2); % y-direction
n2_mask = n2_mask .* (1 - mask);
n2_mask = n2_mask + 2 .* (1 - n2_mask - mask);
ghost = mean(epi(n2_mask==1)) - mean(epi(n2_mask == 2));
signal = median(epi(n2_mask==0));
gsr_y = ghost / signal;

saved_info = [saved_info, gsr_y];

t2 = array2table(saved_info,'VariableNames', {'gm_vol','wm_vol','csf_vol','FD', 'Nspikes','dvars_b4ica', 'std_dvars_b4ica','dvars_afterica', 'std_dvars_afterica','tSNR_b4ica', 'tSNR_afterica', 'gsry'});
writetable(t2, savedfname);

quit;

end

