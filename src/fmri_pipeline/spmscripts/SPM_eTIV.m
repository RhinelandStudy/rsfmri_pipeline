function SPM_eTIV(seg8matfile)


fs = filesep;
clear tmpfile;

disp('****************************************************');
disp('INPUT SETTINGS....');
disp('****************************************************');
disp(strcat('Segmentation:', seg8matfile));
disp('***************************************************');

[fp,fn,ext] = fileparts(seg8matfile);

matlabbatch{1}.spm.util.tvol.matfiles = cellstr(seg8matfile);
matlabbatch{1}.spm.util.tvol.tmax = 3;
matlabbatch{1}.spm.util.tvol.mask = cellstr(strcat(spm('Dir'), fs, 'tpm', fs, 'mask_ICV.nii,1'));
matlabbatch{1}.spm.util.tvol.outf = strcat(fp,fs,'etiv');

tic
spm_jobman('run', matlabbatch);
toc
clear matlabbatch;

end

