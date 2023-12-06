#!/bin/bash

for filename in *.nii; do
  fname=`$FSLDIR/bin/remove_ext ${filename}`
mkdir -p ${fname}_hp.ica
  ${FSLDIR}/bin/fslmaths ${fname} -Tmean ${fname}_hp.ica/mean_func
  ${FSLDIR}/bin/bet ${fname}_hp.ica/mean_func ${fname}_hp.ica/brain -f 0.4 -g 0 -n -m -t
  ${FSLDIR}/bin/fslmaths ${fname} -sub ${fname}_hp.ica/mean_func -bptf 187.1 -1 -add ${fname}_hp.ica/mean_func ${fname}_hp.ica/filtered_func_data
  ${FSLDIR}/bin/melodic -i ${fname}_hp.ica/filtered_func_data -o ${fname}_hp.ica/filtered_func_data.ica -m ${fname}_hp.ica/brain_mask.nii.gz --tr=0.57 -a concat -d -250 --nobet --Oall
mkdir -p ${fname}_hp.ica/reg
  ${FSLDIR}/bin/fslroi ${fname}_hp.ica/filtered_func_data ${fname}_hp.ica/reg/example_func 526 1
cp /home/zengw/Documents/FMRI/Newpip/FIX/Training/ICA/ICA/highres/highres.nii.gz ${fname}_hp.ica/reg/
  ${FSLDIR}/bin/flirt -in ${fname}_hp.ica/reg/highres -ref ${fname}_hp.ica/reg/example_func -omat ${fname}_hp.ica/reg/highres2example_func.mat
mkdir -p ${fname}_hp.ica/mc
cp /home/zengw/Documents/FMRI/Newpip/FIX/Training/ICA/ICA/motion/${fname}.txt ${fname}_hp.ica/mc/prefiltered_func_data_mcf.par
cp ${fname}_hp.ica/filtered_func_data.ica/mask.nii.gz ${fname}_hp.ica/
done
