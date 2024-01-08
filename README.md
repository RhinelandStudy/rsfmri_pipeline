# rsfmri_pipeline

This repository constains the Nipype wrapper for the  RestingState functional MRI processing pipeline using RestingState BOLD scans.

This pipeline corresponds to standard SPM process and includes a specifically trained FMRIB's ICA-based X-noisefier (FIX 1.06) classifier trained with 80 participants from the Rhineland Study. If you use this wrapper please cite:

Griffanti, L. et al. ICA-based artefact removal and accelerated fMRI acquisition for improved resting state network imaging. NeuroImage 95, 232–247 (2014).

```
@article{salimi2014automatic,
  title={Automatic denoising of functional MRI data: combining independent component analysis and hierarchical fusion of classifiers},
  author={Salimi-Khorshidi, Gholamreza and Douaud, Gwena{\"e}lle and Beckmann, Christian F and Glasser, Matthew F and Griffanti, Ludovica and Smith, Stephen M},
  journal={Neuroimage},
  volume={90},
  pages={449--468},
  year={2014},
  publisher={Elsevier}
}
```

The QC pipeline including the carpet plot based on:

Power, J. D. A simple but useful way to assess fMRI scan qualities. Neuroimage 154, 150–158 (2017).

```
@article{power2017simple,
  title={A simple but useful way to assess fMRI scan qualities},
  author={Power, Jonathan D},
  journal={Neuroimage},
  volume={154},
  pages={150--158},
  year={2017},
  publisher={Elsevier}
}
```


## Build docker image

```bash

docker build -t rsfmri_pipeline -f docker/Dockerfile .


```

## Or pull from docker hub

```bash
docker pull dznerheinlandstudie/rheinlandstudie:rsfmri_pipeline
```

## Run pipeline:

### Using docker
The pipeline can be run with docker by running the container as follows:


```bash

docker run --rm -v /path/to/input_scans:/input \
                 -v /path/to/work_folder:/work \
                 -v /path/to/output:/output \
        dznerheinlandstudie/rheinlandstudie:rsfmri_pipeline \
        run_fmri_pipeline \
        -s /input \
        --subjects test_subject_01 \
        -w /work \
        -o /output \ 
        -p 2 -m 2 \
        -d /path/to/TrainWts.RData \
        --fixth 10

```

The command line options are described briefly if the pipeline is started with only ```-h``` option.

### Using Singulraity

The pipeline can be run with Singularity by running the singularity image as follows:

```bash


singularity build rsfmri_pipeline.sif docker://dznerheinlandstudie/rheinlandstudie:rsfmri_pipeline
```

When the singularit image is created, then it can be run as follows:

```bash

singularity run -e -B /path/to/input_scans:/input \
                   -B /path/to/work:/work \
                   -B /path/to/output:/output \
                   -B /path/to/trainwts:/trainwts \
       rsfmri_pipeline.sif 
    run_fmri_pipeline \
     -s /input \
     -w /work \
     -o /output \
     -p 2 \
     -m 2 \
     -d /trainwts/TrainWts.RData \
     --fixth 10 \
     --subjects test_subject_01

```
