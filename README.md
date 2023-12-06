# rsfmri_pipeline
RestingState functional MRI processing pipeline using RestingState BOLD scans.


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
        -p 2 -t 2 -m 2 \
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
                   -B /path/to/inputdata:/input \
                   -B /path/to/work:/work \
                   -B /path/to/output:/output \
                   -B /path/to/trainwts:/trainwts \
       rsfmri_pipeline.sif 
    run_fmri_pipeline \
     -p 2 \
     -m 2 \
     -d /trainwts/TrainWts.RData \
     --fixth 10 \
     --subjects test_subject_01

```
