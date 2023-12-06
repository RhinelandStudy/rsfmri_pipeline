"""
#Rhineland Study fMRI pre-processing pipeline
#rs_fmri_pipeline_v2: RestingState fMRI data processing using SPM12, FSL with nipype workflow engine.
"""
import os
import sys
from glob import glob

if os.path.exists('MANIFEST'): os.remove('MANIFEST')

def main(**extra_args):
    from setuptools import setup
    setup(name='fmri_pipeline',
        version='2.0.0',
        description='collection of scripts and routines to process rsFMRI data',
        long_description="""RhinelandStudy processing for rsfMRI RestingState scans.""",
        author='Hweeling Lee',
        author_email='Hwee-Ling.Lee@dzne.de',
        license='MIT',
        packages= ['fmri_pipeline', 'fmri_pipeline.spmscripts'],
        package_data = {'fmri_pipeline':
                ['FC_templates/*.csv',
                 'FC_templates/*.xlsx',
                 'FC_templates/*.nii',
                 'train_weights/*.RData',
                 'spmscripts/*.m',
                 'spmscripts/*.sh'
                ]},
        include_package_data=True,
        entry_points={ 'console_scripts': ['run_fmri_pipeline=fmri_pipeline.run_fmri_pipeline:main',
                                           'run_vbm_pipeline=fmri_pipeline.run_vbm_pipeline:main']},
        classifiers = [c.strip() for c in """\
            Development Status :: v2
            Intended Audience :: Developers
            Intended Audience :: Science/Research
            Operating System :: OS Independent
            Programming Language :: Python
            Topic :: Software Development
            """.splitlines() if len(c.split()) > 0],
          maintainer = 'RheinlandStudy MRI/RS-IT group, DZNE',
          maintainer_email = 'mohammad.shahid@dzne.de',
          install_requires=['nipype', 'pycrypto','psycopg2','pyxnat'],
          zip_safe=False,
          **extra_args
          )

if __name__ == "__main__":
    main()
