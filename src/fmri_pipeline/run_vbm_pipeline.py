#!/usr/bin/env python

"""
Copyright 2023 Population Health Sciences, German Center for Neurodegenerative Diseases (DZNE)
Licensed under the Apache License, Version 2.0 (the "License"); 
you may not use this file except in compliance with the License. 
You may obtain a copy of the License at  http://www.apache.org/licenses/LICENSE-2.0 
Unless required by applicable law or agreed to in writing, software distributed under the 
License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and limitations under the License.
"""

from .vbm_pipeline import create_preprocess_wf

import nipype.pipeline.engine as pe
from nipype import config, logging
import logging as lgng

import os, sys, glob
import argparse
from itertools import chain

   
    
def get_pp_wf(scans_dir, work_dir, 
                         outputdir, subject_ids,
                         spm, mcr, wfname='preprocess_workflow'):
    
    ppwf = pe.Workflow(name=wfname)
    
    ppwf = create_preprocess_wf(scans_dir=scans_dir,
                              work_dir=work_dir,
                              subj_ids=subject_ids,
                              spm_exe_script=spm,
                              mcr_path=mcr,
                              name=wfname)
    ppwf.inputs.inputnode.outputdir = outputdir
    return ppwf

    
def main():
    """
    Command line wrapper for preprocessing data
    """
    parser = argparse.ArgumentParser(description='Run fMRI preprocess pipeline for RestingState '\
                                     '(bold) data.',
                                     epilog='Example-1: {prog} -s '\
                                     '~/data/scans -w '\
                                     '~/data/work  -o '\
                                     '~/data/outputdir -p 2 -t 2 '\
                                     '--subjects subj1 subj2 '\
                                     '\nExample-2: (All subjects): {prog} --scansdir ~/data/scans -o '\
                                     ' ~/data/outputdir -w ~/data/workdir -p [NProc] '\
                                     '\n\n'
                                     .format(prog=os.path.basename\
                                             (sys.argv[0])),\
                                     formatter_class=argparse.\
                                     RawTextHelpFormatter)

    parser.add_argument('-s', '--scansdir', help='Scans directory where scans'\
                        ' for each subjects (saved in image-uuid sub folders) are found.',required=True)
    parser.add_argument('-w', '--workdir', help='Work directory where data' \
                        ' is processed (intermediate files).', required=True)
    parser.add_argument('-o', '--outputdir', help='Output directory where ' \
                        'results will be saved.', required=True)

    parser.add_argument('--subjects', help='One or more subject IDs'\
                        '(space separated). If not specified, all subjects are'\
                        ' considered from  --scansdir.', \
                        default=None, required=False, nargs='+', action='append')

    parser.add_argument('-b', '--debug', help='debug mode', action='store_true')
    parser.add_argument('-p', '--processes', help='parallel processes', \
                        default=1, type=int)
    parser.add_argument('-n', '--name', help='Pipeline workflow name', 
                        default='preprocess_workflow')
    parser.add_argument('--spm',help='SPM12 executable shell script (/opt/spm12/run_spm12.sh)',\
                        required=False, default='/opt/spm12_r7771_mcrv97/run_spm12.sh')
    parser.add_argument('--mcr',help='Matlab Compiler Runtime (/opt/MCR2019b/v97)',\
                        required=False, default='/opt/MCR2019b/v97/')
 
    args = parser.parse_args()

    
    scans_dir = os.path.abspath(os.path.expandvars(args.scansdir))
    if not os.path.exists(scans_dir):
        raise IOError("Scans directory does not exist.")
        
    subject_ids = []

    if args.subjects:
        subject_ids = list(chain.from_iterable(args.subjects))
    else:
        subject_ids = glob.glob(scans_dir.rstrip('/') + '/*')
        subject_ids = [os.path.basename(s.rstrip('/')) for s in subject_ids]

    print("Creating preprocessing pipeline workflow...")
    work_dir = os.path.abspath(os.path.expandvars(args.workdir))
    out_dir =  os.path.abspath(os.path.expandvars(args.outputdir))

    if not os.path.exists(work_dir):
        os.makedirs(work_dir)

    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    
    config.update_config({
        'logging': {'log_directory': args.workdir, 'log_to_file': True},
        'execution': {'job_finished_timeout' : 15,
                      'poll_sleep_duration' : 15,
                      'hash_method' : 'content',
                      'local_hash_check' : False,
                      'stop_on_first_crash':False,
                      'crashdump_dir': args.workdir,
                      'crashfile_format': 'txt'
                       },
                       })

    #config.enable_debug_mode()
    logging.update_logging(config)

   
    preproc_wf = get_pp_wf(scans_dir=scans_dir,
                           work_dir=work_dir,
                           outputdir=out_dir,
                           subject_ids=subject_ids,
                           spm=args.spm,
                           mcr=args.mcr,
                           wfname='preproc_workflow')
     
    # Visualize workflow
    if args.debug:
        preproc_wf.write_graph(graph2use='colored', simple_form=True)

    preproc_wf.run(
                            plugin='MultiProc', 
                            plugin_args={'n_procs' : args.processes}
                    )
    
    print('Done rsfMRI preprocess pipeline!!!')

if __name__ == '__main__':
    sys.exit(main())
