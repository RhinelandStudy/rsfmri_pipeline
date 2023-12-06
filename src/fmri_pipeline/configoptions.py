#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 13 15:18:10 2017

@author: shahidm
"""

from os.path import realpath, join, abspath, dirname


# defaults
SCRIPT_PATH = dirname(realpath(__file__))

BOLD_MASK = abspath(join(SCRIPT_PATH, 'FC_templates','BOLD_mask.nii'))

