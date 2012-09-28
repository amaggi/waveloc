#!/usr/bin/env python
# encoding: utf-8

import os, sys, optparse, glob
import numpy as np
import logging


logging.basicConfig(level=logging.DEBUG, format='%(levelname)s : %(asctime)s : %(message)s')

data_dir=os.path.abspath('../../data/LAquila/Data')
nll_file='nll_ijen.in'

nl=open(nll_file)
lines=nl.readlines()
nl.close()

stalines=[line for line in lines if line.split().[0]=]


all_files=glob.glob(os.path.join(data_dir,'2009*','*.sac'))
all_short_files=[os.path.basename(filename) for filename in all_files]
sta_list=set([filename.split('.')[0]) for filename in all_short_files])

print sta


