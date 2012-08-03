#!/usr/bin/env python
# encoding: utf-8

import os, sys, optparse, glob
import numpy as np
import logging
from obspy.core import read, utcdatetime, trace
from OP_waveforms import Waveform
from kurtogram import Fast_Kurtogram

logging.basicConfig(level=logging.INFO, format='%(levelname)s : %(asctime)s : %(message)s')

p=optparse.OptionParser()
p.add_option('--data_file', '-d', action='store',help="data filename")
p.add_option('--output_file','-o', action='store',help="output filename")
p.add_option('--starttime',action='store',help="start time for data e.g. 2010-10-14T00:00:00.0Z")
p.add_option('--endtime',action='store',help="end time for data e.g. 2010-10-14T10:00:00.0Z")
p.add_option('--verbose','-v',action='store_true',help='print debugging information to stout')
(options,arguments)=p.parse_args()

tdeb=utcdatetime.UTCDateTime(options.starttime)
tfin=utcdatetime.UTCDateTime(options.endtime)

if options.verbose:
  logging.info('\n\
  Input filename  = %s\n\
  Output basename = %s\n\
  Start time = %s\n\
  End time   = %s\n' % (options.data_file, options.output_file, options.starttime, options.endtime))

# read waveform between time limits
wf=Waveform()
wf.read_from_file(options.data_file,starttime=tdeb,endtime=tfin)
dt=wf.delta
x=wf.values
print(wf.stream)

# set up parameters for kurtogram analysis
N=len(x)
N2=np.log2(N)-7
nlevel=int(np.fix(N2))
c,flower,fupper = Fast_Kurtogram(x, nlevel,options.verbose,Fs=1/dt,opt2=1)

logging.info("Frequency band for best kurtogram : %.2f Hz - %.2f Hz"%(flower,fupper))

filt_name="filt_%s"%options.output_file
kurt_name="kurt_%s"%options.output_file
wf.bp_filter(flower,fupper,rmean=True,taper=True)
wf.write_to_file_filled(filt_name,format='MSEED')
wf.process_kurtosis(100*dt,recursive=True,post_taper=True)
wf.write_to_file_filled(kurt_name,format='MSEED')

