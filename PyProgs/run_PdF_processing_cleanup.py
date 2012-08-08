#!/usr/bin/env python
# encoding: utf-8

import os, sys, optparse, glob

import numpy as np


from OP_waveforms import *

#get waveloc path from environment
base_path=os.getenv('WAVELOC_PATH')

p = optparse.OptionParser()
p.add_option('--datadir',action='store',help="data subdirectory")
p.add_option('--dataglob',action='store',help="data glob")
p.add_option('--starttime',action='store',help="start time for data e.g. 2010-10-14T00:00:00.0Z")
p.add_option('--endtime',action='store',help="end time for data e.g. 2010-10-14T10:00:00.0Z")
p.add_option('--short_data_length',action='store',help="isolated segments of data shorter than this (in seconds) will be removed before processing e.g. 600 ")
p.add_option('--c1',action='store',help="low frequency corner of band pass filter ")
p.add_option('--c2',action='store',help="high frequency corner of band pass filter ")
p.add_option('--kwin',action='store',help="length of kurtosis window (seconds)")

(options,arguments)=p.parse_args()


data_dir="%s/data/%s"%(base_path,options.datadir)

files=glob.glob(data_dir + os.sep + options.dataglob)


filter_c1=np.float(options.c1)
filter_c2=np.float(options.c2)
kurt_window=np.float(options.kwin)

# start and end time to process
start_time=utcdatetime.UTCDateTime(options.starttime)
end_time=utcdatetime.UTCDateTime(options.endtime)

short_data_length=np.float(options.short_data_length)

# for file in files : 
wf=Waveform()
for file in files:
  wf.read_from_file(file, starttime=start_time, endtime=end_time)
  wf.cleanup_traces(short_data_length=short_data_length)
  wf.bp_filter(filter_c1, filter_c2,rmean=True,taper=True)
  wf.write_to_file_filled("%s.filt.sac"%file,format='SAC',fill_value=0)
  wf.process_kurtosis(kurt_window,post_taper=True)
  wf.write_to_file_filled("%s.filt_kurt.sac"%file,format='SAC',fill_value=0)
  
