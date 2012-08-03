#!/usr/bin/env python
# encoding: utf-8

import os, sys, optparse, glob

import numpy as np


from OP_waveforms import *

#get waveloc path from environment
base_path=os.getenv('WAVELOC_PATH_PDF')

p = optparse.OptionParser()
p.add_option('--datadir',action='store',help="data subdirectory")

(options,arguments)=p.parse_args()


data_dir="%s/data/%s"%(base_path,options.datadir)

# list of files to process
#data_glob="FLR*1day"
#data_glob="RVL*1day"
#data_glob="[FHS]*1day"
data_glob="*1day"
files=glob.glob(data_dir + os.sep + data_glob)

#parametres de filtrage et kurtosis 
filter_c1=4.0
filter_c2=10.0
kurt_window=3.0

# start and end time to process
start_time=utcdatetime.UTCDateTime("2010-10-14T00:00:00.0Z")
#start_time=utcdatetime.UTCDateTime("2010-10-14T03:50:00.0Z")
#end_time=utcdatetime.UTCDateTime("2010-10-14T05:00:00.0Z")
#start_time=utcdatetime.UTCDateTime("2010-10-14T09:00:00.0Z")
end_time=utcdatetime.UTCDateTime("2010-10-14T16:00:00.0Z")
#end_time=utcdatetime.UTCDateTime("2010-10-14T04:00:00.0Z")

# for file in files : 
wf=Waveform()
for file in files:
  wf.read_from_file(file, starttime=start_time, endtime=end_time)
  wf.bp_filter(filter_c1, filter_c2)
  wf.write_to_file_filled("%s.filt.sac"%file,format='SAC')
  wf.process_kurtosis(kurt_window)
  wf.write_to_file_filled("%s.filt_kurt.sac"%file,format='SAC')
  
