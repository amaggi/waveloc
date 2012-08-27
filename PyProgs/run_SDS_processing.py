#!/usr/bin/env python
# encoding: utf-8

import os, sys, optparse, glob
import numpy as np
from OP_waveforms import *
import logging

logging.basicConfig(level=logging.DEBUG, format='%(levelname)s : %(asctime)s : %(message)s')

#get waveloc path from environment
base_path=os.getenv('WAVELOC_PATH')

p = optparse.OptionParser()
p.add_option('--datadir',action='store',help="data subdirectory")
p.add_option('--net_list',action='store',help="list of network codes (e.g. \"BE,G\") ")
p.add_option('--sta_list',action='store',help="list of station names (e.g. \"STA1,STA2\") ")
p.add_option('--comp_list',action='store',help="list of component names (e.g. \"HHZ,LHZ\") ")
p.add_option('--starttime',action='store',help="start time for data e.g. 2010-10-14T00:00:00.0Z")
p.add_option('--endtime',action='store',help="end time for data e.g. 2010-10-14T10:00:00.0Z")

p.add_option('--c1',action='store',help="low frequency corner of band pass filter ")
p.add_option('--c2',action='store',help="high frequency corner of band pass filter ")
p.add_option('--kwin',action='store',help="length of kurtosis window (seconds)")

(options,arguments)=p.parse_args()


data_dir=os.path.join(base_path,'data',options.datadir)

filter_c1=np.float(options.c1)
filter_c2=np.float(options.c2)
kurt_window=np.float(options.kwin)

# start and end time to process
start_time=utcdatetime.UTCDateTime(options.starttime)
end_time=utcdatetime.UTCDateTime(options.endtime)

logging.debug(options.net_list)
logging.debug(options.sta_list)
logging.debug(options.comp_list)

for net in options.net_list.split(','):
  for sta in options.sta_list.split(','):
    for comp in options.comp_list.split(','):
      full_path=os.path.join(data_dir,net,sta,"%s.D"%comp)
      logging.debug("Full path : %s"%full_path)
      if os.path.exists(full_path):
        filt_filename=os.path.join(data_dir,"%s.%s.%s.%s.filt.mseed"%(start_time.isoformat(),net,sta,comp))
        wf=Waveform()
        try:
          wf.read_from_SDS(data_dir,net,sta,comp, starttime=start_time, endtime=end_time)

          wf.bp_filter(filter_c1, filter_c2,rmean=True,taper=True)
          #wf.resample(50.0)

          wf.write_to_file_filled(filt_filename,format='MSEED',fill_value=0)
          logging.debug(filt_filename)

          kurt_filename=os.path.join(data_dir,"%s.%s.%s.%s.filt_kurt.mseed"%(start_time.isoformat(),net,sta,comp))
          wf.process_kurtosis(kurt_window,post_taper=True)
          wf.write_to_file_filled(kurt_filename,format='MSEED',fill_value=0)
          logging.debug(kurt_filename)
  
        except UserWarning:
          logging.info('No data within time limits for %s %s %s'%(net,sta,comp))
