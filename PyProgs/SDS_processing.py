#!/usr/bin/env python
# encoding: utf-8

import os, sys, optparse, glob
import numpy as np
from OP_waveforms import *
import logging


def do_SDS_processing_setup_and_run(datadir='TEST',net_list="",sta_list="",comp_list="",starttime=None,endtime=None,resample=False,fs=None,c1=None,c2=None,kwin=None,krec=False,kderiv=False):

  #get waveloc path from environment
  base_path=os.getenv('WAVELOC_PATH')

  data_dir=os.path.join(base_path,'data',datadir)

  filter_c1=c1
  filter_c2=c2
  kurt_window=kwin

  # start and end time to process
  start_time=utcdatetime.UTCDateTime(starttime)
  end_time=utcdatetime.UTCDateTime(endtime)

  logging.debug('Network list = %s'%net_list)
  logging.debug('Station list = %s'%sta_list)
  logging.debug('Component list = %s'%comp_list)

  net_list=net_list.split(',')
  sta_list=sta_list.split(',')
  comp_list=comp_list.split(',')

  # loop over data
  for net in net_list:
    for sta in sta_list:
      for comp in comp_list:
        full_path=os.path.join(data_dir,net,sta,"%s.D"%comp)
        logging.debug("Full path : %s"%full_path)
        if os.path.exists(full_path):
          filt_filename=os.path.join(data_dir,"%s.%s.%s.%s.filt.mseed"%(start_time.isoformat(),net,sta,comp))
          logging.info("Processing to create %s" % (filt_filename))

          wf=Waveform()
          try:
            wf.read_from_SDS(data_dir,net,sta,comp, starttime=start_time, endtime=end_time)
            wf.bp_filter(filter_c1, filter_c2,rmean=True,taper=True)
            if resample:
              wf.resample(np.float(fs))
            wf.write_to_file_filled(filt_filename,format='MSEED',fill_value=0)

            kurt_filename=os.path.join(data_dir,"%s.%s.%s.%s.filt_kurt.mseed"%(start_time.isoformat(),net,sta,comp))
            logging.info("Processing to create %s " % (kurt_filename))
            wf.process_kurtosis(kurt_window,recursive=krec,pre_taper=True, post_taper=True)
            wf.write_to_file_filled(kurt_filename,format='MSEED',fill_value=0)

            if kderiv:
              kurt_grad_filename=os.path.join(data_dir,"%s.%s.%s.%s.filt_kurt_grad.mseed"%(start_time.isoformat(),net,sta,comp))
              logging.info("Processing to create %s " % (kurt_grad_filename))
              wf.take_positive_derivative(pre_taper=True,post_taper=True)
              wf.write_to_file_filled(kurt_grad_filename,format='MSEED',fill_value=0)
  
          except UserWarning:
            logging.info('No data within time limits for %s %s %s'%(net,sta,comp))



if __name__ == '__main__' :

  logging.basicConfig(level=logging.DEBUG, format='%(levelname)s : %(asctime)s : %(message)s')

  p = optparse.OptionParser()
  p.add_option('--datadir',action='store',help="data subdirectory")
  p.add_option('--net_list',action='store',help="list of network codes (e.g. \"BE,G\") ")
  p.add_option('--sta_list',action='store',help="list of station names (e.g. \"STA1,STA2\") ")
  p.add_option('--comp_list',action='store',help="list of component names (e.g. \"HHZ,LHZ\") ")
  p.add_option('--starttime',action='store',help="start time for data e.g. 2010-10-14T00:00:00.0Z")
  p.add_option('--endtime',action='store',help="end time for data e.g. 2010-10-14T10:00:00.0Z")
  p.add_option('--resample',action='store_true',default=False, help="resample data")
  p.add_option('--fs',action='store',default=False, help="resample frequency")

  p.add_option('--c1',action='store',type='float',help="low frequency corner of band pass filter ")
  p.add_option('--c2',action='store',type='float',help="high frequency corner of band pass filter ")
  p.add_option('--kwin',action='store',type='float',help="length of kurtosis window (seconds)")
  p.add_option('--krec',action='store_true',default=False, help="use recursive kurtosis calculation (faster but less precise)")
  p.add_option('--kderiv',action='store_true',default=False, help="use derivative of kurtosis")

  (options,arguments)=p.parse_args()

  do_SDS_processing_setup_and_run(
     datadir=options.datadir,
     net_list=options.net_list,
     sta_list=options.sta_list,
     comp_list=options.comp_list,
     starttime=options.starttime,
     endtime=options.endtime,
     resample=options.resample,
     fs=options.fs,
     c1=options.c1,
     c2=options.c2,
     kwin=options.kwin,
     krec=options.krec,
     kderiv=options.kderiv)

