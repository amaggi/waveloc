#!/usr/bin/env python
# encoding: utf-8

import os, sys, optparse, glob
import numpy as np
from OP_waveforms import *
import logging


def do_SDS_processing_setup_and_run(opdict):

  #get waveloc path from environment
  base_path=os.getenv('WAVELOC_PATH')

  data_dir=os.path.join(base_path,'data',opdict['datadir'])

  filter_c1=opdict['c1']
  filter_c2=opdict['c2']
  kurt_window=opdict['kwin']

  # start and end time to process
  start_time=utcdatetime.UTCDateTime(opdict['starttime'])
  end_time=utcdatetime.UTCDateTime(opdict['endtime'])

  logging.debug('Network list = %s'%opdict['net_list'])
  logging.debug('Station list = %s'%opdict['sta_list'])
  logging.debug('Component list = %s'%opdict['comp_list'])

  net_list=opdict['net_list'].split(',')
  sta_list=opdict['sta_list'].split(',')
  comp_list=opdict['comp_list'].split(',')

  # loop over data
  for net in net_list:
    for sta in sta_list:
      for comp in comp_list:
        full_path=os.path.join(data_dir,net,sta,"%s.D"%comp)
        logging.debug("Full path : %s"%full_path)
        if os.path.exists(full_path):

          filt_filename=os.path.join(data_dir,"%s.%s.%s.%s.filt.mseed"%(start_time.isoformat(),net,sta,comp))
          logging.debug("Processing to create %s" % (filt_filename))
          wf=Waveform()
          try:
            wf.read_from_SDS(data_dir,net,sta,comp, starttime=start_time, endtime=end_time)
            wf.bp_filter(filter_c1, filter_c2,rmean=True,taper=True)
            if opdict['resample']:
              wf.resample(opdict['fs'])
            wf.write_to_file_filled(filt_filename,format='MSEED',fill_value=0)

            kurt_filename=os.path.join(data_dir,"%s.%s.%s.%s.filt_kurt.mseed"%(start_time.isoformat(),net,sta,comp))
            logging.debug("Processing to create %s" % (kurt_filename))
            wf.process_kurtosis(kurt_window,recursive=opdict['krec'],pre_taper=True, post_taper=True)
            wf.write_to_file_filled(kurt_filename,format='MSEED',fill_value=0)

            if opdict['kderiv']:
              kurt_grad_filename=os.path.join(data_dir,"%s.%s.%s.%s.filt_kurt_grad.mseed"%(start_time.isoformat(),net,sta,comp))
              logging.debug("Processing to create %s" % (kurt_grad_filename))
              wf.take_positive_derivative(pre_taper=True,post_taper=True)
              wf.write_to_file_filled(kurt_grad_filename,format='MSEED',fill_value=0)
  
          except UserWarning:
            logging.info('No data within time limits for %s %s %s'%(net,sta,comp))



if __name__ == '__main__' :

  from options import WavelocOptions

  logging.basicConfig(level=logging.INFO, format='%(levelname)s : %(asctime)s : %(message)s')

  wo = WavelocOptions()
  args=wo.p.parse_args()

  wo.set_all_arguments(args)
  wo.verify_SDS_processing_options()

  do_SDS_processing_setup_and_run(wo.opdict)
