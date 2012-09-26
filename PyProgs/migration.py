#!/usr/bin/env python
# encoding: utf-8

import os, sys, optparse, h5py

import numpy as np

from OP_waveforms import *

from grids_paths import *
from time import time, sleep
from sub_PdF_waveloc import do_migration_loop_continuous
from NllGridLib import read_stations_file,read_hdr_file
from hdf5_grids import get_interpolated_time_grids
import logging

def do_migration_setup_and_run(opdict):

  base_path=opdict['base_path']
  verbose=opdict['verbose']
  runtime=opdict['time']

  # stations
  stations_filename=os.path.join(base_path,'lib',opdict['stations'])
  stations=read_stations_file(stations_filename)

  # output directory
  output_dir=os.path.join(base_path,'out',opdict['outdir'])
  stack_dir=os.path.join(output_dir,'stack')

  # data
  data_dir=os.path.join(base_path,'data',opdict['datadir'])
  data_glob=opdict['gradglob']
  data_files=glob.glob(os.path.join(data_dir,data_glob))
  if len(data_files)==0: 
    logging.error('No data files found for %s and %s'%(data_dir,data_glob))
    raise UserWarning

  # grids
  grid_filename_base=os.path.join(base_path,'lib',opdict['time_grid'])
  search_grid_filename=os.path.join(base_path,'lib',opdict['search_grid'])
  grid_info=read_hdr_file(search_grid_filename)
  time_grids=get_interpolated_time_grids(opdict)

  #start and end times
  starttime=opdict['starttime']
  endtime=opdict['endtime']
  data_length=opdict['data_length']
  data_overlap=opdict['data_overlap']

  initial_start_time=utcdatetime.UTCDateTime(starttime)
  initial_end_time=initial_start_time+data_length

  final_end_time=utcdatetime.UTCDateTime(endtime)

  time_shift_secs=data_length-data_overlap


  ######### FOR EACH TIME SPAN - DO MIGRATION #############

  # start loop over time
  start_time=initial_start_time
  end_time=initial_end_time

  if runtime:
    t_ref=time()  

  while (start_time < final_end_time):

    # read data
    logging.info("Reading data  : %s - %s."%(start_time.isoformat(), end_time.isoformat()))
    data,delta=read_data_compatible_with_time_dict(data_files,time_grids,start_time,end_time)

    # do migration if have enough data (3 is bare minimum)
    if len(data.keys())>=3:
      logging.info("Migrating data : %s - %s."%(start_time.isoformat(), end_time.isoformat()))
      do_migration_loop_continuous(opdict, data, delta, start_time, end_time, grid_info, time_grids)
    elif len(data.keys())==0:
      logging.warn('No data found between %s and %s.'%(start_time.isoformat(),end_time.isoformat()))
    else:
      logging.warn('Insufficient data found between %s and %s.'%(start_time.isoformat(),end_time.isoformat()))
      
    # Reset the start and end times to loop again
    start_time=start_time+time_shift_secs
    end_time=end_time+time_shift_secs

  if runtime:
    t=time()-t_ref
    logging.info("Time for migrating all time slices : %.2f s\n" % (t))




if __name__ == '__main__':

  from options import WavelocOptions
  logging.basicConfig(level=logging.INFO, format='%(levelname)s : %(asctime)s : %(message)s')

  wo = WavelocOptions()
  args=wo.p.parse_args()

  wo.set_all_arguments(args)
  wo.verify_migration_options()

  do_migration_setup_and_run(wo.opdict)


