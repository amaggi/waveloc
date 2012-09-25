#!/usr/bin/env python
# encoding: utf-8

import os, sys, optparse

import numpy as np


from OP_waveforms import *

from grids_paths import *
from time import time, sleep
from sub_PdF_waveloc import do_migration_loop_continuous
from NllGridLib import read_stations_file
from hdf5_grids import get_interpolated_time_grids
import logging

def do_migration_setup_and_run(opdict):

  base_path=opdict['base_path']
  verbose=opdict['verbose']
  runtime=opdict['time']
  

  # set variables from command line options
  # grid
  grid_filename_base=os.path.join(base_path,'lib',opdict['time_grid'])
  search_grid_filename=os.path.join(base_path,'lib',opdict['search_grid'])

  # stations
  stations_filename=os.path.join(base_path,'lib',opdict['stations'])

  # output directory
  output_dir=os.path.join(base_path,'out',opdict['outdir'])
  stack_dir=os.path.join(output_dir,'stack')

  # data
  data_dir=os.path.join(base_path,'data',opdict['datadir'])
  data_glob=opdict['gradglob']


  if verbose: 
    print ""
    print "Input parameters:"
    print "-----------------"
    print "Grid        = %s"%(grid_filename_base)
    print "Stations    = %s"%(stations_filename)
    print "Data        = %s"%(os.path.join(data_dir,data_glob))
    print ""
    print "Output parameters:"
    print "-----------------"
    print "Out dir     = %s"%(output_dir)
  

  #raise UserWarning ('Stop here')
  #######################################################################
  #                       START PROCESSING 
  #######################################################################

  if verbose: 
    print ""
    print "----------------"
    print "START PROCESSING"
    print "----------------"
    print ""

  # Create Obspy streams for output

  #  ***** reading station file ******

  if verbose:
    logging.info("Reading station file")

  if runtime:
    t_ref=time()

  stations=read_stations_file(stations_filename)

  if runtime:
    t=time()-t_ref
    logging.info("Time for reading %d stations from file : %.4f s\n" % (len(stations.keys()),t))

  datafile_list=glob.glob(os.path.join(data_dir,data_glob))


  ########### DEAL WITH START AND END TILES ############
  starttime=opdict['starttime']
  endtime=opdict['endtime']
  data_length=opdict['data_length']
  data_overlap=opdict['data_overlap']

  initial_start_time=utcdatetime.UTCDateTime(starttime)
  initial_end_time=initial_start_time+data_length

  final_end_time=utcdatetime.UTCDateTime(endtime)

  time_shift_secs=data_length-data_overlap



  ######### INTERPOLATE TRAVEL TIMES #############

  # The time grid will contain as array values just the travel-times needed 
  # (interpolated from the full NLL files) so we can free up the memory as soon as possible

  if verbose:
    logging.info("Extracting useful travel-times")

  if runtime:
    t_ref=time()  

  time_grids=get_interpolated_time_grids(opdict)

  if runtime:
    t=time()-t_ref
    logging.info("Time for retreiving %d travel time grids : %.2f s\n" % (len(time_grids.keys()),t))


  # READ DATA

  # start loop over time
  start_time=initial_start_time
  end_time=initial_end_time

  if runtime:
    t_ref=time()  

  while (start_time < final_end_time):

    do_migration_loop_continuous(start_time, end_time, data_dir, output_dir, data_glob, search_grid_filename, time_grids, verbose, runtime)

    # Reset the start and end times to loop again
    start_time=start_time+time_shift_secs
    end_time=end_time+time_shift_secs

  if runtime:
    t=time()-t_ref
    logging.info("Time for migrating : %.2f s\n" % (t))




if __name__ == '__main__':

  from options import WavelocOptions
  logging.basicConfig(level=logging.INFO, format='%(levelname)s : %(asctime)s : %(message)s')

  wo = WavelocOptions()
  args=wo.p.parse_args()

  wo.set_all_arguments(args)
  wo.verify_migration_options()

  do_migration_setup_and_run(wo.opdict)


