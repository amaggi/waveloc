#!/usr/bin/env python
# encoding: utf-8

import os, sys, optparse

import numpy as np


from OP_waveforms import *

from grids_paths import *
from time import time, sleep
from sub_PdF_waveloc import do_migration_loop_continuous
import Queue, threading
import logging

def do_migration_setup_and_run(base_path="",runtime=False, verbose=False, twoD=False, time_grid="", search_grid="", stations="", outdir="", datadir="", dataglob="", starttime="", endtime="", data_length=None, data_overlap=None, load_ttimes_buf=False):

  # set variables from command line options
  # grid
  logging.debug(base_path)
  logging.debug(time_grid)
  grid_filename_base=os.path.join(base_path,'aux',time_grid)
  search_grid_filename=os.path.join(base_path,'aux',search_grid)

  # stations
  stations_filename=os.path.join(base_path,'aux',stations)

  # output directory
  output_dir=os.path.join(base_path,'out',outdir)
  stack_dir=os.path.join(outdir,'stack')
  if not os.path.exists(stack_dir):
    os.makedirs(stack_dir)

  # data
  data_dir=os.path.join(base_path,'data',datadir)
  data_glob=dataglob


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

  sta=StationList()
  sta.read_from_file(stations_filename)

  if runtime:
    t=time()-t_ref
    logging.info("Time for reading %d stations from file : %.4f s\n" % (sta.nsta,t))

  datafile_list=glob.glob(os.path.join(data_dir,data_glob))

  cha=ChannelList()
  cha.populate_from_station_list_and_data_files(sta,datafile_list)


  ########### DEAL WITH START AND END TILES ############

  initial_start_time=utcdatetime.UTCDateTime(starttime)
  initial_end_time=initial_start_time+np.float(data_length)

  final_end_time=utcdatetime.UTCDateTime(endtime)

  time_shift_secs=data_length-data_overlap



  ######### INTERPOLATE TRAVEL TIMES #############

  # The time grid will contain as array values just the travel-times needed 
  # (interpolated from the full NLL files) so we can free up the memory as soon as possible

  if verbose:
    logging.info("Extracting useful travel-times")

  if runtime:
    t_ref=time()  

  time_grid=QDTimeGrid()
  time_grid.read_NLL_hdr_file(search_grid_filename)
  if twoD:
    time_grid.populate_from_2D_time_grids(grid_filename_base,cha)
  else:
    time_grid.populate_from_time_grids(grid_filename_base,cha,load_ttimes_buf)

  if runtime:
    t=time()-t_ref
    n_times=time_grid.nx*time_grid.ny*time_grid.nz+cha.ncha
    logging.info("Time for extracting and saving %dx%dx%dx%d=%d travel-times : %.2f s\n" % (time_grid.nx,time_grid.ny,time_grid.nz,cha.ncha,n_times,t))


  # READ DATA

  # start loop over time
  start_time=initial_start_time
  end_time=initial_end_time

  if runtime:
    t_ref=time()  

  while (start_time < final_end_time):

    do_migration_loop_continuous(start_time, end_time, data_dir, output_dir, data_glob, search_grid_filename, time_grid, verbose, runtime)

    # Reset the start and end times to loop again
    start_time=start_time+time_shift_secs
    end_time=end_time+time_shift_secs

  if runtime:
    t=time()-t_ref
    n_times=time_grid.nx*time_grid.ny*time_grid.nz+cha.ncha
    logging.info("Time for migrating : %.2f s\n" % (t))




if __name__ == '__main__':

  logging.basicConfig(level=logging.INFO, format='%(levelname)s : %(asctime)s : %(message)s')


  #get waveloc path from environment
  base_path=os.getenv('WAVELOC_PATH')

  # setup arguments for command line parsing
  p = optparse.OptionParser()
  p.add_option('--time','-t',action='store_true',help='print timing information to stout')
  p.add_option('--verbose','-v',action='store_true',default=False,help='print debugging information to stout')
  p.add_option('--2D',action='store_true',default=False,dest='twoD',help='use 2D time grids')
  p.add_option('--time_grid',action='store',help="time grid basename e.g. belgium.P (found in $WAVELOC_PATH/aux)")
  p.add_option('--search_grid',action='store',help="search grid e.g. grid.500m.search.hdr (found in $WAVELOC_PATH/aux)")
  p.add_option('--stations','-s',action='store',help='station list (found in $WAVELOC_PATH/aux)')
  p.add_option('--outdir','-o',action='store',help="output subdirectory")
  p.add_option('--datadir',action='store',help="data subdirectory")
  p.add_option('--dataglob',action='store',help="data glob")
  p.add_option('--starttime',action='store',help="start time for data e.g. 2010-10-14T00:00:00.0Z")
  p.add_option('--endtime',action='store',help="end time for data e.g. 2010-10-14T10:00:00.0Z")
  p.add_option('--data_length',action='store',type='float',help="length in seconds for data segments to analyse (e.g. 630)")
  p.add_option('--data_overlap',action='store',type='float',help="length in seconds for overlapping data segments (e.g. 30)")
  p.add_option('--load_ttimes_buf',action='store_true',default=False,dest='load_buf',help='load pre-calculated travel-times for the search grid from file')

  # parse command line
  (options,arguments)=p.parse_args()

  if options.outdir==None: raise UserWarning("No output subdirectory supplied")
  if options.datadir==None: raise UserWarning("No data subdirectory supplied")
  if options.dataglob==None: raise UserWarning("No data glob supplied")

  base_path=os.getenv('WAVELOC_PATH')

  do_migration_setup_and_run(base_path=base_path,
    runtime=options.time,
    verbose=options.verbose,
    twoD=options.twoD,
    time_grid=options.time_grid,
    search_grid=options.search_grid,
    stations=options.stations,
    outdir=options.outdir,
    datadir=options.datadir,
    dataglob=options.dataglob,
    starttime=options.starttime,
    endtime=options.endtime,
    data_length=options.data_length,
    data_overlap=options.data_overlap,
    load_ttimes_buf=options.load_buf)


