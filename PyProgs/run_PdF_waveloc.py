#!/usr/bin/env python
# encoding: utf-8

import os, sys, optparse

import numpy as np


from OP_waveforms import *

from grids_paths import *
from time import time, sleep
from sub_PdF_waveloc import do_migration_loop_continuous

import cProfile

import logging
logging.basicConfig(level=logging.DEBUG, format='%(levelname)s : %(asctime)s : %(message)s')


def main(options):
  # set variables from command line options
  # grid
  grid_filename_base=os.path.join(base_path,'lib',options.time_grid)
  search_grid_filename=os.path.join(base_path,'lib',options.search_grid)
  # stations
  stations_filename=os.path.join(base_path,'lib',options.stations)
  # output directory
  if options.outdir==None: raise UserWarning("No output subdirectory supplied")
  output_dir=os.path.join(base_path,'out',options.outdir)
  if not os.path.exists(output_dir):
    os.makedirs(output_dir)
  # data
  if options.datadir==None: raise UserWarning("No data subdirectory supplied")
  if options.dataglob==None: raise UserWarning("No data glob supplied")
  data_dir=os.path.join(base_path,'data',options.datadir)
  data_glob=options.dataglob


  if options.verbose: 
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

  if options.verbose:
    print "Reading station file"

  if options.time:
    t_ref=time()

  sta=StationList()
  sta.read_from_file(stations_filename)

  if options.time:
    t=time()-t_ref
    print "Time for reading %d stations from file : %.4f s\n" % (sta.nsta,t)

  datafile_list=glob.glob(os.path.join(data_dir,data_glob))

  cha=ChannelList()
  cha.populate_from_station_list_and_data_files(sta,datafile_list)


  ######### INTERPOLATE TRAVEL TIMES #############

  # The time grid will contain as array values just the travel-times needed 
  # (interpolated from the full NLL files) so we can free up the memory as soon as possible

  if options.verbose:
    print "Extracting useful travel-times"

  if options.time:
    t_ref=time()  

  time_grid=QDTimeGrid()
  time_grid.read_NLL_hdr_file(search_grid_filename)
  time_grid.populate_from_time_grids(grid_filename_base,output_dir,cha)

  if options.time:
    t=time()-t_ref
    print "Time for extracting and saving %dx%dx%dx%d travel-times : %.2f s\n" % (time_grid.nx,time_grid.ny,time_grid.nz,cha.ncha,t)



  # READ DATA

  # start loop over time
  initial_start_time=utcdatetime.UTCDateTime(options.starttime)
  initial_end_time=initial_start_time+np.float(options.data_length)

  final_end_time=utcdatetime.UTCDateTime(options.endtime)

  time_shift_secs=np.float(options.data_length)-np.float(options.data_overlap)

  start_time=initial_start_time
  end_time=initial_end_time

  while (start_time < final_end_time):

    # DO MIGRATION !!
    do_migration_loop_continuous(start_time, end_time, data_dir, output_dir, data_glob, search_grid_filename, time_grid, options.verbose, options.time)


    # Reset the start and end times to loop again
    start_time=start_time+time_shift_secs
    end_time=end_time+time_shift_secs



#get waveloc path from environment
base_path=os.getenv('WAVELOC_PATH')

# lists of pre-calculated grids 
time_grids=['Slow_len.100m.P','belgium.P','ijen.P']
search_grids=['grid.500m.search.hdr','grid.Taisne.search.hdr','grid.belgium.search.hdr','grid.ijen.search.hdr']

# setup arguments for command line parsing
p = optparse.OptionParser()
p.add_option('--time','-t',action='store_true',help='print timing information to stout')
p.add_option('--verbose','-v',action='store_true',help='print debugging information to stout')
p.add_option('--time_grid',action='store',type='choice',choices=time_grids,help="time grid %s"%(time_grids))
p.add_option('--search_grid',action='store',type='choice',choices=search_grids,help="search grid %s"%(search_grids))
p.add_option('--stations','-s',action='store',default='channels_HHZ.dat',help='station list (found in $WAVELOC_PATH/lib)')
p.add_option('--outdir','-o',action='store',help="output subdirectory")
p.add_option('--datadir',action='store',help="data subdirectory")
p.add_option('--dataglob',action='store',help="data glob")
p.add_option('--starttime',action='store',help="start time for data e.g. 2010-10-14T00:00:00.0Z")
p.add_option('--endtime',action='store',help="end time for data e.g. 2010-10-14T10:00:00.0Z")
p.add_option('--data_length',action='store',help="length in seconds for data segments to analyse (e.g. 630)")
p.add_option('--data_overlap',action='store',help="length in seconds for overlapping data segments (e.g. 30)")
p.add_option('--profile',action='store_true',help="Do profiling")



# parse command line
(options,arguments)=p.parse_args()

if options.profile:
  cProfile.run('main(options)','run_waveloc.prof')
else:
  main(options)
