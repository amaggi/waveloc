#!/usr/bin/env python
# encoding: utf-8

import os, sys

import numpy as np
import numexpr as ne
import carray as ca

from OP_waveforms import *

from grids_paths import *
from time import time, sleep
from scipy import weave
from scipy.weave import converters
import logging

def do_innermost_migration_loop(start_time, end_time, data, time_grid, delta, search_grid_filename,options_verbose=False, options_time=False):

  if options_time:
    t_ref=time()  

  time_dict=time_grid.buf[0]
  time_ids=time_dict.keys()

  data_ids=data.keys()
      
  # create quick and dirty integer versions of the kurtosed data for stacking purposes
  # integer data = data  truncated to integer value and stored in 16 bits

  #logging.info("\nCreating 16 bit integer version of the data for faster stacking...")

  integer_data={}
  for key,wf in data.iteritems():
    integer_data[key]=wf.values

  if options_time:
    t=time()-t_ref
    logging.info("Time for reading %d data streams with %d points : %.4f s\n" % (len(data.keys()),wf.npts,t))

 ######### DO THE MIGRATION #############

  logging.info("Stacking shifted time series... ")

  if options_time:
    t_ref=time()  

  (n_buf, norm_stack_len, stack_shift_time, stack_grid) = migrate_4D_stack(integer_data, delta, search_grid_filename, time_grid)
  stack_start_time=start_time-stack_shift_time

  logging.debug("Stack geographical dimension = %d"%n_buf)
  logging.debug("Stack time extent = %d points = %.2f s"%(norm_stack_len, norm_stack_len*delta))
  logging.debug("Start time of stack (wrt start time of data)= %.2f s"%(-stack_shift_time))
  logging.debug("Start time of stack %s"%(stack_start_time.isoformat()))

  if options_time:
    t=time()-t_ref
    logging.info("Time for stacking and saving %d stacks, each of extent %d points : %.2f s\n" % (n_buf,norm_stack_len,t))
 
  return n_buf, norm_stack_len, stack_shift_time, stack_start_time, stack_grid

def do_inner_migration_loop(start_time, end_time, data, time_grid, delta, search_grid_filename, options_verbose=False, options_time=False):

 
  (n_buf, norm_stack_len, stack_shift_time, stack_start_time, stack_grid)= do_innermost_migration_loop(start_time, end_time, data, time_grid, delta, search_grid_filename,options_verbose, options_time)

  ###### Extract maximum of stack #######

  logging.info("Extracting maximum of stack")

  if options_time:
    t_ref=time()  

  
  # magic matrix manipulations (see relevant inotebook)
  #ib_max=[stack_grid[:,:,:,it].argmax for it in range(norm_stack_len)]
  
  max_val=stack_grid.max(0).max(0).max(0)
  max_x=stack_grid.max(2).max(1).argmax(0)
  max_y=stack_grid.max(2).max(0).argmax(0)
  max_z=stack_grid.max(1).max(0).argmax(0)

  dx=time_grid.dx
  dy=time_grid.dy
  dz=time_grid.dz
  x_orig=time_grid.x_orig
  y_orig=time_grid.y_orig
  z_orig=time_grid.z_orig

  #go from indexes to coordinates
  max_x=ne.evaluate('max_x*dx+x_orig')
  max_y=ne.evaluate('max_y*dy+y_orig')
  max_z=ne.evaluate('max_z*dz+z_orig')

  # iterate over stack
#  for itime in range(norm_stack_len):
#    time_slice=stack_grid[:,:,:,itime].flatten()
#    ib_max=np.argmax(time_slice)
#    max_val[itime]=time_slice[ib_max]
#    ix,iy,iz=time_grid.get_ix_iy_iz(ib_max)
#    max_x[itime]=ix*time_grid.dx+time_grid.x_orig
#    max_y[itime]=iy*time_grid.dy+time_grid.y_orig
#    max_z[itime]=iz*time_grid.dz+time_grid.z_orig

  if options_time:
    t=time()-t_ref
    logging.info("Time for extracting maxima from %d slices : %.2f s\n" % (norm_stack_len,t))


  return(max_val,max_x,max_y,max_z,stack_start_time,norm_stack_len,stack_grid)
 

def do_write_stack_files(max_val,max_x,max_y,max_z,delta,stack_start_time,norm_stack_len,output_dir,stack_basename):
 
  ######## Transform results into waveform ###########

  logging.info("Writing stack files as seismograms...")

  # create headers
  stats_stack={'network': 'UV', 'station' : 'STACK', 'location': '', 'channel': 'STK', 'npts' : norm_stack_len, 'sampling_rate': 1/delta, 'mseed': {'dataquality': 'R'}}
  stats_x={'network': 'UV', 'station' : 'XYZ', 'location': '', 'channel': 'X', 'npts' : norm_stack_len, 'sampling_rate': 1/delta, 'mseed': {'dataquality': 'R'}}
  stats_y={'network': 'UV', 'station' : 'XYZ', 'location': '', 'channel': 'Y', 'npts' : norm_stack_len, 'sampling_rate': 1/delta, 'mseed': {'dataquality': 'R'}}
  stats_z={'network': 'UV', 'station' : 'XYZ', 'location': '', 'channel': 'Z', 'npts' : norm_stack_len, 'sampling_rate': 1/delta, 'mseed': {'dataquality': 'R'}}

  # set time
  stats_stack['starttime']=stack_start_time
  stats_x['starttime']=stack_start_time
  stats_y['starttime']=stack_start_time
  stats_z['starttime']=stack_start_time

  # create traces
  st_stack=Stream([Trace(data=max_val, header=stats_stack)])
  st_x=Stream([Trace(data=max_x, header=stats_x)])
  st_y=Stream([Trace(data=max_y, header=stats_y)])
  st_z=Stream([Trace(data=max_z, header=stats_z)])

  stack_file=os.path.join(output_dir,'stack',"%s_max_%s.mseed"%(stack_basename,stack_start_time.isoformat()))
  x_file=os.path.join(output_dir,'stack',"%s_x_%s.mseed"%(stack_basename,stack_start_time.isoformat()))
  y_file=os.path.join(output_dir,'stack',"%s_y_%s.mseed"%(stack_basename,stack_start_time.isoformat()))
  z_file=os.path.join(output_dir,'stack',"%s_z_%s.mseed"%(stack_basename,stack_start_time.isoformat()))

  st_stack.write(stack_file, format='MSEED')
  st_x.write(x_file, format='MSEED')
  st_y.write(y_file, format='MSEED')
  st_z.write(z_file, format='MSEED')

def do_write_grids(stack_grid,time_step_sec,delta,norm_stack_len,stack_start_time,output_dir,grid_basename):

  # write grids every time_step_sec seconds
  time_step=int(floor((time_step_sec/delta)))
  itimes=numpy.arange(0,norm_stack_len,time_step)

  for itime in itimes:
    slicetime=stack_start_time+itime*delta
    timestamp=slicetime.isoformat()
    grid_file=os.path.join(output_dir,'grid',"%s_%s.dat"%(grid_basename,timestamp))
    stack_grid.write_grid_timeslice(itime=itime,filename=grid_file)


def do_migration_loop_continuous(start_time, end_time, data_dir, output_dir, data_glob, search_grid_filename, time_grid, options_verbose, options_time):


  logging.info("Processing time slice %s"%start_time.isoformat())

  time_dict=time_grid.buf[0]

  # read data into a dictionary

  logging.info("Reading processed data into dictionary")

  if options_time:
    t_ref=time()

  # Read in the data
  files=glob.glob(os.path.join(data_dir,data_glob))
  data=read_data_compatible_with_time_dict(files,time_dict,start_time,end_time)

  # Set the global variable delta (dt for all the seismograms)
  try:
    delta=data.values()[0].delta
  except IndexError:
    raise UserWarning("File list empty - check --dataglob option")

  # DO MIGRATION
  (max_val,max_x,max_y,max_z,stack_start_time,norm_stack_len,stack_grid)=do_inner_migration_loop(start_time, end_time, data, time_grid, delta, search_grid_filename,options_verbose, options_time)

  # WRITE STACK FILES 
  do_write_stack_files(max_val,max_x,max_y,max_z,delta,stack_start_time,norm_stack_len,output_dir,'stack')

  # WRITE GRID FILES 

  #logging.info("Writing grids...")

  #do_write_grids(stack_grid,0.5,delta,norm_stack_len,stack_start_time,output_dir,'grid')


  # clean_up big memory
  del(stack_grid)


def do_migration_loop_reloc(start_time, end_time, output_dir, kurtosis_filenames, search_grid_filename, time_grid, options_verbose, options_time):


  logging.info("Processing time slice %s"%start_time.isoformat())

  time_dict=time_grid.buf[0]

  # read data into a dictionary

  logging.info("Reading processed data into dictionary")

  if options_time:
    t_ref=time()

  data=read_data_compatible_with_time_dict(kurtosis_filenames,time_dict,start_time,end_time)

  # Set the global variable delta (dt for all the seismograms)
  try:
    delta=data.values()[0].delta
  except IndexError:
    raise UserWarning("File list empty - check --dataglob option")

  # DO MIGRATION
  (max_val,max_x,max_y,max_z,stack_start_time,norm_stack_len,stack_grid)=do_inner_migration_loop(start_time, end_time, data, time_grid, delta, search_grid_filename,options_verbose, options_time)

  # WRITE STACK FILES 
  do_write_stack_files(max_val,max_x,max_y,max_z,delta,stack_start_time,norm_stack_len,output_dir,'reloc_stack')

  # WRITE GRID FILES 

  #logging.info("Writing grids...")

  #do_write_grids(stack_grid,0.1,delta,norm_stack_len,stack_start_time,output_dir,'reloc_grid')


  # clean_up big memory
  del(stack_grid)


def do_migration_loop_plot(start_time, end_time, o_time, grid_dir, kurtosis_filenames, search_grid_filename, time_grid):


  logging.info("Processing time slice %s"%start_time.isoformat())

  time_dict=time_grid.buf[0]

  # read data into a dictionary

  logging.info("Reading processed data into dictionary")


  data=read_data_compatible_with_time_dict(kurtosis_filenames,time_dict,start_time,end_time)

  # Set the global variable delta (dt for all the seismograms)
  try:
    delta=data.values()[0].delta
  except IndexError:
    raise UserWarning("File list empty - check --dataglob option")

  # DO MIGRATION
  (max_val,max_x,max_y,max_z,stack_start_time,norm_stack_len,stack_grid)=do_inner_migration_loop(start_time, end_time, data, time_grid, delta, search_grid_filename)

  # WRITE GRID FILES 

  logging.info("Writing plot grid...")

  grid_filename=do_write_grid_at_time(stack_grid,o_time,delta,norm_stack_len,stack_start_time,grid_dir,'plot_grid')


  # clean_up big memory
  del(stack_grid)

  return grid_filename

def do_write_grid_at_time(stack_grid,o_time,delta,norm_stack_len,stack_start_time,grid_dir,grid_basename):

  itime=np.int(floor((o_time-stack_start_time)/delta))
  timestamp=o_time.isoformat()
  grid_file=os.path.join(grid_dir,"%s_%s.dat"%(grid_basename,timestamp))
  stack_grid.write_grid_timeslice(itime=itime,filename=grid_file)

  return grid_file
