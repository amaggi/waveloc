#!/usr/bin/env python
# encoding: utf-8

import os, sys, h5py, tempfile

import numpy as np
#import numexpr as ne
#import carray as ca

from OP_waveforms import *

from hdf5_grids import migrate_4D_stack, extract_max_values
from NllGridLib import read_hdr_file
from time import time, sleep
from scipy import weave
from scipy.weave import converters
import logging

def do_innermost_migration_loop(start_time, end_time, data, time_grids, delta, grid_info,options_verbose=False, options_time=False):

  logging.warn('Deprecated - do not use')

  ######### DO THE MIGRATION #############

  logging.info("Stacking shifted time series... ")

  if options_time:
    t_ref=time()  

  (n_buf, norm_stack_len, stack_shift_time, stack_grid) = migrate_4D_stack(data, delta, grid_info, time_grids)
  stack_start_time=start_time-stack_shift_time

  logging.debug("Stack geographical dimension = %d"%n_buf)
  logging.debug("Stack time extent = %d points = %.2f s"%(norm_stack_len, norm_stack_len*delta))
  logging.debug("Start time of stack (wrt start time of data)= %.2f s"%(-stack_shift_time))
  logging.debug("Start time of stack %s"%(stack_start_time.isoformat()))

  if options_time:
    t=time()-t_ref
    logging.info("Time for stacking and saving %d stacks, each of extent %d points : %.2f s\n" % (n_buf,norm_stack_len,t))
 
  return n_buf, norm_stack_len, stack_shift_time, stack_start_time, stack_grid

def do_inner_migration_loop(start_time, end_time, data, time_grids, delta, grid_info, options_verbose=False, options_time=False):

  
  logging.warn('Deprecated - do not use')
 
  (n_buf, norm_stack_len, stack_shift_time, stack_start_time, stack_grid)= do_innermost_migration_loop(start_time, end_time, data, time_grids, delta, grid_info,options_verbose, options_time)

  ###### Extract maximum of stack #######


  logging.info("Extracting maximum of stack")

  if options_time:
    t_ref=time()  

  max_val,max_x,max_y,max_z=extract_max_values(stack_grid,grid_info)


  if options_time:
    t=time()-t_ref
    logging.info("Time for extracting maxima from %d slices : %.2f s\n" % (norm_stack_len,t))


  return(max_val,max_x,max_y,max_z,stack_start_time,norm_stack_len,stack_grid)
 

def do_write_stack_files(max_val,max_x,max_y,max_z,delta,stack_start_time,norm_stack_len,output_dir,stack_basename):
 
  logging.warn('Deprecated - do not use')
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

def do_write_hdf5_stack_files(max_val,max_x,max_y,max_z,delta,stack_start_time,norm_stack_len,output_dir,stack_basename):
 
  ######## Transform results into waveform ###########

  logging.info("Writing stack files as hdf5 files ...")

  # create headers
  stats_stack={'network': 'UV', 'station' : 'STACK', 'location': '', 'channel': 'STK', 'npts' : norm_stack_len, 'sampling_rate': 1/delta}
  stats_x={'network': 'UV', 'station' : 'XYZ', 'location': '', 'channel': 'X', 'npts' : norm_stack_len, 'sampling_rate': 1/delta}
  stats_y={'network': 'UV', 'station' : 'XYZ', 'location': '', 'channel': 'Y', 'npts' : norm_stack_len, 'sampling_rate': 1/delta}
  stats_z={'network': 'UV', 'station' : 'XYZ', 'location': '', 'channel': 'Z', 'npts' : norm_stack_len, 'sampling_rate': 1/delta}

  # set time
  stats_stack['starttime']=stack_start_time.isoformat()
  stats_x['starttime']=stack_start_time.isoformat()
  stats_y['starttime']=stack_start_time.isoformat()
  stats_z['starttime']=stack_start_time.isoformat()

  stack_file=os.path.join(output_dir,'stack',"%s_all_%s.hdf5"%(stack_basename,stack_start_time.isoformat()))
  f=h5py.File(stack_file,'w')
  f.create_dataset('max_val',data=max_val)
  f.create_dataset('max_x',data=max_x)
  f.create_dataset('max_y',data=max_y)
  f.create_dataset('max_z',data=max_z)

  f.close()



def do_write_grids(stack_grid,time_step_sec,delta,norm_stack_len,stack_start_time,output_dir,grid_basename):

  logging.warn('Deprecated - do not use')
  # write grids every time_step_sec seconds
  time_step=int(floor((time_step_sec/delta)))
  itimes=numpy.arange(0,norm_stack_len,time_step)

  for itime in itimes:
    slicetime=stack_start_time+itime*delta
    timestamp=slicetime.isoformat()
    grid_file=os.path.join(output_dir,'grid',"%s_%s.dat"%(grid_basename,timestamp))
    stack_grid.write_grid_timeslice(itime=itime,filename=grid_file)

 

def do_migration_loop_reloc(start_time, end_time, output_dir, kurtosis_filenames, grid_info, time_grids, options_verbose, options_time):


  logging.info("Processing time slice %s"%start_time.isoformat())


  # read data into a dictionary

  logging.info("Reading processed data into dictionary")

  if options_time:
    t_ref=time()

  data=read_data_compatible_with_time_dict(kurtosis_filenames,time_grids,start_time,end_time)

  # Set the global variable delta (dt for all the seismograms)
  try:
    delta=data.values()[0].delta
  except IndexError:
    raise UserWarning("File list empty - check --dataglob option")

  # DO MIGRATION
  (max_val,max_x,max_y,max_z,stack_start_time,norm_stack_len,stack_grid)=do_inner_migration_loop(start_time, end_time, data, time_grids, delta, grid_info,options_verbose, options_time)

  # WRITE STACK FILES 
  #do_write_stack_files(max_val,max_x,max_y,max_z,delta,stack_start_time,norm_stack_len,output_dir,'reloc_stack')
  do_write_hdf5_stack_files(max_val,max_x,max_y,max_z,delta,stack_start_time,norm_stack_len,output_dir,'reloc_stack')

  # WRITE GRID FILES 

  #logging.info("Writing grids...")

  #do_write_grids(stack_grid,0.1,delta,norm_stack_len,stack_start_time,output_dir,'reloc_grid')


  # clean_up big memory
  del(stack_grid)


def do_migration_loop_plot(start_time, end_time, o_time, grid_dir, kurtosis_filenames, grid_info, time_grids, write=False):

  logging.info("Processing time slice %s"%start_time.isoformat())


  # read data into a dictionary

  logging.info("Reading processed data into dictionary")


  data=read_data_compatible_with_time_dict(kurtosis_filenames,time_grids,start_time,end_time)

  # Set the global variable delta (dt for all the seismograms)
  try:
    delta=data.values()[0].delta
  except IndexError:
    raise UserWarning("File list empty - check --dataglob option")

  # DO MIGRATION
  (n_buf, norm_stack_len, stack_shift_time, stack_start_time, stack_grid)=do_innermost_migration_loop(start_time, end_time, data, time_grids, delta, grid_info)
  logging.info(o_time)
  logging.info(stack_start_time)
  logging.info(stack_shift_time)

  # WRITE GRID FILES 

  logging.info("Writing plot grid...")
  timestamp=o_time.isoformat()
  grid_file=os.path.join(grid_dir,"%s_%s.dat"%('plot_grid',timestamp))

  if write : stack_grid[:,:,:,0:norm_stack_len].tofile(grid_file)

  # set up information
  stack_grid_info={}
  stack_grid_info['dat_file']=grid_file
  stack_grid_info['grid_shape']=stack_grid[:,:,:,0:norm_stack_len].shape
  stack_grid_info['grid_spacing']=time_grid.dx,time_grid.dy,time_grid.dz,delta
  stack_grid_info['grid_orig']=time_grid.x_orig,time_grid.y_orig,time_grid.z_orig
  stack_grid_info['stack_shift_time']=stack_shift_time
  stack_grid_info['stack_starttime']=stack_start_time
  stack_grid_info['stack_otime']=o_time
 
  logging.info(grid_info)

  # WRITE INFO FILE
  info_file=os.path.join(grid_dir,"%s_%s.info"%('plot_gridinfo',timestamp))
  f=open(info_file,'w')
  f.write(str(grid_info))


  # clean_up big memory
#  del(stack_grid)

  return stack_grid_info, stack_grid[:,:,:,0:norm_stack_len]

def do_write_grid_at_time(stack_grid,o_time,delta,norm_stack_len,stack_start_time,grid_dir,grid_basename):

  itime=np.int(floor((o_time-stack_start_time)/delta))
  timestamp=o_time.isoformat()
  grid_file=os.path.join(grid_dir,"%s_%s.dat"%(grid_basename,timestamp))
  stack_grid.write_grid_timeslice(itime=itime,filename=grid_file)

  return grid_file
