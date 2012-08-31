#!/usr/bin/env python
# encoding: utf-8

import os, sys, optparse, glob, logging
import numpy as np
from obspy.core import *
from obspy.signal import *
from filters import smooth
from locations_trigger import filter_max_stack, number_good_kurtosis_for_location
from grids_paths import StationList, ChannelList, QDTimeGrid, QDGrid
from OP_waveforms import Waveform, read_data_compatible_with_time_dict
from sub_PdF_waveloc import do_innermost_migration_loop
from integrate4D import *
from plot_mpl import *


def trigger_detections(st_max,loc_level):

  logging.debug('Detecting using loc_level = %.2f'%loc_level)
  detections=[]
  for tr in st_max:
    trigs=trigger.triggerOnset(tr.data,loc_level,loc_level)
    logging.debug('Found %d triggers.'%len(trigs))

    for trig in trigs:
      t_start=tr.stats.starttime+trig[0]*tr.stats.delta
      t_end=tr.stats.starttime+trig[1]*tr.stats.delta
      detections.append((t_start,t_end))

  return detections

def compute_stats_from_4Dgrid(opdict,starttime,endtime):
  
  base_path=opdict['base_path']
  datadir=opdict['datadir']
  outdir=opdict['outdir']
  # directories
  lib_path = os.path.join(base_path,'lib')
  data_path= os.path.join(base_path,'data',datadir)
  out_path = os.path.join(base_path,'out',outdir)
  grid_path= os.path.join(base_path,'out',outdir,'grid')
  loc_path = os.path.join(base_path,'out',outdir,'loc')
  fig_path = os.path.join(base_path,'out',outdir,'fig')


  # files
  search_grid=opdict['search_grid']
  time_grid=opdict['time_grid']
  stations=opdict['stations']
  kurt_glob=opdict['kurtglob']
  grad_glob=opdict['gradglob']

  hdr_file = os.path.join(lib_path,search_grid)
  grid_filename_base=os.path.join(lib_path,time_grid)
  stations_filename= os.path.join(lib_path,stations)
  search_grid_filename= hdr_file
  kurt_files=glob.glob(os.path.join(data_path,kurt_glob))
  grad_files=glob.glob(os.path.join(data_path,grad_glob))
  kurt_files.sort()
  grad_files.sort()

  # set up information for migration
  sta=StationList()
  sta.read_from_file(stations_filename)

  cha=ChannelList()
  cha.populate_from_station_list_and_data_files(sta,kurt_files)

  time_grid=QDTimeGrid()
  time_grid.read_NLL_hdr_file(hdr_file)
  time_grid.populate_from_time_grids(grid_filename_base,cha,out_path,load_buf=True)

  max_grid_time=0.0
  for point in time_grid.buf:
    max_time_tmp=max(point.values())
    max_grid_time=max(max_grid_time,max_time_tmp)
  
  # get grid geometry information
  dummy_grid=QDGrid()
  dummy_grid.read_NLL_hdr_file(hdr_file)
  (nx,ny,nz)=(dummy_grid.nx, dummy_grid.ny, dummy_grid.nz)
  (dx,dy,dz)=(dummy_grid.dx, dummy_grid.dy, dummy_grid.dz)
  (x_orig,y_orig,z_orig)=(dummy_grid.x_orig, dummy_grid.y_orig, dummy_grid.z_orig)

  # TODO
  # clean up data (preselect on signal to noise ratios etc)

  logging.info("Processing time slice %s - %s "%(starttime.isoformat(),endtime.isoformat()))
 
  endtime=endtime+2*max_grid_time
  starttime=starttime-2*max_grid_time
  logging.info("Expanding to time slice %s - %s to account for propagation time"%(starttime.isoformat(),endtime.isoformat()))

  time_dict=time_grid.buf[0]

  # read data into a dictionary

  logging.info("Reading processed data into dictionary")
  data=read_data_compatible_with_time_dict(grad_files, time_dict, starttime, endtime)

  # Set the global variable delta (dt for all the seismograms)
  try:
    delta=data.values()[0].delta
  except IndexError:
    raise UserWarning("File list empty - check --dataglob option")

  # DO MIGRATION
  n_buf, norm_stack_len, stack_shift_time, stack_start_time, stack_grid = do_innermost_migration_loop(starttime, endtime, data, time_grid, delta, search_grid_filename)

  # set up integration limits
  x=np.arange(nx)*dx+x_orig
  y=np.arange(ny)*dy+y_orig
  z=np.arange(nz)*dz+z_orig
  t=np.arange(norm_stack_len)*delta
  
  #logging.debug('Expected shape of time axis : %s, actual shape : %s.  Shape of 4D grid : %s.'%(norm_stack_len, x.shape, stack_grid.shape))

  # normalize grid for first probability density calculation
  stack_grid_int=compute_integral4D(stack_grid[:,:,:,0:norm_stack_len],x,y,z,t)
  stack_grid_norm=stack_grid[:,:,:,0:norm_stack_len] / stack_grid_int

 
  # integrate normalized grid over all space dimensions to get marginal over time
  prob_t = si.trapz(si.trapz(si.trapz(stack_grid_norm,x=x,axis=0),x=y,axis=0),x=z,axis=0)
  exp_t = si.trapz(t*prob_t,x=t,axis=0)
  var_t = si.trapz((t-exp_t)*(t-exp_t)*prob_t,x=t,axis=0)
  sigma_t = np.sqrt(var_t)
  it_exp=int(round(exp_t/delta))
  nt_sigma=int(round(sigma_t/delta))
  it_left=it_exp-nt_sigma
  it_right=it_exp+nt_sigma
  t_slice=t[it_left:it_right]
 
  exp_x,exp_y,exp_z,exp_t,cov_matrix,prob_dict = compute_expected_coordinates4D(stack_grid[:,:,:,it_exp-nt_sigma:it_exp+nt_sigma],x,y,z,t_slice,return_2Dgrids=True)
  sigma_x=np.sqrt(cov_matrix[0,0])
  sigma_y=np.sqrt(cov_matrix[1,1])
  sigma_z=np.sqrt(cov_matrix[2,2])
  sigma_t=np.sqrt(cov_matrix[3,3])

  # turn origin time into utcDateTime object
  exp_otime = stack_start_time + exp_t


  loc=(exp_otime,sigma_t,exp_x,sigma_x,exp_y,sigma_y,exp_z,sigma_z)

  # do plotting
  fig_name=os.path.join(fig_path,'fig_st_mpl_%s'%exp_otime.isoformat())
  plot_probloc_mpl(prob_dict,[x,y,z,t_slice],fig_name)

  return loc

def do_locations_prob_setup_and_run(opdict):

  base_path=opdict['base_path']
  outdir=opdict['outdir']
  datadir=opdict['datadir']
  kurtglob=opdict['kurtglob']
  # set up some paths
  data_path= os.path.join(base_path,'data',datadir)
  stack_path=os.path.join(base_path,'out',outdir,'stack')
  loc_path=  os.path.join(base_path,'out',outdir,'loc')
  loc_filename=os.path.join(loc_path,'locations_prob.dat')

  kurt_files=glob.glob(os.path.join(data_path,kurtglob))
  kurt_files.sort()

  logging.info("Path for stack files : %s"%stack_path)
  logging.info("Path for loc files : %s"%loc_path)
  logging.info("Location file : %s"%loc_filename)

  logging.info("Merging and filtering stack files ...")
  st_max=read(os.path.join(stack_path,"stack_max*"))
  st_max.write(os.path.join(stack_path,"combined_stack_max.mseed"),format='MSEED')

  st_max_filt=filter_max_stack(st_max,1.0)
  st_max_filt.write(os.path.join(stack_path,"combined_stack_max_filt.mseed"),format='MSEED')

  loclevel=opdict['loclevel']
  detection_list=trigger_detections(st_max_filt,loclevel)

  logging.debug(detection_list)

  locs=[]
  for (starttime,endtime) in detection_list:
    
    loc=compute_stats_from_4Dgrid(opdict,starttime,endtime)
    locs.append(loc)


  # write to file
  loc_file=open(loc_filename,'w')

  snr_limit=opdict['snr_limit']
  sn_time=opdict['sn_time']
  n_kurt_min=opdict['n_kurt_min']


  n_ok=0
  for (exp_otime,sigma_t,exp_x,sigma_x,exp_y,sigma_y,exp_z,sigma_z) in locs:
    if number_good_kurtosis_for_location(kurt_files,exp_otime,snr_limit,sn_time) > n_kurt_min:
      logging.info("PROB DENSITY : Time %s s pm %.2fs, x=%.4f pm %.4f, y=%.4f pm %.4f, z=%.4f pm %.4f"%(exp_otime.isoformat(),sigma_t, exp_x, sigma_x,exp_y,sigma_y,exp_z,sigma_z))
      loc_file.write("PROB DENSITY : Time %s s pm %.2fs, x=%.4f pm %.4f, y=%.4f pm %.4f, z=%.4f pm %.4f\n"%(exp_otime.isoformat(),sigma_t, exp_x, sigma_x,exp_y,sigma_y,exp_z,sigma_z))
      n_ok=n_ok+1
  loc_file.close()
  logging.info('Wrote %d locations to file %s.'%(n_ok,loc_filename))


if __name__=='__main__':

  from options import WavelocOptions
  logging.basicConfig(level=logging.INFO, format='%(levelname)s : %(asctime)s : %(message)s')

  wo = WavelocOptions()
  args=wo.p.parse_args()

  wo.set_all_arguments(args)
  wo.verify_location_options()

  do_locations_prob_setup_and_run(wo.opdict)


