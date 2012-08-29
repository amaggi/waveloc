#!/usr/bin/env python
# encoding: utf-8

import os, sys, optparse, glob, logging
import numpy as np
from obspy.core import *
from obspy.signal import *
from filters import smooth
from locations_trigger import filter_max_stack, number_good_kurtosis_for_location
from grids_paths import StationList, ChannelList, QDTimeGrid, QDGrid
from OP_waveforms import Waveform
from sub_PdF_waveloc import do_innermost_migration_loop
from integrate4D import *


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
  aux_path = os.path.join(base_path,'aux')
  data_path= os.path.join(base_path,'data',datadir)
  grid_path= os.path.join(base_path,'out',outdir,'grid')
  loc_path = os.path.join(base_path,'out',outdir,'loc')


  # files
  search_grid=opdict['search_grid']
  time_grid=opdict['time_grid']
  stations=opdict['stations']
  kurt_glob=opdict['kurtglob']
  grad_glob=opdict['gradglob']

  hdr_file = os.path.join(aux_path,search_grid)
  grid_filename_base=os.path.join(aux_path,time_grid)
  stations_filename= os.path.join(aux_path,stations)
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
  time_grid.populate_from_time_grids(grid_filename_base,cha,load_buf=True)

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

  data={}

  for filename in grad_files:
    wf=Waveform()
    try:
      # read will return UserWarning if there is no data within start and end time
      # will pad blanks with zeros if required (no tapering applied, as kurtosis files are already correctly tapered to zero)
      wf.read_from_file(filename,starttime=starttime,endtime=endtime,pad_value=0)
      wf_id="%s.%s"%(wf.station,wf.comp)
      # if all is ok, and we have a corresponding time id, add data to dictionary
      if time_dict.has_key(wf_id):
        data[wf_id]=wf
      else:
        logging.info('Station %s not present in time_grid.  Ignoring station...'%wf_id)
    except UserWarning,msg:
      # for any UserWarning, ignore data
      logging.error("No data data found between limits for file %s. Ignore station."%filename)

  # Set the global variable delta (dt for all the seismograms)
  try:
    delta=data.values()[0].delta
  except IndexError:
    raise UserWarning("File list empty - check --dataglob option")

  # DO MIGRATION
  n_buf, norm_stack_len, stack_shift_time, stack_start_time, stack_grid = do_innermost_migration_loop(starttime, endtime, data, time_grid, delta, search_grid_filename)

  # set up integration limits
  x0=np.arange(nx)*dx
  x1=np.arange(ny)*dy
  x2=np.arange(nz)*dz
  xt=np.arange(norm_stack_len)*delta
  
  logging.debug('Expected shape of time axis : %s, actual shape : %s.  Shape of 4D grid : %s.'%(norm_stack_len, xt.shape, stack_grid.shape))
  exp_x0,exp_x1,exp_x2,exp_xt = compute_expected_coordinates4D(stack_grid[:,:,:,0:norm_stack_len],x0,x1,x2,xt)

  exp_x0 += x_orig
  exp_x1 += y_orig
  exp_x2 += z_orig
  exp_otime = stack_start_time + exp_xt

  logging.info('Located event at %s, x = %.3f, y = %.3f, z = %.3f'%(exp_otime,exp_x0,exp_x1,exp_x2))

def do_locations_prob_setup_and_run(opdict):

  base_path=opdict['base_path']
  outdir=opdict['outdir']
  # set up some paths
  stack_path=os.path.join(base_path,'out',outdir,'stack')
  loc_path=os.path.join(base_path,'out',outdir,'loc')
  loc_filename=os.path.join(loc_path,'locations_prob.dat')

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

  for (starttime,endtime) in detection_list:
    
    compute_stats_from_4Dgrid(opdict,starttime,endtime)

if __name__=='__main__':

  from options import WavelocOptions
  logging.basicConfig(level=logging.INFO, format='%(levelname)s : %(asctime)s : %(message)s')

  wo = WavelocOptions()
  args=wo.p.parse_args()

  wo.set_all_arguments(args)
  wo.verify_location_options()

  do_locations_prob_setup_and_run(wo.opdict)


