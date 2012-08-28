#!/usr/bin/env python
# encoding: utf-8

import os, sys, optparse, glob, logging
from obspy.core import *
from obspy.signal import *
from filters import smooth
from locations_trigger import filter_max_stack, number_good_kurtosis_for_location


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

def compute_stats_from_4Dgrid(base_path,outdir,datadir,search_grid,time_grid,stations,data_glob,kurt_glob,grad_glob,starttime,endtime):
  
  # directories
  aux_path = os.path.join(base_path,'aux')
  data_path= os.path.join(base_path,'data',datadir)
  grid_path= os.path.join(base_path,'out',outdir,'grid')
  loc_path = os.path.join(base_path,'out',outdir,'loc')

  if not os.path.exists(grid_path):
    os.makedirs(grid_path)
  if not os.path.exists(loc_path):
    os.makedirs(loc_path)

  # files
  hdr_file = os.path.join(aux_path,search_grid)
  grid_filename_base=os.path.join(aux_path,time_grid)
  stations_filename= os.path.join(aux_path,stations)
  search_grid_filename= os.path.join(aux_path,search_grid)
  data_files=glob.glob(os.path.join(data_path,data_glob))
  kurt_files=glob.glob(os.path.join(data_path,kurt_glob))
  grad_files=glob.glob(os.path.join(data_path,grad_glob))
  data_files.sort()
  kurt_files.sort()
  grad_files.sort()

  # set up information for migration
  sta=StationList()
  sta.read_from_file(stations_filename)

  cha=ChannelList()
  cha.populate_from_station_list_and_data_files(sta,data_files)

  time_grid=QDTimeGrid()
  time_grid.read_NLL_hdr_file(hdr_file)
  time_grid.populate_from_time_grids(grid_filename_base,cha,load_buf=True)
  
  # get grid geometry information
  dummy_grid=QDGrid()
  dummy_grid.read_NLL_hdr_file(hdr_file)
  (nx,ny,nz)=(dummy_grid.nx, dummy_grid.ny, dummy_grid.nz)
  (dx,dy,dz)=(dummy_grid.dx, dummy_grid.dy, dummy_grid.dz)
  (x_orig,y_orig,z_orig)=(dummy_grid.x_orig, dummy_grid.y_orig, dummy_grid.z_orig)

  # TODO
  # clean up data (preselect on signal to noise ratios etc)

  logging.info("Processing time slice %s - %s "%starttime.isoformat(),endtime.isoformat())

  time_dict=time_grid.buf[0]

  # read data into a dictionary

  logging.info("Reading processed data into dictionary")

  data={}

  for filename in grad_files:
    wf=Waveform()
    try:
      # read will return UserWarning if there is no data within start and end time
      # will pad blanks with zeros if required (no tapering applied, as kurtosis files are already correctly tapered to zero)
      wf.read_from_file(filename,starttime=start_time,endtime=end_time,pad_value=0)
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
  xt=np.arange(norm_stack_len)*delta
  x0=np.arange(nx)*dx
  x1=np.arange(ny)*dy
  x2=np.arange(nz)*dz
  
  print xt.shape(), stack_grid.shape()


def do_locations_prob_setup_and_run(base_path="",outdir="",loclevel=None,datadir="",dataglob="",snr_limit=None,sn_time=None,n_kurt_min=None,search_grid=""):

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

  detection_list=trigger_detections(st_max_filt,loclevel)

  print detection_list

  for (starttime,endtime) in detection_list:
    
    compute_stats_from_4Dgrid(base_path,outdir,datadir,search_grid,time_grid,stations,data_glob,kurt_glob,grad_glob,starttime,endtime)

if __name__=='__main__':

  logging.basicConfig(level=logging.DEBUG, format='%(levelname)s : %(asctime)s : %(message)s')
  base_path=os.getenv('WAVELOC_PATH')

  # Read command line

  p = optparse.OptionParser()
  p.add_option('--outdir', '-o', action='store', help='output subdirectory in which the stack directory is found')
  p.add_option('--loclevel', action='store', default=100, type='float',help='trigger stack level for locations (e.g. 100) ')
  p.add_option('--datadir',action='store',help="data subdirectory")
  p.add_option('--dataglob',action='store',help="data glob")
  p.add_option('--snr_limit',action='store',default=10.0, type='float',help="signal_to_noise level for kurtosis acceptance")
  p.add_option('--sn_time',action='store',default=10.0, type='float',help="time over which to calculate the signal_to_noise ratio for kurtosis acceptance")
  p.add_option('--n_kurt_min',action='store',default=4, type='int',help="min number of good kurtosis traces for a location")
  p.add_option('--search_grid',action='store',help="search grid %s")
  p.add_option('--stations','-s',action='store',default='channels_HHZ.dat',help='station list (found in $WAVELOC_PATH/aux)')
  p.add_option('--data_glob',action='store',help="data glob")
  p.add_option('--kurt_glob',action='store',help="kurtosis glob")
  p.add_option('--grad_glob',action='store',help="kurtosis gradient glob")
  p.add_option('--snr_limit',action='store',default=10.0, help="signal_to_noise level for kurtosis acceptance")
  p.add_option('--sn_time',action='store',default=10.0, help="time over which to calculate the signal_to_noise ratio for kurtosis acceptance")

  (options,arguements)=p.parse_args()

  do_locations_prob_setup_and_run(base_path=base_path,outdir=options.outdir,loclevel=options.loclevel,datadir=options.datadir,dataglob=options.dataglob,snr_limit=options.snr_limit,sn_time=options.sn_time,n_kurt_min=options.n_kurt_min)
