#!/usr/bin/env python
# encoding: utf-8

import os, sys, optparse, h5py
import glob, logging
import numpy as np

from OP_waveforms import Waveform
from obspy.core import read, UTCDateTime
from locations_trigger import read_locs_from_file
from NllGridLib import read_stations_file,read_hdr_file
from hdf5_grids import get_interpolated_time_grids
from migration import do_migration_loop_continuous
from OP_waveforms import read_data_compatible_with_time_dict
from integrate4D import compute_expected_coordinates4D, \
                        compute_expected_coordinates3D

def read_prob_locs_from_file(filename):
  from obspy.core import utcdatetime

  locs=[]

  f=open(filename,'r')
  lines=f.readlines()
  f.close()

  for line in lines:
    loc={}

    loc['o_time']=utcdatetime.UTCDateTime(line.split()[5])
    loc['o_err']=np.float(line.split()[8])
    loc['x_mean']=np.float(line.split()[11])
    loc['x_sigma']=np.float(line.split()[13])
    loc['y_mean']=np.float(line.split()[16])
    loc['y_sigma']=np.float(line.split()[18])
    loc['z_mean']=np.float(line.split()[21])
    loc['z_sigma']=np.float(line.split()[23])

    locs.append(loc)

  return locs
 
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

def do_locations_prob_setup_and_run(opdict,space_only=True):

  # get / set info
  base_path=opdict['base_path']

  locfile=os.path.join(base_path,'out',opdict['outdir'],'loc','locations.dat')
  locfile_prob=os.path.join(base_path,'out',opdict['outdir'],'loc','locations_prob.dat')
  f_prob=open(locfile_prob,'w')

  # if locfile does not exist then make it by running trigger location
  if not os.path.exists(locfile):
    logging.info('No location found at %s.  Running trigger location first...'%locfile)
    do_locations_trigger_setup_and_run(opdict)

  # directories
  grid_dir=os.path.join(base_path,'out',opdict['outdir'],'grid')
  output_dir=os.path.join(base_path,'out',opdict['outdir'])

  # data files
  data_dir=os.path.join(base_path,'data',opdict['datadir'])
  data_glob=opdict['dataglob']
  kurt_glob=opdict['kurtglob']
  grad_glob=opdict['gradglob']
  data_files=glob.glob(os.path.join(data_dir,data_glob))
  kurt_files=glob.glob(os.path.join(data_dir,kurt_glob))
  grad_files=glob.glob(os.path.join(data_dir,grad_glob))
  data_files.sort()
  kurt_files.sort()
  grad_files.sort()

  # stations
  stations_filename=os.path.join(base_path,'lib',opdict['stations'])
  stations=read_stations_file(stations_filename)

  # grids
  grid_filename_base=os.path.join(base_path,'lib',opdict['time_grid'])
  search_grid_filename=os.path.join(base_path,'lib',opdict['search_grid'])

  # read time grid information
  time_grids=get_interpolated_time_grids(opdict)

  # read locations
  locs=read_locs_from_file(locfile)


  # iterate over locations
  for loc in locs:

    # create the appropriate grid on the fly

    # generate the grids
    o_time=loc['o_time']
    if space_only:
        start_time=o_time
        end_time  =o_time
    else:
        start_time=o_time-3*loc['o_err_left']
        end_time=o_time+3*loc['o_err_right']

    # make a buffer for migration
    start_time_migration = start_time - 10.0
    end_time_migration   =   end_time + 10.0

    # re-read grid info to ensure clean copy
    grid_info=read_hdr_file(search_grid_filename)
 
    # read data
    grad_dict,delta = read_data_compatible_with_time_dict(grad_files,
          time_grids, start_time_migration, end_time_migration)

    # do migration (all metadata on grid is added to grid_info)
    do_migration_loop_continuous(opdict, grad_dict, delta,
          start_time_migration, grid_info, time_grids, keep_grid=True)


    # integrate to get the marginal probability density distributions

    # get required info
    grid_starttime=grid_info['start_time']
    nx,ny,nz,nt=grid_info['grid_shape']
    dx,dy,dz,dt=grid_info['grid_spacing']
    x_orig,y_orig,z_orig=grid_info['grid_orig']

    # we are only interested in the time around the origin time of the event
    it_left  = np.int(np.round((start_time - grid_starttime)/dt))
    it_right = np.int(np.round((end_time   - grid_starttime)/dt))
    it_true  = np.int(np.round((o_time     - grid_starttime)/dt))
    nt=(it_right-it_left)+1

    # set up integration axes (wrt reference)
    x=np.arange(nx)*dx
    y=np.arange(ny)*dy
    z=np.arange(nz)*dz
    if not space_only:
      t=np.arange(nt)*dt

    # open the grid file
    grid_filename=grid_info['dat_file']
    f=h5py.File(grid_filename,'r')
    stack_grid=f['stack_grid']

    # extract the portion of interest (copy data)
    if space_only:
        stack_3D=np.empty((nx,ny,nz))
        stack_3D[:] = stack_grid[:,it_true].reshape(nx,ny,nz)
    else:
        stack_4D=np.empty((nx,ny,nz,nt))
        stack_4D[:] = stack_grid[:,it_left:it_right+1].reshape(nx,ny,nz,nt)

    # close the grid file
    f.close()

    # Get expected values (normalizes grid internally)
    if space_only:
        exp_x, exp_y, exp_z, cov_matrix, prob_dict = \
            compute_expected_coordinates3D(stack_3D,x,y,z,return_2Dgrids=True)
    else:
        exp_x, exp_y, exp_z, exp_t, cov_matrix, prob_dict = \
            compute_expected_coordinates4D(stack_4D,x,y,z,t,return_2Dgrids=True)
    
    # save the marginals to a hdf5 file in loc subdirectory

    # put reference location back
    exp_x = exp_x + x_orig
    exp_y = exp_y + y_orig
    exp_z = exp_z + z_orig
    if space_only:
        exp_t = o_time
    else:
        exp_t = start_time + exp_t

    # extract uncertainties from covariance matrix
    if space_only:
        sig_x,sig_y,sig_z = np.sqrt(np.diagonal(cov_matrix))
        sig_t = (loc['o_err_left']+loc['o_err_right'])/2.
    else:
        sig_x,sig_y,sig_z,sig_t = np.sqrt(np.diagonal(cov_matrix))



    # write the expected values to a plain text locations file

    f_prob.write("PROB DENSITY : T = %s s pm %.2f s, x= %.4f pm %.4f km, \
y= %.4f pm %.4f km, z= %.4f pm %.4f km\n" % (exp_t.isoformat(), sig_t, \
      exp_x, sig_x, exp_y, sig_y, exp_z, sig_z))

  # close location file
  f_prob.close()



if __name__=='__main__':

  from options import WavelocOptions
  logging.basicConfig(level=logging.INFO, format='%(levelname)s : %(asctime)s : %(message)s')

  wo = WavelocOptions()
  args=wo.p.parse_args()

  wo.set_all_arguments(args)
  wo.verify_location_options()

  do_locations_prob_setup_and_run(wo.opdict)


