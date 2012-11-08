import os, sys, optparse, h5py
import logging, glob
import numpy as np

from OP_waveforms import Waveform
from obspy.core import read, UTCDateTime
from locations_trigger import read_locs_from_file
from locations_prob import read_prob_locs_from_file
from NllGridLib import read_stations_file,read_hdr_file
from hdf5_grids import get_interpolated_time_grids
from migration import do_migration_loop_continuous
from OP_waveforms import read_data_compatible_with_time_dict
from plot_mpl import plotLocationGrid, plotLocationWaveforms, plotProbLoc

def do_plotting_setup_and_run(opdict,plot_wfm=True,plot_grid=True):

  # get / set info
  base_path=opdict['base_path']

  locfile=os.path.join(base_path,'out',opdict['outdir'],'loc','locations.dat')
  stackfile=os.path.join(base_path,'out',opdict['outdir'],'stack','combined_stack_all.hdf5')
  grid_dir=os.path.join(base_path,'out',opdict['outdir'],'grid')
  output_dir=os.path.join(base_path,'out',opdict['outdir'])

  data_dir=os.path.join(base_path,'data',opdict['datadir'])

  data_glob=opdict['dataglob']
  data_files=glob.glob(os.path.join(data_dir,data_glob))
  data_files.sort()

  kurt_glob=opdict['kurtglob']
  kurt_files=glob.glob(os.path.join(data_dir,kurt_glob))
  kurt_files.sort()
  mig_files=kurt_files

  if opdict['kderiv']:
    grad_glob=opdict['gradglob']
    grad_files=glob.glob(os.path.join(data_dir,grad_glob))
    grad_files.sort()
    mig_files=grad_files

    if opdict['gauss']:
      gauss_glob=opdict['gaussglob']
      gauss_files=glob.glob(os.path.join(data_dir,gauss_glob))
      gauss_files.sort()
      mig_files=gauss_files

  figdir=os.path.join(base_path,'out',opdict['outdir'],'fig')

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

  # open stack file
  f_stack=h5py.File(stackfile,'r')
  max_val=f_stack['max_val_smooth']
  stack_start_time=UTCDateTime(max_val.attrs['start_time'])
  
  for loc in locs:
    # generate the grids
    o_time=loc['o_time']
    start_time=o_time-opdict['plot_tbefore']
    end_time=o_time+opdict['plot_tafter']

    # re-read grid info to ensure clean copy
    grid_info=read_hdr_file(search_grid_filename)
    nx=grid_info['nx']
    ny=grid_info['ny']
    nz=grid_info['nz']
    dx=grid_info['dx']
    dy=grid_info['dy']
    dz=grid_info['dz']

    x=loc['x_mean']
    y=loc['y_mean']
    z=loc['z_mean']
    # get the corresponding travel-times for time-shifting
    ttimes={}
    for sta in time_grids.keys():
        ttimes[sta]=time_grids[sta].value_at_point(x,y,z)

    tshift_migration=max(ttimes.values())

    start_time_migration=start_time-tshift_migration
    end_time_migration=end_time+tshift_migration

    if plot_grid:
      logging.info('Plotting grid for location %s'%o_time.isoformat())
      # TODO implement a rough estimation of the stack shift based on propagation time across the whole network

      # read data
      mig_dict,delta = read_data_compatible_with_time_dict(mig_files,
            time_grids, start_time_migration, end_time_migration)

      # do migration
      do_migration_loop_continuous(opdict, mig_dict, delta,
            start_time_migration, grid_info, time_grids, keep_grid=True)

      # plot
      plotLocationGrid(loc,grid_info,figdir,opdict['plot_otime_window'])

    if plot_wfm:

      logging.info('Plotting waveforms for location %s'%o_time.isoformat())

      # get the index of the location
#      ix=np.int(np.round((loc['x_mean']-grid_info['x_orig'])/dx))
#      iy=np.int(np.round((loc['y_mean']-grid_info['y_orig'])/dy))
#      iz=np.int(np.round((loc['z_mean']-grid_info['z_orig'])/dz))
#      ib= ix*ny*nz + iy*nz + iz

      # read data
      data_dict,delta = read_data_compatible_with_time_dict(data_files,
            time_grids, start_time_migration, end_time_migration)
      grad_dict,delta = read_data_compatible_with_time_dict(grad_files,
            time_grids, start_time_migration, end_time_migration)

      # cut desired portion out of data
      for sta in data_dict.keys():
          tmp=data_dict[sta]
          istart=np.int(np.round(
              (start_time + ttimes[sta] - start_time_migration) / delta))
          iend=istart + np.int(np.round(
              (opdict['plot_tbefore'] + opdict['plot_tafter'])  / delta))
          # sanity check in case event is close to start or end of data
          if istart < 0 : istart=0
          if iend   > len(tmp) : iend = len(tmp)
          data_dict[sta]=tmp[istart:iend]
          # do slice
          tmp=grad_dict[sta]
          grad_dict[sta]=tmp[istart:iend]

      # retrieve relevant portion of stack max
      istart=np.int(np.round(
          (o_time - opdict['plot_tbefore'] -stack_start_time) / delta))
      iend=istart + np.int(np.round(
          (opdict['plot_tbefore'] + opdict['plot_tafter'])  / delta))
      # sanity check in case event is close to start or end of data
      if istart < 0 : 
          start_time = start_time + np.abs(istart)*dt
          istart=0
      if iend   > len(max_val) : iend = len(max_val)
      # do slice
      stack_wfm=max_val[istart:iend]

      # plot
      plotLocationWaveforms(loc,start_time,delta,data_dict,grad_dict,stack_wfm,figdir)

  f_stack.close()


def do_probloc_plotting_setup_and_run(opdict):

  # get / set info
  base_path=opdict['base_path']
  space_only=opdict['probloc_spaceonly']

  locfile=os.path.join(base_path,'out',opdict['outdir'],'loc','locations.dat')
  problocfile=os.path.join(base_path,'out',opdict['outdir'],'loc','locations_prob.dat')
  problocgrid=os.path.join(base_path,'out',opdict['outdir'],'loc','locations_prob.hdf5')
  output_dir=os.path.join(base_path,'out',opdict['outdir'])


  figdir=os.path.join(base_path,'out',opdict['outdir'],'fig')

  # read locations
  locs=read_locs_from_file(locfile)
  prob_locs=read_prob_locs_from_file(problocfile)

  # open hdf5 file
  f = h5py.File(problocgrid,'r')

  # for each loc
  for i in xrange(len(locs)):
    loc      =      locs[i]
    prob_loc = prob_locs[i]

    otime=prob_loc['o_time']
    # hdf5 group name is prob loc origin time as string
    grp = f[prob_loc['o_time'].isoformat()]

    # do plotting
    plotProbLoc(grp,prob_loc,loc,figdir, space_only)

 
  # close hdf5 file
  f.close()

if __name__ == '__main__':

  from options import WavelocOptions
  logging.basicConfig(level=logging.INFO, format='%(levelname)s : %(asctime)s : %(message)s')

  wo = WavelocOptions()
  args=wo.p.parse_args()

  wo.set_all_arguments(args)
  wo.verify_plotting_options()

  do_plotting_setup_and_run(wo.opdict)


