import os, sys, optparse
import logging, glob
import numpy as np

from OP_waveforms import Waveform
from obspy.core import read
from locations_trigger import read_locs_from_file
from NllGridLib import read_stations_file,read_hdr_file
from hdf5_grids import get_interpolated_time_grids
from migration import do_migration_loop_continuous
from OP_waveforms import read_data_compatible_with_time_dict
from plot_mpl import plotLocationGrid

def do_plotting_setup_and_run(opdict):

  # get / set info
  base_path=opdict['base_path']

  locfile=os.path.join(base_path,'out',opdict['outdir'],'loc','locations.dat')
  grid_dir=os.path.join(base_path,'out',opdict['outdir'],'grid')
  output_dir=os.path.join(base_path,'out',opdict['outdir'])

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
  
  #for loc in locs[0:1]:
  for loc in locs:
    # generate the grids
    print loc
    o_time=loc['o_time']
    start_time=o_time-opdict['plot_tbefore']
    end_time=o_time+opdict['plot_tafter']

    # TODO implement a rough estimation of the stack shift based on propagation time across the whole network
    start_time_migration=start_time-10.0
    end_time_migration=end_time+10.0


#    # make dictionary of station names with snr on kurtosis ratios
#    # TODO Fix this to estimate K-time from o_time and propagation time to station
#    snr_start_time=o_time-opdict['sn_time']
#    snr_end_time=o_time+opdict['sn_time']
#    snr_dict={}
#    wf=Waveform()
#    for filename in kurt_files:
#      wf.read_from_file(filename,starttime=snr_start_time,endtime=snr_end_time)
#      snr=wf.get_snr(o_time,snr_start_time,snr_end_time)
#      station_name=wf.trace.stats.station
#      snr_dict[station_name]=snr
 
    # set output filename
#    plot_filename=os.path.join(figdir,"loc_%s.png"%o_time.isoformat())

#    # select grad files for which the snr is > snr_limit
#    grad_files_selected=[]
#    for filename in grad_files:
#      st=read(filename,headonly=True)
#      station_name=st.traces[0].stats.station
#      logging.debug("Checking station %s : snr_value = %.2f"%(station_name,snr_dict[station_name]))
#      if snr_dict[station_name]>np.float(opdict['snr_limit']):
#        grad_files_selected.append(filename)

    # read data
    data,delta = read_data_compatible_with_time_dict(grad_files,
            time_grids, start_time_migration, end_time_migration)

    # re-read grid info to ensure clean copy
    grid_info=read_hdr_file(search_grid_filename)

    # do migration
    do_migration_loop_continuous(opdict, data, delta,
            start_time_migration, grid_info, time_grids, keep_grid=True)

    # plot
    plotLocationGrid(loc,grid_info,figdir)


if __name__ == '__main__':

  from options import WavelocOptions
  logging.basicConfig(level=logging.INFO, format='%(levelname)s : %(asctime)s : %(message)s')

  wo = WavelocOptions()
  args=wo.p.parse_args()

  wo.set_all_arguments(args)
  wo.verify_plotting_options()

  do_plotting_setup_and_run(wo.opdict)


