#!/usr/bin/env python
# encoding: utf-8

import os, sys, optparse, glob
from obspy.core import *
import matplotlib.pyplot as plt
import numpy as np
from OP_waveforms import Waveform
from grids_paths import QDGrid, StationList, ChannelList, QDTimeGrid
from sub_PdF_waveloc import do_migration_loop_plot
from PIL import Image
from plot_slice_mayavi import plot_slice_mayavi
from NLL_IO import qd_read_hyp_file, qd_read_picks_from_hyp_file
import logging


logging.basicConfig(level=logging.DEBUG, format='%(levelname)s : %(asctime)s : %(message)s')

# get path
base_path=os.getenv('WAVELOC_PATH')

# Read command line


p=optparse.OptionParser()
p.add_option('--run_mayavi', action='store_true', help="use mayavi to re-plot grid slices")
p.add_option('--outdir', '-o', action='store', help='output subdirectory')
p.add_option('--datadir', '-d', action='store', help='data subdirectory')
p.add_option('--2D',action='store_true',default=False,dest='twoD',help='use 2D time grids')
p.add_option('--loc_picks_dir', action='store', help='subdirectory for location picks')
p.add_option('--search_grid',action='store',help="search grid %s")
p.add_option('--stations','-s',action='store',default='channels_HHZ.dat',help='station list (found in $WAVELOC_PATH/lib)')
p.add_option('--data_glob',action='store',help="data glob")
p.add_option('--kurt_glob',action='store',help="kurtosis glob")
p.add_option('--grad_glob',action='store',help="kurtosis gradient glob")
p.add_option('--time_grid',action='store',help="time grid")
p.add_option('--hyp_glob',action='store',help="hypocenter pick glob")
#p.add_option('--reloc', action='store_true', default=False, help='apply to relocated events')
p.add_option('--max_stack', action='store', default=200, help='stack value at which to saturate color scale')
p.add_option('--snr_limit',action='store',default=10.0, help="signal_to_noise level for kurtosis acceptance")
p.add_option('--sn_time',action='store',default=10.0, help="time over which to calculate the signal_to_noise ratio for kurtosis acceptance")

(options,arguments)=p.parse_args()

#do_reloc=options.reloc
max_stack_value=np.float(options.max_stack)

if options.loc_picks_dir==None or options.hyp_glob==None:
  do_hyp = False
else:
  do_hyp = True

lib_path=os.path.join(base_path,'lib')
stack_path=os.path.join(base_path,'out',options.outdir,'stack')
grid_path= os.path.join(base_path,'out',options.outdir,'grid')
out_path= os.path.join(base_path,'out',options.outdir)

stack_file=os.path.join(stack_path,"combined_stack_max_filt.mseed")
loc_path=os.path.join(base_path,'out',options.outdir,'loc')
loc_filename=os.path.join(loc_path,'locations.dat')

data_path=os.path.join(base_path,'data',options.datadir)
hdr_file= os.path.join(lib_path,options.search_grid)
grid_filename_base=os.path.join(lib_path,options.time_grid)
stations_filename=os.path.join(lib_path,options.stations)

if do_hyp:
  hyp_path=os.path.join(base_path,'data',options_loc_picks_dir,'loc_hyp')


data_files=glob.glob(os.path.join(data_path,options.data_glob))
kurt_files=glob.glob(os.path.join(data_path,options.kurt_glob))
grad_files=glob.glob(os.path.join(data_path,options.grad_glob))
data_files.sort()
kurt_files.sort()
grad_files.sort()
if do_hyp:
  hyp_files=glob.glob(os.path.join(hyp_path,options.hyp_glob))

#  ***** reading station file ******

sta=StationList()
sta.read_from_file(stations_filename)


cha=ChannelList()
cha.populate_from_station_list_and_data_files(sta,data_files)

######### INTERPOLATE TRAVEL TIMES #############

# The time grid will contain as array values just the travel-times needed 
# (interpolated from the full NLL files) so we can free up the memory as soon as possible

time_grid=QDTimeGrid()
time_grid.read_NLL_hdr_file(hdr_file)
if options.twoD:
  time_grid.populate_from_2D_time_grids(grid_filename_base,cha)
else:
  time_grid.populate_from_time_grids(grid_filename_base,cha,out_path,load_buf=True)


print "Getting Grid geometry"
dummy_grid=QDGrid()
dummy_grid.read_NLL_hdr_file(hdr_file)
(nx,ny,nz)=(dummy_grid.nx, dummy_grid.ny, dummy_grid.nz)

if do_hyp:
  logging.info("Reading hyp parameters")
  hyp_parameters=[]
  for hyp_file in hyp_files:
    logging.debug("Hyp info from %s"%hyp_file)
    (otime,hypo_x,sigma_x,hypo_y,sigma_y,hypo_z,sigma_z)=qd_read_hyp_file(hyp_file)
    phase_dict=qd_read_picks_from_hyp_file(hyp_file)
    hyp_parameters.append(((otime,hypo_x,sigma_x,hypo_y,sigma_y,hypo_z,sigma_z),phase_dict))

  info.debug(hyp_parameters[0])

# read location file
loc_file=open(loc_filename,'r')
loc_lines=loc_file.readlines()
loc_file.close()

locs=[]
for line in loc_lines:
  words=line.split()
  stack_max=np.float(words[2].split(',')[0])
  stack_time=utcdatetime.UTCDateTime(words[3])
  stack_time_err_left=np.float(words[5])
  stack_time_err_right=np.float(words[8])
  stack_x=np.float(words[11])
  stack_x_err=np.float(words[13])
  stack_y=np.float(words[16])
  stack_y_err=np.float(words[18])
  stack_z=np.float(words[21])
  stack_z_err=np.float(words[23])
  locs.append((stack_max,stack_time,stack_time_err_left,stack_time_err_right,stack_x,stack_x_err,stack_y,stack_y_err,stack_z,stack_z_err))


#for loc in locs[0:4]:
for loc in locs:

  print loc

  stack_time=loc[1]
  stack_time_err_left=loc[2]
  stack_time_err_right=loc[3]
  stack_x=loc[4]
  stack_x_err=loc[5]
  stack_y=loc[6]
  stack_y_err=loc[7]
  stack_z=loc[8]
  stack_z_err=loc[9]
  # set start and end time of plot
  start_time=stack_time-10.0
  end_time=stack_time+20.0

  start_time_migration=stack_time-60.0
  end_time_migration=stack_time+60.0


  # make dictionary of station names with snr ratios
  snr_start_time=stack_time-options.sn_time
  snr_end_time=stack_time+options.sn_time
  snr_dict={}
  wf=Waveform()
  for filename in kurt_files:
    wf.read_from_file(filename,starttime=snr_start_time,endtime=snr_end_time)
    snr=wf.get_snr(stack_time,snr_start_time,snr_end_time)
    station_name=wf.trace.stats.station
    snr_dict[station_name]=snr

  logging.debug("Signal to noise ratios on kurtosis")
  logging.debug(snr_dict)

  # set output filename
  plot_filename=os.path.join(loc_path,"loc_%s.png"%stack_time.isoformat())

  # select grad files for which the snr is > snr_limit
  grad_files_selected=[]
  for filename in grad_files:
    st=read(filename,headonly=True)
    station_name=st.traces[0].stats.station
    logging.debug("Checking station %s : snr_value = %.2f"%(station_name,snr_dict[station_name]))
    if snr_dict[station_name]>np.float(options.snr_limit):
      grad_files_selected.append(filename)

  grid_name=do_migration_loop_plot(start_time_migration,end_time_migration,stack_time,grid_path,grad_files_selected,hdr_file,time_grid)
  # Set png name to give to plot_slice_mayavi
  png_name="%s.png"%grid_name

  if do_hyp:
    # find hyp info
    A=[(np.abs(hyp_tuple[0][0]-stack_time),hyp_tuple) for hyp_tuple in hyp_parameters ]
    A.sort()
    hyp_tuple=A[0][1]
    (otime, hypo_x, sigma_x, hypo_y, sigma_y, hypo_z, sigma_z)=hyp_tuple[0]
    phases=hyp_tuple[1]

  # if the closest picked time is too far away, discard the picked location as no longer relevant
  if not do_hyp or (np.abs(otime-stack_time)>10.0) : 
    hypo_x=None  
    hypo_y=None  
    hypo_z=None  
    phses=None


  # create .png file using mayavi
  if options.run_mayavi:
    print "Creating plot using mayavi"
    (x_data,y_data,z_data)=plot_slice_mayavi(grid_name, png_name, hypo_x, hypo_y, hypo_z, options.search_grid,max_stack_value)
  png=Image.open(png_name)

  stack=read(stack_file, starttime=start_time, endtime=end_time)
  stack_tr=stack.traces[0]

  n_traces=len(data_files)+1

  print "Creating figure..."

  # create figure

  fig = plt.figure(figsize=(12,6))
  fig.subplots_adjust(top=0.90, bottom=0.10, left=0.1, right=0.90)

  # plot data and stack
  ax=fig.add_subplot(n_traces,4,1,title="Data")
  ax.set_axis_off()
  for ifile in range(len(data_files)):
    ax=fig.add_subplot(n_traces,4,4*ifile+1)
    ax.set_axis_off()
    try:
      st=read(data_files[ifile],starttime=start_time, endtime=end_time)
      tr=st.traces[0]
      ax.plot(tr.data)
      pos=list(ax.get_position().bounds)
      fig.text(pos[0] - 0.01, pos[1], tr.stats.station, fontsize=10, horizontalalignment='right')
      # get and plot p_time
      try:
        p_time=phases[tr.stats.station]
        if p_time > start_time and p_time < end_time:
          i,j=ax.get_ylim()
          ax.vlines(np.int(np.floor((p_time-start_time)/tr.stats.delta)), i, j, color='r', lw=1)
      except KeyError:
        pass
      except NameError:
        pass
    except IndexError:
     pass

  ax1=fig.add_subplot(n_traces,4,4*(n_traces-1)+1)
  ax1.set_axis_off()
  ax1.plot(stack_tr.data,'r')
  i,j=ax1.get_ylim()
  ax1.vlines(np.int(np.floor((stack_time-start_time)/stack_tr.stats.delta)), i, j, color='g', lw=2)
  ax1.hlines((i+j)/2, np.int(np.floor((stack_time-start_time-stack_time_err_left)/stack_tr.stats.delta)), np.int(np.floor((stack_time-start_time+stack_time_err_right)/stack_tr.stats.delta)), color='g', lw=2)
  pos=list(ax1.get_position().bounds)
  fig.text(pos[0] - 0.01, pos[1], stack_tr.stats.station, fontsize=10, horizontalalignment='right')

  if do_hyp:
    # put time of hyp
    if otime > start_time and otime < end_time:
      ax1.vlines(np.int(np.floor((otime-start_time)/stack_tr.stats.delta)), i, j, color='b', lw=2)

  
  # plot kurtosis and stack
  ax=fig.add_subplot(n_traces,4,2,title="Kurtosis")
  ax.set_axis_off()
  for ifile in range(len(kurt_files)):
    ax=fig.add_subplot(n_traces,4,4*ifile+2)
    ax.set_axis_off()
    try:
      st=read(kurt_files[ifile],starttime=start_time, endtime=end_time)
      tr=st.traces[0]
      ax.plot(tr.data)
      # get and plot p_time
      try:
        p_time=phases[tr.stats.station]
        if p_time > start_time and p_time < end_time:
          i,j=ax.get_ylim()
          ax.vlines(np.int(np.floor((p_time-start_time)/tr.stats.delta)), i, j, color='r', lw=1)
      except KeyError:
        pass
      except NameError:
        pass
    except IndexError:
     pass
  ax1=fig.add_subplot(n_traces,4,4*(n_traces-1)+2)
  ax1.set_axis_off()
  ax1.plot(stack_tr.data,'r')
  i,j=ax1.get_ylim()
  ax1.vlines(np.int(np.floor((stack_time-start_time)/stack_tr.stats.delta)), i, j, color='g', lw=2)
  ax1.hlines((i+j)/2, np.int(np.floor((stack_time-start_time-stack_time_err_left)/stack_tr.stats.delta)), np.int(np.floor((stack_time-start_time+stack_time_err_right)/stack_tr.stats.delta)), color='g', lw=2)
  if do_hyp:
    # put time of hyp
    if otime > start_time and otime < end_time:
      ax1.vlines(np.int(np.floor((otime-start_time)/stack_tr.stats.delta)), i, j, color='b', lw=2)
 
  # plot isometkurtric mayavi plot through location
  ax2a=fig.add_subplot(2,4,3, title = "Isometric")
  ax2a.set_axis_off()
  im = plt.imshow(png,origin='lower')

  # plot isometric mayavi plot through location
  ax2b=fig.add_subplot(2,4,4, title = "YZ plane")
  x_data.shape=(dummy_grid.nz,dummy_grid.ny)
  plt.imshow(x_data,cmap=plt.cm.jet,vmin=-1,vmax=max_stack_value)
  plt.plot((stack_y-dummy_grid.y_orig)/dummy_grid.dy, (stack_z-dummy_grid.z_orig)/dummy_grid.dz, 'bo')
  try:
    plt.plot((hypo_y-dummy_grid.y_orig)/dummy_grid.dy, (hypo_z-dummy_grid.z_orig)/dummy_grid.dz, 'ro')
  except TypeError:
    pass


  # plot isometric mayavi plot through location
  ax2c=fig.add_subplot(2,4,7, title = "XZ plane")
  y_data.shape=(dummy_grid.nz,dummy_grid.nx)
  plt.imshow(y_data,cmap=plt.cm.jet,vmin=-1,vmax=max_stack_value)
  plt.plot((stack_x-dummy_grid.x_orig)/dummy_grid.dx, (stack_z-dummy_grid.z_orig)/dummy_grid.dz, 'bo')
  try:
    plt.plot((hypo_x-dummy_grid.x_orig)/dummy_grid.dx, (hypo_z-dummy_grid.z_orig)/dummy_grid.dz, 'ro')
  except TypeError:
    pass

  # plot isometric mayavi plot through location
  ax2d=fig.add_subplot(2,4,8, title = "XY plane")
  z_data.shape=(dummy_grid.ny,dummy_grid.nx)
  plt.imshow(z_data,cmap=plt.cm.jet,vmin=-1,vmax=max_stack_value)
  plt.plot((stack_x-dummy_grid.x_orig)/dummy_grid.dx, (stack_y-dummy_grid.y_orig)/dummy_grid.dy, 'bo')
  try:
    plt.plot((hypo_x-dummy_grid.x_orig)/dummy_grid.dx, (hypo_y-dummy_grid.y_orig)/dummy_grid.dy, 'ro')
  except TypeError:
    pass


  plt.savefig(plot_filename)
  logging.info("Saved figure in %s"%plot_filename)
  plt.clf()
