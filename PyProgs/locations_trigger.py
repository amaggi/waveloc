#!/usr/bin/env python
# encoding: utf-8

import os, sys, optparse, glob, h5py
from obspy.core import *
from obspy.signal import *
from OP_waveforms import stream_taper, Waveform
from filters import smooth
import matplotlib.pyplot as plt
import numpy as np
import logging
from OP_waveforms import read_data_compatible_with_time_dict
from hdf5_grids import get_interpolated_time_grids

def plot_location_triggers(trace,trig_start,trig_end,trig_95_start, trig_95_end, show=True):
    df = trace.stats.sampling_rate
    dt = trace.stats.delta
    npts = trace.stats.npts
    t = np.arange(npts) / df
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax1.plot(t[trig_start-100:trig_end+100], trace.data[trig_start-100:trig_end+100], 'k')
    i, j = ax1.get_ylim()
    try:
        ax1.vlines(trig_start / df , i, j, color='r', lw=1, label="Trigger On")
        ax1.vlines(trig_end / df, i, j, color='b', lw=1, label="Trigger Off")
        ax1.vlines(trig_95_start / df , i, j, color='r', lw=2 )
        ax1.vlines(trig_95_end / df , i, j, color='b', lw=2 )
#        ax1.legend()
    except IndexError:
        pass
    fig.suptitle(trace.id)
    fig.canvas.draw()
    if show:
        plt.show()

def filter_max_stack(st_max,corner):
# maximum of stack is too rough - filter it 
  logging.debug("Doing filtering")

  st_filt = st_max.copy()
  st_filt.clear()

  for tr in st_max.traces:
    tr_filt=tr.copy()
    #x_filt=lowpass(tr_filt.data,corner,1/tr.stats.delta,zerophase=True)
    x_filt=smooth(tr_filt.data)
    tr_filt.data=x_filt
    st_filt.append(tr_filt)


  logging.debug("Done!")

  return st_filt

def number_good_kurtosis_for_location(kurt_files,data_files,loc,time_dict,snr_limit=10.0,snr_tr_limit=10.0,sn_time=10.0):

  o_time=loc['o_time']
  stack_x=loc['x_mean']
  stack_y=loc['y_mean']
  stack_z=loc['z_mean']
  # TODO - Fix this to estimate K-time from o_time and propagation time to station
  n_good_kurt=0
  wf=Waveform()

  for ifile in xrange(len(kurt_files)):
    kfilename=kurt_files[ifile]
    dfilename=data_files[ifile]

    st=read(kfilename,headonly=True)
    staname=st.traces[0].stats.station

    traveltime=time_dict[staname].value_at_point(stack_x,stack_y,stack_z)
    start_time=o_time+traveltime-sn_time
    end_time=o_time+traveltime+sn_time
    try:
      wf.read_from_file(kfilename,starttime=start_time,endtime=end_time)
      snr=wf.get_snr(o_time+traveltime,start_time,end_time)

      wf.read_from_file(dfilename,starttime=start_time,endtime=end_time)
      snr_tr=wf.get_snr(o_time+traveltime,start_time,end_time)

      if snr > snr_limit and snr_tr > snr_tr_limit:
        n_good_kurt = n_good_kurt + 1
    except UserWarning:
      logging.info('No data around %s for file %s.'%(o_time.isoformat(),filename))
  return n_good_kurt

def trigger_locations_inner(max_val,max_x,max_y,max_z,left_trig,right_trig,start_time,delta):    

    locs=[]
    trigs=trigger.triggerOnset(np.array(max_val),left_trig,right_trig)

    df=1/delta

    logging.debug('Found %d triggers.'%len(trigs))

    for trig in trigs:

      i_start=trig[0]
      i_end=trig[1]+1
      i_max_trig=np.argmax(max_val[i_start:i_end])+i_start
      max_trig=max_val[i_max_trig]
      max_trig_95=0.95*max_trig
      logging.debug('Max_trig = %.3f, max_trig_95 = %.3f'%(max_trig,max_trig_95))
      trigs_95=trigger.triggerOnset(max_val[i_start:i_end],max_trig_95,max_trig_95)
      for trig_95 in trigs_95:
        if i_max_trig >= trig_95[0]+i_start and i_max_trig <= trig_95[1]+i_start:
          loc_dict={}
          loc_dict['max_trig']=max_trig
          i_start_95=trig_95[0]+i_start
          i_end_95=trig_95[1]+1+i_start
          loc_dict['x_mean']=np.mean(max_x[i_start_95:i_end_95])
          loc_dict['x_sigma']=np.std(max_x[i_start_95:i_end_95])
          loc_dict['y_mean']=np.mean(max_y[i_start_95:i_end_95])
          loc_dict['y_sigma']=np.std(max_y[i_start_95:i_end_95])
          loc_dict['z_mean']=np.mean(max_z[i_start_95:i_end_95])
          loc_dict['z_sigma']=np.std(max_z[i_start_95:i_end_95])

          loc_dict['o_time']=start_time + i_max_trig*delta
          loc_dict['o_err_left']=(i_max_trig-i_start_95)*delta
          loc_dict['o_err_right']=(i_end_95-i_max_trig)*delta

          #locs.append([max_trig,o_time,o_err_left, o_err_right,x_mean,x_sigma,y_mean,y_sigma,z_mean,z_sigma])
          locs.append(loc_dict)

    return locs


def trigger_locations(st_max_filt,st_x,st_y,st_z,left_trig,right_trig):

  logging.warn('Deprecated - do not use')

  locs=[]
  for i_filt in range(st_max_filt.count()):

    tr_filt=st_max_filt.traces[i_filt]
    tr_x=st_x.traces[i_filt]
    tr_y=st_y.traces[i_filt]
    tr_z=st_z.traces[i_filt]

    i_locs=trigger_locations_inner(tr_filt.data,tr_x.data,tr_y.data,tr_z.data,left_trig,right_trig,tr_filt.stats.delta)
    for loc in i_locs:
      # fix-up the origin time
      loc['o_time'] = tr_filt.stats.starttime+loc['o_time']

    locs.extend(i_locs)

  return locs
   
def do_locations_trigger_setup_and_run(opdict):

  base_path=opdict['base_path']
  # parse command line
  data_dir=os.path.join(base_path,'data',opdict['datadir'])
  kurt_files=glob.glob(os.path.join(data_dir,opdict['kurtglob']))
  data_files=glob.glob(os.path.join(data_dir,opdict['dataglob']))
  kurt_files.sort()
  data_files.sort()

  dataglob=opdict['dataglob']

  time_grids=get_interpolated_time_grids(opdict)

  logging.info("Starting log for combine_stacks.")

  out_path=os.path.join(base_path,'out',opdict['outdir'])
  stack_path=os.path.join(out_path,'stack')

  reloc=opdict['reloc']
  if reloc:
    loc_path=os.path.join(out_path,'reloc')
    stack_files=glob.glob(os.path.join(stack_path,'reloc_stack_all*.hdf5'))
    stack_files.sort()
  else:
    loc_path=os.path.join(out_path,'loc')
    stack_files=glob.glob(os.path.join(stack_path,'stack_all*.hdf5'))
    stack_files.sort()

  n_stacks=len(stack_files)
  if n_stacks == 0 :
    raise UserWarning('Empty list of stacks in %s'%(stack_path))

  loc_filename=os.path.join(loc_path,"locations.dat")
  logging.info("Path for stack files : %s"%stack_path)
  logging.info("Path for loc files : %s"%loc_path)
  logging.info("Location file : %s"%loc_filename)

  # DO DATA PREP ACCORDING TO RELOC OR NOT

  logging.info("\nDealing with continuous location, so merging stack files directly ...\n")

  # get basic info from first file
  f_stack = h5py.File(stack_files[0],'r')
  max_val = f_stack['max_val']
  dt = max_val.attrs['dt']
  f_stack.close()

  # get start times  (get first and last times)
  start_times=[]
  end_times=[]
  for fname in stack_files:
    f_stack = h5py.File(fname,'r')
    max_val = f_stack['max_val']
    start_times.append(utcdatetime.UTCDateTime(max_val.attrs['start_time']))
    end_times.append(  utcdatetime.UTCDateTime(max_val.attrs['start_time'])+dt*len(max_val))
    f_stack.close()

  first_start_time = min(start_times)
  last_end_time = max(end_times)

  nt_full=int((last_end_time-first_start_time)/dt)+1


  # create - assume all stacks are of the same length and will be concatenated end to end 
  #          (this will give more than enough space) 
  f = h5py.File(os.path.join(stack_path,'combined_stack_all.hdf5'),'w')
  cmax_val = f.create_dataset('max_val',(nt_full,), 'f', chunks=(nt_full,))
  cmax_x = f.create_dataset('max_x',(nt_full,), 'f', chunks=(nt_full,))
  cmax_y = f.create_dataset('max_y',(nt_full,), 'f', chunks=(nt_full,))
  cmax_z = f.create_dataset('max_z',(nt_full,), 'f', chunks=(nt_full,))

  # concatenate unsmoothed versions of max_val to avoid 
  # problems at file starts and ends
  for i in range(n_stacks):
    f_stack = h5py.File(stack_files[i],'r')
    max_val = f_stack['max_val']
    max_x = f_stack['max_x']
    max_y = f_stack['max_y']
    max_z = f_stack['max_z']

    # get time info for this stack
    nt = len(max_val)
    start_time = utcdatetime.UTCDateTime(max_val.attrs['start_time'])
    ibegin=np.int((start_time-first_start_time)/dt)

    # copy data over into the right place
    cmax_val[ibegin:ibegin+nt] = max_val[:]
    cmax_x[ibegin:ibegin+nt] = max_x[:]
    cmax_y[ibegin:ibegin+nt] = max_y[:]
    cmax_z[ibegin:ibegin+nt] = max_z[:]

    # close the stack
    f_stack.close()

  # create the smoothed version of the max stack
  cmax_val_smooth = f.create_dataset('max_val_smooth',(nt_full,), 'f', chunks=(nt_full,))
  cmax_val_smooth[:] = smooth(np.array(cmax_val),51)

  for name in f:
    dset=f[name]
    dset.attrs['dt']=dt
    dset.attrs['start_time']=first_start_time.isoformat()


  # DO TRIGGERING AND LOCATION
  if opdict['auto_loclevel']:
    loclevel=opdict['snr_loclevel']*np.median(cmax_val_smooth)
    opdict['loclevel']=loclevel
  else:
    loclevel=opdict['loclevel']
  left_trig=loclevel
  right_trig=loclevel

  loc_list=trigger_locations_inner(cmax_val_smooth[:],cmax_x,cmax_y,cmax_z,left_trig,right_trig,first_start_time,dt)
  logging.info('Found %d initial.'%(len(loc_list)))

  # close the stack file
  f.close()

  loc_file=open(loc_filename,'w')

  snr_limit=opdict['snr_limit']
  snr_tr_limit=opdict['snr_tr_limit']
  sn_time=opdict['sn_time']
  n_kurt_min=opdict['n_kurt_min']

  # Header of locations.dat
  loc_file.write('FILTER : %.d - %.d Hz\n'%(opdict['c1'],opdict['c2']))
  loc_file.write('KURTOSIS = window: %.2f s, recurs: %s, grad: %s, gauss: %s\n'%(opdict['kwin'],opdict['krec'],opdict['kderiv'],opdict['gauss']))
  loc_file.write('OPTIONS = reloc: %s\n'%reloc)
  loc_file.write('LOCATION = level: %d, window of analysis: %.2f s, kurtosis snr: %.2f, waveform snr: %.2f, number of stations: %d\n\n'%(loclevel,sn_time,snr_limit,snr_tr_limit,n_kurt_min))

  n_ok=0
  locs=[]
  for loc in loc_list:
    if number_good_kurtosis_for_location(kurt_files,data_files,loc,time_grids,snr_limit,snr_tr_limit,sn_time) > n_kurt_min:
      logging.info("Max = %.2f, %s - %.2fs + %.2f s, x=%.4f pm %.4f km, y=%.4f pm %.4f km, z=%.4f pm %.4f km"%(loc['max_trig'],loc['o_time'].isoformat(),loc['o_err_left'], loc['o_err_right'],loc['x_mean'],loc['x_sigma'],loc['y_mean'],loc['y_sigma'],loc['z_mean'],loc['z_sigma']))
      loc_file.write("Max = %.2f, %s - %.2f s + %.2f s, x= %.4f pm %.4f km, y= %.4f pm %.4f km, z= %.4f pm %.4f km\n"%(loc['max_trig'],loc['o_time'].isoformat(),loc['o_err_left'], loc['o_err_right'],loc['x_mean'],loc['x_sigma'],loc['y_mean'],loc['y_sigma'],loc['z_mean'],loc['z_sigma']))
      n_ok=n_ok+1
      locs.append(loc)
    else:
      logging.info("Not enough kurtosis picks for : Max = %.2f, %s - %.2fs + %.2fs, x=%.4f pm %.4f, y=%.4f pm %.4f, z=%.4f pm %.4f"%(loc['max_trig'],loc['o_time'].isoformat(),loc['o_err_left'], loc['o_err_right'],loc['x_mean'],loc['x_sigma'],loc['y_mean'],loc['y_sigma'],loc['z_mean'],loc['z_sigma']))
  loc_file.close()
  logging.info('Wrote %d locations to file %s.'%(n_ok,loc_filename))

  return locs

def read_locs_from_file(filename):

  from obspy.core import utcdatetime

  locs=[]

  f=open(filename,'r')
  lines=f.readlines()
  f.close()

  nb_line=0

  for line in lines:

    if nb_line >= 5:
      loc={}

      loc['max_trig']=np.float(line.split()[2].split(',')[0])
      loc['o_time']=utcdatetime.UTCDateTime(line.split()[3])
      loc['o_err_left']=np.float(line.split()[5])
      loc['o_err_right']=np.float(line.split()[8])
      loc['x_mean']=np.float(line.split()[11])
      loc['x_sigma']=np.float(line.split()[13])
      loc['y_mean']=np.float(line.split()[16])
      loc['y_sigma']=np.float(line.split()[18])
      loc['z_mean']=np.float(line.split()[21])
      loc['z_sigma']=np.float(line.split()[23])

      locs.append(loc)

    nb_line=nb_line+1

  return locs


def read_header_from_file(filename):

  f=open(filename,'r')
  lines=f.readlines()
  f.close()

  head={}

  head['filt']='%s-%s'%(lines[0].split()[2],lines[0].split()[4])
  head['rec']=lines[1].split()[6][:-1]
  head['grad']=lines[1].split()[8][:-1]
  head['gauss']=lines[1].split()[10]
  head['reloc']=lines[2].split()[3]
  head['detect']=lines[3].split()[3][:-1]

  return head


 
if __name__=='__main__':

  from options import WavelocOptions
  logging.basicConfig(level=logging.INFO, format='%(levelname)s : %(asctime)s : %(message)s')

  wo = WavelocOptions()
  args=wo.p.parse_args()

  wo.set_all_arguments(args)
  wo.verify_location_options()

  do_locations_trigger_setup_and_run(wo.opdict)
