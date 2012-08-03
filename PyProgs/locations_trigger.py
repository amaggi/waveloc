#!/usr/bin/env python
# encoding: utf-8

import os, sys, optparse, glob
from obspy.core import *
from obspy.signal import *
from OP_waveforms import stream_taper, Waveform
import matplotlib.pyplot as plt
import numpy as np
import logging

def plot_location_triggers(trace,trig_start,trig_end,trig_95_start, trig_95_end, show=True):
    df = trace.stats.sampling_rate
    dt = trace.stats.delta
    npts = trace.stats.npts
    t = np.arange(npts, dtype='float32') / df
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
    x_filt=lowpass(tr_filt.data,corner,1/tr.stats.delta,zerophase=True)
    tr_filt.data=x_filt
    st_filt.append(tr_filt)


  logging.debug("Done!")

  return st_filt

def number_good_kurtosis_for_location(kurt_files,o_time,snr_limit=10.0,sn_time=10.0):
  n_good_kurt=0
  wf=Waveform()
  start_time=o_time-sn_time
  end_time=o_time+sn_time 
  for filename in kurt_files:
    try:
      wf.read_from_file(filename,starttime=start_time,endtime=end_time)
      k_max=np.max(wf.trace.data)
      wf.read_from_file(filename,starttime=start_time,endtime=o_time)
      k_mean=np.mean(np.abs(wf.trace.data))
      if k_mean==0.0:
        snr=0
      else:
        snr=k_max/k_mean
      if snr > snr_limit:
        n_good_kurt = n_good_kurt + 1
    except UserWarning:
      logging.info('No data around %s for file %s.'%(o_time.isoformat(),filename))
  return n_good_kurt
    
    

def trigger_locations(st_max_filt,st_x,st_y,st_z,left_trig,right_trig):

  locs=[]
  for i_filt in range(st_max_filt.count()):

    tr_filt=st_max_filt.traces[i_filt]
    tr_x=st_x.traces[i_filt]
    tr_y=st_y.traces[i_filt]
    tr_z=st_z.traces[i_filt]

    trigs=trigger.triggerOnset(tr_filt.data,left_trig,right_trig)

    df=1/tr_filt.stats.delta

    logging.debug('Found %d triggers.'%len(trigs))

    for trig in trigs:

      i_start=trig[0]
      i_end=trig[1]+1
      i_max_trig=np.argmax(tr_filt.data[i_start:i_end])+i_start
      max_trig=tr_filt.data[i_max_trig]
      max_trig_95=0.95*max_trig
      trigs_95=trigger.triggerOnset(tr_filt.data[i_start:i_end],max_trig_95,max_trig_95)
      for trig_95 in trigs_95:
        if i_max_trig > trig_95[0]+i_start and i_max_trig < trig_95[1]+i_start:
          i_start_95=trig_95[0]+i_start
          i_end_95=trig_95[1]+1+i_start
          x_mean=np.mean(tr_x[i_start_95:i_end_95])
          x_sigma=np.sqrt(np.var(tr_x[i_start_95:i_end_95]))
          y_mean=np.mean(tr_y[i_start_95:i_end_95])
          y_sigma=np.sqrt(np.var(tr_y[i_start_95:i_end_95]))
          z_mean=np.mean(tr_z[i_start_95:i_end_95])
          z_sigma=np.sqrt(np.var(tr_z[i_start_95:i_end_95]))

          o_time=tr_filt.stats.starttime+i_max_trig*tr_filt.stats.delta
          o_err_left=(i_max_trig-i_start_95)*tr_filt.stats.delta
          o_err_right=(i_end_95-i_max_trig)*tr_filt.stats.delta

          locs.append((max_trig,o_time,o_err_left, o_err_right,x_mean,x_sigma,y_mean,y_sigma,z_mean,z_sigma))

  return locs
   
#  plot_location_trigg

if __name__=='__main__':

  # get path
  base_path=os.getenv('WAVELOC_PATH_PDF')

  # Read command line

  p = optparse.OptionParser()
  p.add_option('--outdir', '-o', action='store', help='output subdirectory in which the stack directory is found')
  p.add_option('--reloc', action='store_true', default=False, help='apply to relocated events')
  p.add_option('--loclevel', action='store', default=100, help='trigger stack level for locations (e.g. 100) ')
  p.add_option('--datadir',action='store',help="data subdirectory")
  p.add_option('--dataglob',action='store',help="data glob")
  p.add_option('--snr_limit',action='store',default=10.0, help="signal_to_noise level for kurtosis acceptance")
  p.add_option('--sn_time',action='store',default=10.0, help="time over which to calculate the signal_to_noise ratio for kurtosis acceptance")
  p.add_option('--n_kurt_min',action='store',default=4, help="min number of good kurtosis traces for a location")

  (options,arguements)=p.parse_args()

  # parse command line
  out_dir=options.outdir
  do_reloc=options.reloc
  loc_level=np.float(options.loclevel)
  data_dir=os.path.join(base_path,'data',options.datadir)
  kurt_files=glob.glob(os.path.join(data_dir,options.dataglob))
  snr_limit=np.float(options.snr_limit)
  sn_time=np.float(options.sn_time)
  n_kurt_min=np.int(options.n_kurt_min)

  # corner frequency for lowpass filtering of max stack
  corner=1.0 
  left_trig=loc_level
  right_trig=loc_level

  # start logging
  #logfile=base_path + os.sep + 'out'+ os.sep +  out_dir + os.sep + 'combine_stacks.log'
  logging.basicConfig(level=logging.INFO, format='%(levelname)s : %(asctime)s : %(message)s')

  logging.info("Starting log for combine_stacks.")

  stack_path="%s/out/%s/stack"%(base_path,out_dir)

  if do_reloc:
    loc_path="%s/out/%s/reloc"%(base_path,out_dir)
  else:
    loc_path="%s/out/%s/loc"%(base_path,out_dir)

  loc_filename="%s/locations.dat"%(loc_path)
  logging.info("Path for stack files : %s"%stack_path)
  logging.info("Path for loc files : %s"%loc_path)
  logging.info("Location file : %s"%loc_filename)

  # DO DATA PREP ACCORDING TO RELOC OR NOT

  if do_reloc:
    logging.info("\nDealing with relocation, so taper before merging ...\n")

    st_max=read("%s/reloc_stack_max*"%stack_path)
    st_max=stream_taper(st_max)
    st_max.merge(method=1,fill_value='interpolate')
    st_max.write("%s/combined_reloc_stack_max.mseed"%stack_path,format='MSEED')
    st_max_filt=filter_max_stack(st_max,corner)
    st_max_filt.write("%s/combined_reloc_stack_max_filt.mseed"%stack_path,format='MSEED')

    st_x=read("%s/reloc_stack_x*"%stack_path)
    st_x=stream_taper(st_x)
    st_x.merge(method=1,fill_value='interpolate')
    st_x.write("%s/combined_reloc_stack_x.mseed"%stack_path,format='MSEED')

    st_y=read("%s/reloc_stack_y*"%stack_path)
    st_y=stream_taper(st_y)
    st_y.merge(method=1,fill_value='interpolate')
    st_y.write("%s/combined_reloc_stack_y.mseed"%stack_path,format='MSEED')

    st_z=read("%s/reloc_stack_z*"%stack_path)
    st_z=stream_taper(st_z)
    st_z.merge(method=1,fill_value='interpolate')
    st_z.write("%s/combined_reloc_stack_z.mseed"%stack_path,format='MSEED')


  else:

    logging.info("\nDealing with continuous location, so merging stack files directly ...\n")

    st_max=read("%s/stack_max*"%stack_path)
    st_max.merge(method=1,fill_value='interpolate')
    st_max.write("%s/combined_stack_max.mseed"%stack_path,format='MSEED')
    st_max_filt=filter_max_stack(st_max,corner)
    st_max_filt.write("%s/combined_stack_max_filt.mseed"%stack_path,format='MSEED')

    st_x=read("%s/stack_x*"%stack_path)
    st_x.merge(method=1,fill_value='interpolate')
    st_x.write("%s/combined_stack_x.mseed"%stack_path,format='MSEED')

    st_y=read("%s/stack_y*"%stack_path)
    st_y.merge(method=1,fill_value='interpolate')
    st_y.write("%s/combined_stack_y.mseed"%stack_path,format='MSEED')

    st_z=read("%s/stack_z*"%stack_path)
    st_z.merge(method=1,fill_value='interpolate')
    st_z.write("%s/combined_stack_z.mseed"%stack_path,format='MSEED')



  # DO TRIGGERING AND LOCATION

  loc_list=trigger_locations(st_max_filt,st_x,st_y,st_z,left_trig,right_trig)
  logging.info('Found %d initial.'%(len(loc_list)))

  loc_file=open(loc_filename,'w')

  n_ok=0
  for (max_trig,o_time,o_err_left, o_err_right,x_mean,x_sigma,y_mean,y_sigma,z_mean,z_sigma) in loc_list:
    if number_good_kurtosis_for_location(kurt_files,o_time,snr_limit,sn_time) > n_kurt_min:
      logging.info("Max = %.2f, %s - %.2fs + %.2fs, x=%.4f pm %.4f, y=%.4f pm %.4f, z=%.4f pm %.4f"%(max_trig,o_time.isoformat(),o_err_left, o_err_right,x_mean,x_sigma,y_mean,y_sigma,z_mean,z_sigma))
      loc_file.write("Max = %.2f, %s - %.2f s + %.2f s, x= %.4f pm %.4f km, y= %.4f pm %.4f km, z= %.4f pm %.4f km\n"%(max_trig,o_time.isoformat(),o_err_left, o_err_right ,x_mean,x_sigma,y_mean,y_sigma,z_mean,z_sigma))
      n_ok=n_ok+1
    else:
      logging.info("Not enough kurtosis picks for : Max = %.2f, %s - %.2fs + %.2fs, x=%.4f pm %.4f, y=%.4f pm %.4f, z=%.4f pm %.4f"%(max_trig,o_time.isoformat(),o_err_left, o_err_right,x_mean,x_sigma,y_mean,y_sigma,z_mean,z_sigma))
  loc_file.close()
  logging.info('Wrote %d locations to file %s.'%(n_ok,loc_filename))

  
