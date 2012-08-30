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
      snr=wf.get_snr(o_time,start_time,end_time)
      if snr > snr_limit:
        n_good_kurt = n_good_kurt + 1
    except UserWarning:
      logging.info('No data around %s for file %s.'%(o_time.isoformat(),filename))
  return n_good_kurt
    
def trigger_locations_inner(max_val,max_x,max_y,max_z,left_trig,right_trig,delta):    

    locs=[]
    trigs=trigger.triggerOnset(max_val,left_trig,right_trig)

    df=1/delta

    logging.debug('Found %d triggers.'%len(trigs))

    for trig in trigs:

      i_start=trig[0]
      i_end=trig[1]+1
      i_max_trig=np.argmax(max_val.data[i_start:i_end])+i_start
      max_trig=max_val[i_max_trig]
      max_trig_95=0.95*max_trig
      logging.debug('Max_trig = %.3f, max_trig_95 = %.3f'%(max_trig,max_trig_95))
      trigs_95=trigger.triggerOnset(max_val[i_start:i_end],max_trig_95,max_trig_95)
      for trig_95 in trigs_95:
        if i_max_trig >= trig_95[0]+i_start and i_max_trig <= trig_95[1]+i_start:
          i_start_95=trig_95[0]+i_start
          i_end_95=trig_95[1]+1+i_start
          x_mean=np.mean(max_x[i_start_95:i_end_95])
          x_sigma=np.std(max_x[i_start_95:i_end_95])
          y_mean=np.mean(max_y[i_start_95:i_end_95])
          y_sigma=np.std(max_y[i_start_95:i_end_95])
          z_mean=np.mean(max_z[i_start_95:i_end_95])
          z_sigma=np.std(max_z[i_start_95:i_end_95])

          o_time=i_max_trig*delta
          o_err_left=(i_max_trig-i_start_95)*delta
          o_err_right=(i_end_95-i_max_trig)*delta

          locs.append([max_trig,o_time,o_err_left, o_err_right,x_mean,x_sigma,y_mean,y_sigma,z_mean,z_sigma])

    return locs
 

def trigger_locations(st_max_filt,st_x,st_y,st_z,left_trig,right_trig):

  locs=[]
  for i_filt in range(st_max_filt.count()):

    tr_filt=st_max_filt.traces[i_filt]
    tr_x=st_x.traces[i_filt]
    tr_y=st_y.traces[i_filt]
    tr_z=st_z.traces[i_filt]

    i_locs=trigger_locations_inner(tr_filt.data,tr_x.data,tr_y.data,tr_z.data,left_trig,right_trig,tr_filt.stats.delta)
    for loc in i_locs:
      # fix-up the origin time
      loc[1] = tr_filt.stats.starttime+loc[1]

    locs.extend(i_locs)

  return locs
   
def do_locations_trigger_setup_and_run(opdict):

  base_path=opdict['base_path']
  # parse command line
  data_dir=os.path.join(base_path,'data',opdict['datadir'])
  kurt_files=glob.glob(os.path.join(data_dir,opdict['kurtglob']))

  # corner frequency for lowpass filtering of max stack
  corner=1.0 
  loclevel=opdict['loclevel']
  left_trig=loclevel
  right_trig=loclevel

  # start logging
  #logfile=base_path + os.sep + 'out'+ os.sep +  outdir + os.sep + 'combine_stacks.log'

  logging.info("Starting log for combine_stacks.")

  out_path=os.path.join(base_path,'out',opdict['outdir'])
  stack_path=os.path.join(out_path,'stack')

  reloc=opdict['reloc']
  if reloc:
    loc_path=os.path.join(out_path,'reloc')
  else:
    loc_path=os.path.join(out_path,'loc')


  loc_filename=os.path.join(loc_path,"locations.dat")
  logging.info("Path for stack files : %s"%stack_path)
  logging.info("Path for loc files : %s"%loc_path)
  logging.info("Location file : %s"%loc_filename)

  # DO DATA PREP ACCORDING TO RELOC OR NOT

  if reloc:
    logging.info("\nDealing with relocation, so taper before merging ...\n")

    st_max=read(os.path.join(stack_path,"reloc_stack_max*"))
    st_max=stream_taper(st_max)
    st_max.merge(method=1,fill_value='interpolate')
    st_max.write(os.path.join(stack_path,"combined_reloc_stack_max.mseed"),format='MSEED')
    st_max_filt=filter_max_stack(st_max,corner)
    st_max_filt.write(os.path.join(stack_path,"combined_reloc_stack_max_filt.mseed"),format='MSEED')

    st_x=read(os.path.join(stack_path,"reloc_stack_x*"))
    st_x=stream_taper(st_x)
    st_x.merge(method=1,fill_value='interpolate')
    st_x.write(os.path.join(stack_path,"combined_reloc_stack_x.mseed"),format='MSEED')

    st_y=read(os.path.join(stack_path,"reloc_stack_y*"))
    st_y=stream_taper(st_y)
    st_y.merge(method=1,fill_value='interpolate')
    st_y.write(os.path.join(stack_path,"combined_reloc_stack_y.mseed"),format='MSEED')

    st_z=read(os.path.join(stack_path,"reloc_stack_z*"))
    st_z=stream_taper(st_z)
    st_z.merge(method=1,fill_value='interpolate')
    st_z.write(os.path.join(stack_path,"combined_reloc_stack_z.mseed"),format='MSEED')


  else:

    logging.info("\nDealing with continuous location, so merging stack files directly ...\n")

    st_max=read(os.path.join(stack_path,"stack_max*"))
    st_max.merge(method=1,fill_value='interpolate')
    st_max.write(os.path.join(stack_path,"combined_stack_max.mseed"),format='MSEED')
    st_max_filt=filter_max_stack(st_max,corner)
    st_max_filt.write(os.path.join(stack_path,"combined_stack_max_filt.mseed"),format='MSEED')

    st_x=read(os.path.join(stack_path,"stack_x*"))
    st_x.merge(method=1,fill_value='interpolate')
    st_x.write(os.path.join(stack_path,"combined_stack_x.mseed"),format='MSEED')

    st_y=read(os.path.join(stack_path,"stack_y*"))
    st_y.merge(method=1,fill_value='interpolate')
    st_y.write(os.path.join(stack_path,"combined_stack_y.mseed"),format='MSEED')

    st_z=read(os.path.join(stack_path,"stack_z*"))
    st_z.merge(method=1,fill_value='interpolate')
    st_z.write(os.path.join(stack_path,"combined_stack_z.mseed"),format='MSEED')



  # DO TRIGGERING AND LOCATION

  loc_list=trigger_locations(st_max_filt,st_x,st_y,st_z,left_trig,right_trig)
  logging.info('Found %d initial.'%(len(loc_list)))

  loc_file=open(loc_filename,'w')

  snr_limit=opdict['snr_limit']
  sn_time=opdict['sn_time']
  n_kurt_min=opdict['n_kurt_min']

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

 
if __name__=='__main__':

  from options import WavelocOptions
  logging.basicConfig(level=logging.INFO, format='%(levelname)s : %(asctime)s : %(message)s')

  wo = WavelocOptions()
  args=wo.p.parse_args()

  wo.set_all_arguments(args)
  wo.verify_location_options()

  do_locations_setup_and_run(wo.opdict)

 
