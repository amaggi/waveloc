#!/usr/bin/env python
# encoding: utf-8

import os, sys, optparse, glob, logging
from obspy.core import *
from obspy.signal import *
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
p.add_option('--low_pass',action='store',default=1.0, type='float',help="low-pass filter to apply to max_stack for robust primary detection")

(options,arguements)=p.parse_args()

logging.info("Starting log for locations_prob.")

stack_path=os.path.join(base_path,'out',options.outdir,'stack')
loc_path=os.path.join(base_path,'out',options.outdir,'loc')
loc_filename=os.path.join(loc_path,'locations_prob.dat')

logging.info("Path for stack files : %s"%stack_path)
logging.info("Path for loc files : %s"%loc_path)
logging.info("Location file : %s"%loc_filename)

logging.info("Merging stack files and low-pass filtering at %.2fs..."%options.low_pass)
st_max=read(os.path.join(stack_path,"stack_max*"))
st_max.merge(method=1,fill_value='interpolate')
st_max.write(os.path.join(stack_path,"combined_stack_max.mseed"),format='MSEED')
st_max.write(os.path.join(stack_path,"combined_stack_max.sac"),format='SAC')
st_max_filt=filter_max_stack(st_max,options.low_pass)
st_max_filt.write(os.path.join(stack_path,"combined_stack_max_filt.mseed"),format='MSEED')
st_max_filt.write(os.path.join(stack_path,"combined_stack_max_filt.sac"),format='SAC')

detection_list=trigger_detections(st_max_filt,options.loclevel)

print detection_list

