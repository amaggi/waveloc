#!/usr/bin/env python
# encoding: utf-8

import os, sys, optparse, glob
from grids_paths import *
from OP_waveforms import *
from obspy.core import *
import numpy as np
from obspy.signal import *
import matplotlib.pyplot as plt
from time import time, sleep
from scipy import weave
import logging
from sub_PdF_waveloc import do_migration_loop_reloc


# get path
base_path=os.getenv('WAVELOC_PATH')

# Read command line

time_grids=['Slow_len.100m.P']
search_grids=['grid.500m.search.hdr','grid.Taisne.search.hdr']

p=optparse.OptionParser()
p.add_option('--outdir', '-o', action='store', help='output subdirectory')
p.add_option('--datadir', '-d', action='store', help='data subdirectory')
p.add_option('--search_grid',action='store',type='choice',choices=search_grids,help="search grid %s"%(search_grids))
p.add_option('--data_glob',action='store',help="data glob")
p.add_option('--kurt_glob',action='store',help="kurtosis glob")
p.add_option('--stations','-s',action='store',default='channels_HHZ.dat',help='station list (found in $WAVELOC_PATH/lib)')
p.add_option('--time_grid',action='store',type='choice',choices=time_grids,help="time grid %s"%(time_grids))
p.add_option('--snr',action='store',default='7',help="cutoff signal-to-noise ratio for re-location (snr is on kurtosis)")

(options,arguments)=p.parse_args()

options_verbose=True
options_time=True

lib_path=base_path + os.sep + "lib"
output_dir=base_path+os.sep+'out'+os.sep+options.outdir
stack_path=output_dir+os.sep+'stack'
loc_path=output_dir+os.sep+'loc'
reloc_path=output_dir+os.sep+'reloc'
data_dir=base_path+os.sep+'data'+os.sep + options.datadir
hdr_file=lib_path+os.sep+options.search_grid

stations_filename="%s/lib/%s"%(base_path,options.stations)
search_grid_filename="%s/%s"%(lib_path,options.search_grid)
grid_filename_base="%s/%s"%(lib_path,options.time_grid)

# start logging
logfile=output_dir + os.sep + 'reloc_by_snr.log'
#logging.basicConfig(filename=logfile,filemode='w',level=logging.DEBUG, format='%(levelname)s:%(asctime)s:%(message)s')
logging.basicConfig(level=logging.DEBUG, format='%(levelname)s : %(asctime)s : %(message)s')

logging.info("Starting log for reloc_by_snr.\n")

kurt_glob="*kurt.sac"
files=glob.glob(data_dir + os.sep + kurt_glob)
snr_limit=np.float(options.snr)
logging.info("SNR cutoff = %.2f"%snr_limit)

# prepare relocation file for writing
reloc_filename=reloc_path + os.sep + 'locations.dat'
reloc_file=open(reloc_filename,'w')
logging.info("Created relocations file %s"%reloc_filename)

# read complete station file and save content in a dictionary
sta_file=open(stations_filename,'r')
sta_lines=sta_file.readlines()
sta_file.close()
sta_dict={}
for line in sta_lines:
  staname=line.split()[1]
  sta_dict[staname] = line
logging.info("Read stations file %s and found %d stations."%(stations_filename,len(sta_dict)))

#  ***** reading station file ******

sta=StationList()
sta.read_from_file(stations_filename)


cha=ChannelList()
cha.populate_from_station_list(sta,["HHZ"])


######### INTERPOLATE TRAVEL TIMES #############

# The time grid will contain as array values just the travel-times needed 
# (interpolated from the full NLL files) so we can free up the memory as soon as possible

time_grid=QDTimeGrid()
time_grid.read_NLL_hdr_file(search_grid_filename)
time_grid.populate_from_time_grids(grid_filename_base,cha,ouput_dir)


# read location file
loc_filename=loc_path + os.sep + 'locations.dat'
loc_file=open(loc_filename,'r')
loc_lines=loc_file.readlines()
loc_file.close()
logging.info("Read locations file %s and found %d locations.\n"%(loc_filename,len(loc_lines)))

logging.info("Running relocation on %d events.\n"%len(loc_lines))

# make usable list from location parameters
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


#iterate over locations
for loc in locs:
  stack_time=loc[1]
  logging.info("Working on location: %s "%stack_time.isoformat())
  start_time=stack_time-20.0
  end_time=stack_time+30.0

  snr_file=reloc_path + os.sep + "%s_SNR.txt"%(stack_time.isoformat())
  logging.info("SNR file: %s "%snr_file)
  snrf=open(snr_file,'w')

  new_sta_filename=reloc_path + os.sep + "%s_stations.txt"%(stack_time.isoformat())
  logging.info("Station file: %s "%new_sta_filename)
  sta_file=open(new_sta_filename,'w')

  kfiles=[] # This list will contain the filenames that pass the signal to noise test
  logging.info('Number of data files to be read : %d'%len(files))
  for file in files:
    try:
      wf=Waveform()
      wf.read_from_file(file,'SAC',starttime=start_time, endtime=end_time, pad_value=0)
      name=wf.station 
      vals=wf.values
    
      # calculate snr ratio
      mean=np.mean(np.abs(vals))
      try:
        snr=np.max(vals)/mean
      except RuntimeWarning:
        logging.warning('Dividing by zero on SNR calculation (empty trace) : stting SNR to zero.')
        snr=0
      


      # save it for later reference
      snrf.write("%s \t %.2f \n"%(name,snr))

      # if the signal to noise ratio is good
      if snr > snr_limit:
        try:
          # extract corresponding line from station file 
          sta_line=sta_dict[name]
          sta_file.write(sta_line)
          # remember the filename for later
          kfiles.append(file)
        except KeyError:
          logging.error('Cannot find entry for station %s in %s.  Ignoring station.'%(name,stations_filename))
          pass
    except UserWarning:
      logging.error('Reading file %s generated UserWarning - ignoring file.'%file)
      pass

  # close open files
  sta_file.close()
  snrf.close()

  logging.info('Number of files kept (files are ok, and SNR > %.2f) : %d.'%(snr_limit,len(kfiles)))

  
  do_migration_loop_reloc(start_time,end_time,output_dir,kfiles,search_grid_filename,time_grid,options_verbose,options_time)


####

  #st_max=read("%s/stack_max_red_%s.sac"%(stack_path, stack_start_time.isoformat()))
  #st_x=read("%s/stack_x_red_%s.sac"%(stack_path, stack_start_time.isoformat()))
  #st_y=read("%s/stack_y_red_%s.sac"%(stack_path, stack_start_time.isoformat()))
  #st_z=read("%s/stack_z_red_%s.sac"%(stack_path, stack_start_time.isoformat()))

  #tr=st_max.traces[0]
  #logging.debug("Doing filtering")
  #x_filt=lowpass(tr.data,1,1/tr.stats.delta,zerophase=True)
  #st_filt=Stream([Trace(header=tr.stats,data=x_filt)])
  #tr_filt=st_filt.traces[0]
  #tr_filt.downsample(decimation_factor=10, strict_length=False, no_filter=True)
  #st_filt.write("%s/stack_max_red_filt_%s.sac"%(stack_path, stack_start_time.isoformat()), format='SAC')
  #logging.debug("Done!")
  #df=1/tr_filt.stats.delta

  #trigs=trigger.triggerOnset(tr_filt.data,50,50)

  #for trig in trigs:
  #  i_start=trig[0]
  #  i_end=trig[1]+1
  #  i_max_trig=np.argmax(tr_filt.data[i_start:i_end])+i_start
  #  max_trig=tr_filt.data[i_max_trig]
  #  max_trig_95=0.95*max_trig
  #  trigs_95=trigger.triggerOnset(tr_filt.data[i_start:i_end],max_trig_95,max_trig_95)
  #  for trig_95 in trigs_95:
  #    if i_max_trig > trig_95[0]+i_start and i_max_trig < trig_95[1]+i_start:
  #      i_start_95=trig_95[0]+i_start
  #      i_end_95=trig_95[1]+1+i_start
  #  x_mean=np.mean(st_x.traces[0][i_start_95:i_end_95])
  #  x_sigma=np.sqrt(np.var(st_x.traces[0][i_start_95:i_end_95]))
  #  y_mean=np.mean(st_y.traces[0][i_start_95:i_end_95])
  #  y_sigma=np.sqrt(np.var(st_y.traces[0][i_start_95:i_end_95]))
  #  z_mean=np.mean(st_z.traces[0][i_start_95:i_end_95])
  #  z_sigma=np.sqrt(np.var(st_z.traces[0][i_start_95:i_end_95]))

  #  o_time=tr_filt.stats.starttime+i_max_trig*tr_filt.stats.delta
  #  o_err_left=(i_max_trig-i_start_95)*tr_filt.stats.delta
  #  o_err_right=(i_end_95-i_max_trig)*tr_filt.stats.delta
 
  #  logging.info("Max = %.2f, %s - %.2fs + %.2fs, x=%.4f pm %.4f, y=%.4f pm %.4f, z=%.4f pm %.4f\n"%(max_trig,o_time.isoformat(),o_err_left, o_err_right,x_mean,x_sigma,y_mean,y_sigma,z_mean,z_sigma))
  #  reloc_file.write("Max = %.2f, %s - %.2f s + %.2f s, x= %.4f pm %.4f km, y= %.4f pm %.4f km, z= %.4f pm %.4f km\n"%(max_trig,o_time.isoformat(),o_err_left, o_err_right ,x_mean,x_sigma,y_mean,y_sigma,z_mean,z_sigma))

