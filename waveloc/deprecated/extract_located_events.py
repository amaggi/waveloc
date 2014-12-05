#!/usr/bin/env python
# encoding: utf-8

import os, sys, optparse, glob
from obspy.core import *
import matplotlib.pyplot as plt
import numpy as np
from NLL_IO import qd_read_hyp_file, qd_read_picks_from_hyp_file


# get path
base_path=os.getenv('WAVELOC_PATH')

# Read command line

p=optparse.OptionParser()
p.add_option('--outdir', '-o', action='store', help='output subdirectory')
p.add_option('--datadir', '-d', action='store', help='data subdirectory')
p.add_option('--data_glob',action='store',help="data glob")
p.add_option('--left_pad', action='store', help='seconds to pad on left of origin time')
p.add_option('--right_pad', action='store', help='seconds to pad on right of origin time')


(options,arguments)=p.parse_args()


lib_path=base_path + os.sep + "lib"
loc_path="%s/out/%s/loc"%(base_path,options.outdir)
data_path="%s/data/%s"%(base_path,options.datadir)
extract_path="%s/out/%s/extract"%(base_path,options.outdir)
if not os.path.exists(extract_path):
  os.makedirs(extract_path)


loc_filename="%s/locations.dat"%(loc_path)

data_files=glob.glob(data_path + os.sep + options.data_glob)
data_files.sort()

#extract_files=[]
#for data_file in data_files:
#  filename=extract_path + os.sep + data_file.split(os.sep)[-1]
#  print filename


station_names=[data_file.split(os.sep)[-1].split('.')[0] for data_file in data_files]

print station_names

for sta in station_names:
  subdir=extract_path + os.sep + sta
  if not os.path.exists(subdir):
    os.makedirs(subdir)


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


for loc in locs[0:1]:
#for loc in locs:

  stack_time=loc[1]
  # set start and end time of plot
  start_time=stack_time-np.float(options.left_pad)
  end_time=stack_time+np.float(options.right_pad)


  for ifile in range(len(data_files)):
    sta=station_names[ifile]
    data_file=data_files[ifile]
    base_file_name=data_file.split(os.sep)[-1]
    outfilename=extract_path + os.sep + sta + os.sep + base_file_name + '_'+stack_time.isoformat()
    print outfilename
    st=read(data_files[ifile],starttime=start_time, endtime=end_time)
    tr=st.traces[0]
    tr.write(outfilename,format='SAC')
      

