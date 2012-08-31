#!/usr/bin/env python
# encoding: utf-8

import os, sys, optparse

import numpy as np


from OP_waveforms import *

from grids_paths import *
from time import time, sleep
from sub_PdF_waveloc import do_migration_loop_continuous
import Queue, threading


import logging
logging.basicConfig(level=logging.INFO, format='%(levelname)s : %(asctime)s : %(message)s')


class myThread (threading.Thread):
  def __init__(self, threadID, name, q):
    threading.Thread.__init__(self)
    self.threadID = threadID
    self.name = name
    self.q = q
  def run(self):
    logging.info("Starting " + self.name)
    process_data(self.name, self.q)
    logging.info("Exiting " + self.name)

def process_data(threadName, q):
  while not exitFlag:
    queueLock.acquire()
    logging.debug('%s acquired lock'%threadName)
    if not workQueue.empty():

      logging.debug('%s going to start work'%threadName)
      #(start_time, end_time, data_dir, output_dir, data_glob, search_grid_filename, time_grid, options.verbose, options.time) = q.get()
      (start_time, end_time, data_dir, output_dir, data_glob, search_grid_filename, options.verbose, options.time) = q.get()
      queueLock.release()
      logging.debug('%s released lock'%threadName)

      logging.info("%s doing migration on slice %s - %s " % (threadName, start_time.isoformat(), end_time.isoformat()))
      do_migration_loop_continuous(start_time, end_time, data_dir, output_dir, data_glob, search_grid_filename, time_grid, options.verbose, options.time)

    else:
      queueLock.release()
      logging.debug('%s released lock'%threadName)

    sleep(1)




#get waveloc path from environment
base_path=os.getenv('WAVELOC_PATH')

# setup arguments for command line parsing
p = optparse.OptionParser()
p.add_option('--nthreads','-n', action='store',type='int', default=1,help="number of threads")
p.add_option('--time','-t',action='store_true',help='print timing information to stout')
p.add_option('--verbose','-v',action='store_true',default=False,help='print debugging information to stout')
p.add_option('--2D',action='store_true',default=False,dest='twoD',help='use 2D time grids')
p.add_option('--time_grid',action='store',help="time grid basename e.g. belgium.P (found in $WAVELOC_PATH/lib)")
p.add_option('--search_grid',action='store',help="search grid e.g. grid.500m.search.hdr (found in $WAVELOC_PATH/lib)")
p.add_option('--stations','-s',action='store',help='station list (found in $WAVELOC_PATH/lib)')
p.add_option('--outdir','-o',action='store',help="output subdirectory")
p.add_option('--datadir',action='store',help="data subdirectory")
p.add_option('--dataglob',action='store',help="data glob")
p.add_option('--starttime',action='store',help="start time for data e.g. 2010-10-14T00:00:00.0Z")
p.add_option('--endtime',action='store',help="end time for data e.g. 2010-10-14T10:00:00.0Z")
p.add_option('--data_length',action='store',help="length in seconds for data segments to analyse (e.g. 630)")
p.add_option('--data_overlap',action='store',help="length in seconds for overlapping data segments (e.g. 30)")
p.add_option('--load_ttimes_buf',action='store_true',default=False,dest='load_buf',help='load pre-calculated travel-times for the search grid from file')



# parse command line
(options,arguments)=p.parse_args()

# set variables from command line options
# grid
logging.debug(base_path)
logging.debug(options.time_grid)
grid_filename_base=os.path.join(base_path,'lib',options.time_grid)
search_grid_filename=os.path.join(base_path,'lib',options.search_grid)
# stations
stations_filename=os.path.join(base_path,'lib',options.stations)
# output directory
if options.outdir==None: raise UserWarning("No output subdirectory supplied")
output_dir=os.path.join(base_path,'out',options.outdir)
if not os.path.exists(output_dir):
  os.makedirs(output_dir)
# data
if options.datadir==None: raise UserWarning("No data subdirectory supplied")
if options.dataglob==None: raise UserWarning("No data glob supplied")
data_dir=os.path.join(base_path,'data',options.datadir)
data_glob=options.dataglob


if options.verbose: 
  print ""
  print "Input parameters:"
  print "-----------------"
  print "Grid        = %s"%(grid_filename_base)
  print "Stations    = %s"%(stations_filename)
  print "Data        = %s"%(os.path.join(data_dir,data_glob))
  print ""
  print "Output parameters:"
  print "-----------------"
  print "Out dir     = %s"%(output_dir)
  

#raise UserWarning ('Stop here')
#######################################################################
#                       START PROCESSING 
#######################################################################

print ""
print "----------------"
print "START PROCESSING"
print "----------------"
print ""

# Create Obspy streams for output

#  ***** reading station file ******

if options.verbose:
  logging.info("Reading station file")

if options.time:
  t_ref=time()

sta=StationList()
sta.read_from_file(stations_filename)

if options.time:
  t=time()-t_ref
  logging.info("Time for reading %d stations from file : %.4f s\n" % (sta.nsta,t))

datafile_list=glob.glob(os.path.join(data_dir,data_glob))

cha=ChannelList()
cha.populate_from_station_list_and_data_files(sta,datafile_list)


########### DEAL WITH START AND END TILES ############

initial_start_time=utcdatetime.UTCDateTime(options.starttime)
initial_end_time=initial_start_time+np.float(options.data_length)

final_end_time=utcdatetime.UTCDateTime(options.endtime)

time_shift_secs=np.float(options.data_length)-np.float(options.data_overlap)



######### INTERPOLATE TRAVEL TIMES #############

# The time grid will contain as array values just the travel-times needed 
# (interpolated from the full NLL files) so we can free up the memory as soon as possible

if options.verbose:
  logging.info("Extracting useful travel-times")

if options.time:
  t_ref=time()  

time_grid=QDTimeGrid()
time_grid.read_NLL_hdr_file(search_grid_filename)
if options.twoD:
  time_grid.populate_from_2D_time_grids(grid_filename_base,cha)
else:
  time_grid.populate_from_time_grids(grid_filename_base,cha,output_dir,options.load_buf)

if options.time:
  t=time()-t_ref
  n_times=time_grid.nx*time_grid.ny*time_grid.nz+cha.ncha
  logging.info("Time for extracting and saving %dx%dx%dx%d=%d travel-times : %.2f s\n" % (time_grid.nx,time_grid.ny,time_grid.nz,cha.ncha,n_times,t))

# number of time slices for queue dimension (estimate large)
number_of_time_slices=np.int(np.fix((final_end_time - initial_start_time)/time_shift_secs))+10

# set up threads
threadList=["Thread-%d"%i for i in range(options.nthreads)]
queueLock = threading.Lock()
workQueue = Queue.Queue(number_of_time_slices)
threads = []
threadID = 0

exitFlag = 0

# Create new threads
for tName in threadList:
  threadID += 1
  thread = myThread(threadID, tName, workQueue)
  thread.start()
  threads.append(thread)


# Fill the queue
queueLock.acquire()

# READ DATA

# start loop over time
start_time=initial_start_time
end_time=initial_end_time

while (start_time < final_end_time):

  # ADD MIGRATION WORK TO QUEUE !!
  logging.info('Added time-slice %s - %s to work-queue'%(start_time.isoformat(),end_time.isoformat()))
  #workQueue.put((start_time, end_time, data_dir, output_dir, data_glob, search_grid_filename, time_grid, options.verbose, options.time))
  workQueue.put((start_time, end_time, data_dir, output_dir, data_glob, search_grid_filename, options.verbose, options.time))


  # Reset the start and end times to loop again
  start_time=start_time+time_shift_secs
  end_time=end_time+time_shift_secs


# end of add-work loop
queueLock.release()

# Wait for queue to empty
while not workQueue.empty():
  pass

# Notify threads it's time to exit
exitFlag = 1

# Wait for all threads to complete
for t in threads:
  t.join()
logging.info("Exiting Main Thread")




