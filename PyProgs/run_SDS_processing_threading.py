#!/usr/bin/env python
# encoding: utf-8

import os, sys, optparse, glob
import numpy as np
from OP_waveforms import *
import logging
import Queue, threading, time

class myThread (threading.Thread):
  def __init__(self, threadID, name, q):
    threading.Thread.__init__(self)
    self.threadID = threadID
    self.name = name
    self.q = q
  def run(self):
    print "Starting " + self.name
    process_data(self.name, self.q)
    print "Exiting " + self.name

def process_data(threadName, q):
  while not exitFlag:
    queueLock.acquire()
    if not workQueue.empty():

      (data_dir,net,sta,comp,starttime,endtime) = q.get()
      queueLock.release()

      filt_filename=os.path.join(data_dir,"%s.%s.%s.%s.filt.mseed"%(start_time.isoformat(),net,sta,comp))
      logging.info("%s processing to create %s" % (threadName, filt_filename))

      wf=Waveform()
      try:
        wf.read_from_SDS(data_dir,net,sta,comp, starttime=start_time, endtime=end_time)

        wf.bp_filter(filter_c1, filter_c2,rmean=True,taper=True)
        if options.resample:
          wf.resample(np.float(options.fs))

        wf.write_to_file_filled(filt_filename,format='MSEED',fill_value=0)

        kurt_filename=os.path.join(data_dir,"%s.%s.%s.%s.filt_kurt.mseed"%(start_time.isoformat(),net,sta,comp))
        logging.info("%s processing to create %s " % (threadName, kurt_filename))
        wf.process_kurtosis(kurt_window,recursive=options.krec,pre_taper=True, post_taper=True)
        wf.write_to_file_filled(kurt_filename,format='MSEED',fill_value=0)
        if options.kderiv:
          kurt_grad_filename=os.path.join(data_dir,"%s.%s.%s.%s.filt_kurt_grad.mseed"%(start_time.isoformat(),net,sta,comp))
          logging.info("%s processing to create %s " % (threadName, kurt_grad_filename))
          wf.take_positive_derivative(pre_taper=True,post_taper=True)
          wf.write_to_file_filled(kurt_grad_filename,format='MSEED',fill_value=0)
  
      except UserWarning:
        logging.info('No data within time limits for %s %s %s'%(net,sta,comp))

    else:
      queueLock.release()

    time.sleep(1)



logging.basicConfig(level=logging.DEBUG, format='%(levelname)s : %(asctime)s : %(message)s')


#get waveloc path from environment
base_path=os.getenv('WAVELOC_PATH_PDF')

p = optparse.OptionParser()
p.add_option('--nthreads','-n', action='store',type='int', default=1,help="number of threads")
p.add_option('--datadir',action='store',help="data subdirectory")
p.add_option('--net_list',action='store',help="list of network codes (e.g. \"BE,G\") ")
p.add_option('--sta_list',action='store',help="list of station names (e.g. \"STA1,STA2\") ")
p.add_option('--comp_list',action='store',help="list of component names (e.g. \"HHZ,LHZ\") ")
p.add_option('--starttime',action='store',help="start time for data e.g. 2010-10-14T00:00:00.0Z")
p.add_option('--endtime',action='store',help="end time for data e.g. 2010-10-14T10:00:00.0Z")
p.add_option('--resample',action='store_true',default=False, help="resample data")
p.add_option('--fs',action='store',default=False, help="resample frequency")

p.add_option('--c1',action='store',type='float',help="low frequency corner of band pass filter ")
p.add_option('--c2',action='store',type='float',help="high frequency corner of band pass filter ")
p.add_option('--kwin',action='store',type='float',help="length of kurtosis window (seconds)")
p.add_option('--krec',action='store_true',default=False, help="use recursive kurtosis calculation (faster but less precise)")
p.add_option('--kderiv',action='store_true',default=False, help="use derivative of kurtosis")

(options,arguments)=p.parse_args()

data_dir=os.path.join(base_path,'data',options.datadir)

# if typing of options does not work
#filter_c1=np.float(options.c1)
#filter_c2=np.float(options.c2)
#kurt_window=np.float(options.kwin)

filter_c1=options.c1
filter_c2=options.c2
kurt_window=options.kwin

# start and end time to process
start_time=utcdatetime.UTCDateTime(options.starttime)
end_time=utcdatetime.UTCDateTime(options.endtime)

logging.debug(options.net_list)
logging.debug(options.sta_list)
logging.debug(options.comp_list)

net_list=options.net_list.split(',')
sta_list=options.sta_list.split(',')
comp_list=options.comp_list.split(',')

n_max_work=len(net_list)*len(sta_list)*len(comp_list)

# set up threads
threadList=["Thread-%d"%i for i in range(options.nthreads)]
queueLock = threading.Lock()
workQueue = Queue.Queue(n_max_work)
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

for net in net_list:
  for sta in sta_list:
    for comp in comp_list:
      full_path=os.path.join(data_dir,net,sta,"%s.D"%comp)
      logging.debug("Full path : %s"%full_path)
      if os.path.exists(full_path):
        # call add info to queue
        workQueue.put((data_dir,net,sta,comp,start_time,end_time))

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

