#!/usr/bin/env python
# encoding: utf-8

"""
Classes pertaining to IO on NLL type files

Created by Alessia Maggi and Alberto Michelini.  

"""
import os, sys, glob 
import numpy as np
from obspy.core import utcdatetime

def qd_read_hyp_file(filename):
  f=open(filename, 'r')
  lines=f.readlines()
  f.close()
  
  for line in lines:
    words=line.split()
    try:
      if words[0]=='HYPOCENTER':
        hypo_x=np.float(words[2])
        hypo_y=np.float(words[4])
        hypo_z=np.float(words[6])
      if words[0]=='GEOGRAPHIC':
        year=np.int(words[2])
        month=np.int(words[3])
        day=np.int(words[4])
        hour=np.int(words[5])
        minute=np.int(words[6])
        seconds=np.float(words[7])
        otime=utcdatetime.UTCDateTime(year,month,day,hour,minute,seconds)
      if words[0]=='STATISTICS':
        sigma_x=np.sqrt(np.float(words[8]))
        sigma_y=np.sqrt(np.float(words[14]))
        sigma_z=np.sqrt(np.float(words[18]))
    except IndexError:
      pass

  return (otime, hypo_x, sigma_x, hypo_y, sigma_y, hypo_z, sigma_z)

def qd_read_picks_from_hyp_file(filename):
  f=open(filename, 'r')
  lines=f.readlines()
  f.close()
  
  
  for iline in range(len(lines)):
    line=lines[iline]
    words=line.split()
    if words[0]=='PHASE':
      iline_phase=iline
      break
  
  phases={}
  for line in lines[iline+1:]:
    words=line.split()
    try:
      if words[4]=='P':
        station=words[0]
        year=np.int(words[6][0:4])
        month=np.int(words[6][4:6])
        day=np.int(words[6][6:8])
        hour=np.int(words[7][0:2])
        minute=np.int(words[7][2:4])
        seconds=np.float(words[8])
        ptime=utcdatetime.UTCDateTime(year,month,day,hour,minute,seconds)
        phases[station]=ptime
    except IndexError:
      pass
     
  return phases

