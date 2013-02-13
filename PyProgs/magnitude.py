#!/usr/bin/env python
# encoding: utf-8

import sys,glob
import numpy as np
from obspy.xseed import Parser
from obspy.signal import pazToFreqResp,estimateMagnitude
import matplotlib.pyplot as plt
from OP_waveforms import *
from locations_trigger import read_locs_from_file,read_header_from_file
from NllGridLib import read_stations_file

# =======================================================
def read_paz(verbose=False):
  files=glob.glob('/home/nadege/Desktop/Dataless/*.dataless')
  files.sort()
  paz={}

  for file in files:
    p = Parser(file)
    blk=p.blockettes

    for j in range(len(blk[50])):

      mult = len(blk[58])/len(blk[52])

      sta=blk[50][j].station_call_letters
      paz[sta]={}

      if j >= 2:
        j=j+1

      for i in range(j*3,len(blk[52])):
        channel=blk[52][i].channel_identifier
        if channel != 'HDF' and channel != 'HDT' and channel != 'HDA':
          paz[sta][channel]={}
          paz[sta][channel]['poles']=np.array(blk[53][i].real_pole)+1j*np.array(blk[53][i].imaginary_pole)
          paz[sta][channel]['zeros']=np.array(blk[53][i].real_zero)+1j*np.array(blk[53][i].imaginary_zero)
          paz[sta][channel]['gain']=blk[53][i].A0_normalization_factor
          paz[sta][channel]['sensitivity']=blk[58][(i+1)*mult-1].sensitivity_gain
        if 'HHE' in paz[sta].keys() and 'HHN' in paz[sta].keys() and 'HHZ' in paz[sta].keys():
          break

      if verbose:
        h,f=pazToFreqResp(paz[sta]['HHZ']['poles'],paz[sta]['HHZ']['zeros'],paz[sta]['HHZ']['gain'],0.01,16384,freq=True)
        fig=plt.figure()
        fig.set_facecolor('white')
        plt.loglog(f,np.abs(h))
        plt.title(sta)
        plt.xlabel('Frequency (Hz)')
        plt.ylabel('Amplitude')
        plt.show()

  return paz
# =======================================================
def fill_values(vals,tdeb,data_glob,data_dir,comp):
  if comp == 'HHE' or comp == 'HHN':
    data_dir=os.path.join(data_dir,comp)
  data_files=glob.glob(os.path.join(data_dir,data_glob))
  data_files.sort()

  for datafile in data_files:
    wf=Waveform()
    wf.read_from_file(datafile)
    vals[wf.station][comp]=wf.values
    tdeb[wf.station][comp]=wf.starttime
    dt=wf.delta

  return vals,tdeb,dt
# =======================================================
def do_comp_mag(opdict):

  base_path=opdict['base_path']
  verbose=opdict['verbose']

  paz=read_paz(verbose)
  #verbose=True

  # output directory
  output_dir=os.path.join(base_path,'out',opdict['outdir'])

  # data
  data_dir=os.path.join(base_path,'data',opdict['datadir'])
  data_glob=opdict['dataglob']

  # location file
  locdir=os.path.join(base_path,'out',opdict['outdir'],'loc')
  locfile=os.path.join(locdir,'locations.dat')
  locs=read_locs_from_file(locfile)
  header=read_header_from_file(os.path.join(locdir,'locations.dat'))
  snr_wf=np.float(header['wf snr'])

  # Stations
  stations_filename=os.path.join(base_path,'lib',opdict['stations'])
  stations=read_stations_file(stations_filename)

  vals,tdeb={},{}
  for sta in sorted(stations):
    vals[sta]={}
    tdeb[sta]={}

  cha_list=['HHZ']
  vals,tdeb,dt=fill_values(vals,tdeb,data_glob,data_dir,'HHZ')
  if os.path.isdir(os.path.join(data_dir,'HHE')):
    vals,tdeb,dt=fill_values(vals,tdeb,data_glob,data_dir,'HHE')
    cha_list.append('HHE')
  if os.path.isdir(os.path.join(data_dir,'HHN')):
    vals,tdeb,dt=fill_values(vals,tdeb,data_glob,data_dir,'HHN')
    cha_list.append('HHN')

  new_file=open(os.path.join(locdir,'locations_mag.dat'),'w')

  mags=[]
  for loc in locs:
    stack_time=loc['o_time']
    loc_x=loc['x_mean']
    loc_y=loc['y_mean']
    loc_z=-loc['z_mean']

    ml=[]
    for sta in sorted(stations):
      if vals[sta]:
        paz_list,p2p_amp,tspan=[],[],[]
        h_dist=np.sqrt((loc_x-stations[sta]['x'])**2+(loc_y-stations[sta]['y'])**2+(loc_z-stations[sta]['elev'])**2)
        for cha in cha_list:
          istart=int(round(stack_time-0.5-tdeb[sta][cha])*1./dt)
          iend=int(round(stack_time+5.5-tdeb[sta][cha])*1./dt)

          x=vals[sta][cha][istart:iend+1]

          if x.any() and np.max(x)/np.mean(np.abs(x)) > snr_wf:
            max_amp=np.max(x)
            i_max_amp=np.argmax(x)
            min_amp=np.abs(np.min(x))
            i_min_amp=np.argmin(x)
            
            paz_list.append(paz[sta][cha])
            tspan.append(np.abs(i_max_amp-i_min_amp)*dt)
            p2p_amp.append(max_amp+min_amp)

            if verbose:
              fig=plt.figure()
              fig.set_facecolor('white')
              plt.plot(x)
              plt.plot(i_min_amp,x[i_min_amp],'ro')
              plt.plot(i_max_amp,x[i_max_amp],'ro')
              plt.title('%s,%s'%(sta,cha))
              plt.show()

        if paz_list:
          mag=estimateMagnitude(paz_list,p2p_amp,tspan,h_dist)
          ml.append(mag)

    new_file.write("Max = %.2f, %s - %.2f s + %.2f s, x= %.4f pm %.4f km, y= %.4f pm %.4f km, z= %.4f pm %.4f km, ml= %.2f pm %.2f\n"%(loc['max_trig'],loc['o_time'].isoformat(),loc['o_err_left'], loc['o_err_right'],loc['x_mean'],loc['x_sigma'],loc['y_mean'],loc['y_sigma'],loc['z_mean'],loc['z_sigma'],np.mean(ml),np.std(ml)))

    mags.append(np.mean(ml))

  new_file.close()
  print np.max(mags),locs[np.argmax(mags)]['o_time']

  fig=plt.figure()
  fig.set_facecolor('white')
  plt.hist(mags,25)
  plt.xlabel('Magnitude')

  r=np.arange(-1.2,1.8,0.1)
  N=[]
  for i in r:
    N.append(len(np.where(mags >= i)[0]))

  i1=np.max(np.where(np.log(N) == np.max(np.log(N))))
  i2=np.argmin(np.log(N))
  p=np.polyfit(r[i1:i2],np.log(N)[i1:i2],deg=1)
  print "b-value:",-p[0] 

  fig=plt.figure()
  fig.set_facecolor('white')
  plt.plot(r,np.log(N))
  plt.plot(r[i1:i2],np.polyval(p,r[i1:i2]),'r')
  #plt.yscale('log')
  plt.xlabel('Magnitude')
  plt.ylabel('log N')
  plt.title('Gutenberg Richter law')
  plt.show()
# =======================================================
if __name__ == '__main__' :

  from options import WavelocOptions
  logging.basicConfig(level=logging.INFO, format='%(levelname)s : %(asctime)s : %(message)s')

  wo = WavelocOptions()
  args=wo.p.parse_args()

  wo.set_all_arguments(args)

  do_comp_mag(wo.opdict)
