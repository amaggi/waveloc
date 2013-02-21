#!/usr/bin/env python
# encoding: utf-8

import sys,glob
import numpy as np
from obspy.xseed import Parser
from obspy.signal import pazToFreqResp,estimateMagnitude
import matplotlib.pyplot as plt
from OP_waveforms import *
from locations_trigger import read_locs_from_file,read_header_from_file, write_header_options
from NllGridLib import read_stations_file

# =======================================================
def read_paz(files):
  """
  Read dataless and extract poles and zeros
  """
  paz={}

  for file in files:
    p = Parser(file)
    blk=p.blockettes

    for j in range(len(blk[50])):

      mult = len(blk[58])/len(blk[52])

      sta=blk[50][j].station_call_letters
      paz[sta]={}

      for i in range(j*3,len(blk[52])):
        channel=blk[52][i].channel_identifier
        paz[sta][channel]={}
        paz[sta][channel]['poles']=np.array(blk[53][i].real_pole)+1j*np.array(blk[53][i].imaginary_pole)
        paz[sta][channel]['zeros']=np.array(blk[53][i].real_zero)+1j*np.array(blk[53][i].imaginary_zero)
        paz[sta][channel]['gain']=blk[53][i].A0_normalization_factor
        paz[sta][channel]['sensitivity']=blk[58][(i+1)*mult-1].sensitivity_gain

  return paz
# =======================================================
def fill_values(vals,tdeb,data_glob,data_dir,comp):
  """
  Create a dictionnary with all values for each channel of each station
  """
  data_files=glob.glob(os.path.join(data_dir,data_glob))
  if len(data_files) == 0:
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
def bvalue(mag,r):
  """
  Compute the b-value by a simple linear fitting
  """
  N=[]
  for i in r:
    N.append(len(np.where(mag >= i)[0]))

  i1=np.min(np.where(np.log10(N) <= np.max(np.log10(N))-0.01)[0])
  i2=np.argmin(np.log10(N))
  p=np.polyfit(r[i1:i2],np.log10(N)[i1:i2],deg=1)

  return p,np.log10(N),i1,i2
# =======================================================
def do_comp_mag(opdict):

  base_path=opdict['base_path']
  verbose=opdict['verbose']

  # dataless
  dataless_glob=glob.glob(os.path.join(base_path,'lib',opdict['dataless']))
  dataless_glob.sort()

  # output directory
  output_dir=os.path.join(base_path,'out',opdict['outdir'])

  # data
  data_dir=os.path.join(base_path,'data',opdict['datadir'])
  data_glob=opdict['dataglob']

  # location file
  locdir=os.path.join(base_path,'out',opdict['outdir'],'loc')
  locfile=os.path.join(locdir,'locations.dat')
  locs=read_locs_from_file(locfile)
  opdict=read_header_from_file(locfile,opdict)
  snr_wf=np.float(opdict['snr_tr_limit'])

  # Stations
  stations_filename=os.path.join(base_path,'lib',opdict['stations'])
  stations=read_stations_file(stations_filename)

  paz=read_paz(dataless_glob)

  vals,tdeb={},{}
  for sta in sorted(stations):
    vals[sta]={}
    tdeb[sta]={}

  cha_list=opdict['comp_list']
  for cha in cha_list:
    vals,tdeb,dt=fill_values(vals,tdeb,data_glob,data_dir,cha)

  new_file=open(locfile,'w')
  write_header_options(new_file,opdict)

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

    if ml:
      mags.append(np.mean(ml))

  new_file.close()

  r=np.arange(-3,3,0.1)
  p,logN,i1,i2 = bvalue(mags,r)
  print "b-value:",-p[0] 


  fig=plt.figure(figsize=(10,5))
  fig.set_facecolor('white')
  ax1 = fig.add_subplot(121)
  ax1.hist(mags,25)
  ax1.set_xlabel('Magnitude')

  ax2 = fig.add_subplot(122,title='Gutenberg Richter law')
  ax2.plot(r,logN)
  ax2.plot(r[i1:i2],np.polyval(p,r[i1:i2]),'r')
  ax2.set_xlabel('Magnitude')
  ax2.set_ylabel('log N')
  plt.show()
# =======================================================
if __name__ == '__main__' :

  from options import WavelocOptions
  logging.basicConfig(level=logging.INFO, format='%(levelname)s : %(asctime)s : %(message)s')

  wo = WavelocOptions()
  args=wo.p.parse_args()

  wo.set_all_arguments(args)
  wo.verify_magnitude_options()

  do_comp_mag(wo.opdict)
