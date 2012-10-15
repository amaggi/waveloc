#!/usr/bin/env python
# encoding: utf-8

import os, h5py
import numpy as np

from obspy.core import utcdatetime
import time
from mayavi import mlab

from locations_trigger import read_locs_from_file
from correlation import BinaryFile
from NllGridLib import read_stations_file
from hdf5_grids import *

####################################################################################
# Compute theoretical traveltimes and arrival times
def traveltimes(x,y,z,t_orig,stations,time_grids):
  t_th={}
  arr_times={}
  for staname in stations.keys():
    t_th[staname]=[]
    arr_times[staname]=[]
    for i in range(len(x)):
      t_th[staname].append(time_grids[staname].value_at_point(x[i],y[i],z[i]))  # traveltime
      arr_times[staname].append(utcdatetime.UTCDateTime(t_orig[i])+t_th[staname][i]) # arrival time
  return t_th, arr_times
# ----------------------------------------------------------------------------------------
# Compute partial derivatives
def partial_deriv(coord,ev,tth):
  norm=(coord[0]-ev[0])**2+(coord[1]-ev[1])**2+(coord[2]-ev[2])**2
  dpx=-(coord[0]-ev[0])*tth/norm
  dpy=-(coord[1]-ev[1])*tth/norm
  dpz=-(coord[2]-ev[2])*tth/norm
  return [dpx,dpy,dpz]
# ----------------------------------------------------------------------------------------
# Fill G (partial derivatives), d (double differences) and W (weights)
def fill_matrix(cluster,x,y,z,t_orig,stations,t_th,t_arr,coeff,delay,threshold):
  G,W,d=[],[],[]
  N=len(cluster)
  nline,num=0,0
  for staname in stations.keys():
    grid_id="%s.HHZ"%staname
    coord=[stations[staname]['x'],stations[staname]['y'],-stations[staname]['elev']]
    for n in range(N):
      ev1=[x[n],y[n],z[n]]
      e1=cluster[n]
      #dp1=time_grid_deriv.value_at_point(xev1,yev1,zev1,grid_id,type=1)
      dp1 = partial_deriv(coord,ev1,t_th[staname][n]) 
      
      for nn in range(n+1,N):
        e2=cluster[nn]
        if delay[staname][e1-1][e2-1]!='NaN' and coeff[staname][e1-1][e2-1] >= threshold:
          # fill G
          G.append(np.zeros(4*N))
          ev2=[x[nn],y[nn],z[nn]]
          #dp2=time_grid_deriv.value_at_point(xev2,yev2,zev2,grid_id,type=1)
          dp2 = partial_deriv(coord,ev2,t_th[staname][nn])
          dp2=[-elt for elt in dp2]

          G[nline][4*n:4*n+4]=dp1+[1]
          G[nline][4*nn:4*nn+4]=dp2+[-1]
          nline+=1

          # fill d
          obs=delay[staname][e1-1][e2-1]+(t_orig[n]-t_orig[nn])
          theo=t_arr[staname][n]-t_arr[staname][nn]
          d.append(obs-theo)

          # fill W
          W.append(coeff[staname][e1-1][e2-1])
    num+=1
  return G, d, W
# ----------------------------------------------------------------------------------------
# Centroid constraint: sum(delta m)=0
def centroid_constraint(G,d,W):
  ncol=np.size(G,1)
  nline=np.size(G,0)
  for i in range(4):
    G.append(np.zeros(ncol))
    G[nline][i::4]=1
    W.append(0.5)
    nline+=1
  d.extend(np.zeros(4))
  W=np.diag(W)
  return np.matrix(G),np.transpose(np.matrix(d)),np.matrix(W)
# ----------------------------------------------------------------------------------------
# Invert the problem and compute the matrix m
def inversion(G,d,W):
  Gt=np.transpose(G)
  Winv=W.getI()
  a=Gt*Winv*G
  ainv=a.getI()
  m=ainv*Gt*Winv*d
  return m
# ----------------------------------------------------------------------------------------
# Extract the coordinates of the events of a given cluster
def coord_cluster(cluster,locs):
  xini,yini,zini,zini_ph,to_ini=[],[],[],[],[]
  for ind in cluster:
    xini.append(locs[ind-1]['x_mean'])
    yini.append(locs[ind-1]['y_mean'])
    zini.append(locs[ind-1]['z_mean'])
    zini_ph.append(-locs[ind-1]['z_mean']) # positive z axis upwards
    to_ini.append(locs[ind-1]['o_time'])
  return xini, yini, zini, zini_ph, to_ini
# ----------------------------------------------------------------------------------------
# Plot old and new locations
def plot_events(cluster,locs,stations,x,y,z,i,threshold,nbmin,area):
  # Stations coordinates
  xsta,ysta,zsta=[],[],[]
  for sta in stations.keys():
    xsta.append(stations[sta]['x'])
    ysta.append(stations[sta]['y'])
    zsta.append(stations[sta]['elev'])

  # Initial hypocentral parameters
  xini, yini, zini, zini_ph, to_ini = coord_cluster(cluster[i],locs)
  z_ph=[-elt for elt in z]

  s=mlab.figure(i,bgcolor=(1,1,1),fgcolor=(0,0,0),size=(1000,900))
  mlab.clf()
  s=mlab.points3d(xini,yini,zini_ph,color=(1,1,0),scale_factor=0.2) # yellow : initial locations
  s=mlab.points3d(xsta,ysta,zsta,color=(1,0,0),scale_factor=0.05,mode='cube')
  s=mlab.points3d(x,y,z_ph,color=(0,1,1),scale_factor=0.2) # cyan : new locations
  s=mlab.axes(extent=area,color=(0,0,0)) # axe des z positif vers le haut
  s=mlab.outline(extent=area,color=(0,0,0))
  s=mlab.title("cluster=%s, threshold=%s, nbmin=%s"%(i,threshold,nbmin),height=0.1,size=0.35,color=(0,0,0))
  mlab.show()

####################################################################################
def do_double_diff(x,y,z,to,stations,coeff,delay,cluster,threshold,t_th, arr_times):
    N=len(cluster)

    # Fill G, d and W
    G, d, W=fill_matrix(cluster,x,y,z,to,stations,t_th,arr_times,coeff,delay,threshold)

    # Centroid constraint : add 4 lines to G, d and W
    G,d,W=centroid_constraint(G,d,W)

    # Inversion
    m=inversion(G,d,W)

    for i in range(N):
      x[i]=x[i]+m[4*i,0]
      y[i]=y[i]+m[4*i+1,0]
      z[i]=z[i]+m[4*i+2,0]
      to[i]=utcdatetime.UTCDateTime(to[i])+m[4*i+3,0]

    return x,y,z,to
####################################################################################
def do_double_diff_setup_and_run(opdict):

  base_path=opdict['base_path']
  verbose=opdict['verbose']
  dd_loc=opdict['dd_loc']

  # Station
  stations_filename=os.path.join(base_path,'lib',opdict['stations'])
  stations=read_stations_file(stations_filename)

  # Output directory
  output_dir=os.path.join(base_path,'out',opdict['outdir'])

  # Location file
  locdir=os.path.join(base_path,'out',opdict['outdir'],'loc')
  loc_filename=os.path.join(locdir,'locations.dat')
  locs=read_locs_from_file(loc_filename)


  # ----------------------------------------------------------------------------------------
  # search grid
  search_grid_filename=os.path.join(base_path,'lib',opdict['search_grid'])
  # traveltimes grid
  grid_filename_base=os.path.join(base_path,'lib',opdict['time_grid'])
  grid_info=read_hdr_file(search_grid_filename)
  time_grids=get_interpolated_time_grids(opdict)
  
  # Extract the UTM coordinates of the area of study
  xstart=grid_info['x_orig']
  xend=xstart+grid_info['nx']*grid_info['dx']
  ystart=grid_info['y_orig']
  yend=ystart+grid_info['ny']*grid_info['dy']
  zend=-grid_info['z_orig']
  zstart=-(-zend+grid_info['nz']*grid_info['dz'])
  area=[xstart,xend,ystart,yend,zstart,zend]
  # ----------------------------------------------------------------------------------------
  nbmin=int(opdict['nbsta'])
  threshold=float(opdict['clus'])

  # Correlation, time delay and cluster files
  corr_file=os.path.join(locdir,opdict['corr'])
  cfile=BinaryFile(corr_file)
  coeff=cfile.read_binary_file()

  delay_file=os.path.join(locdir,opdict['delay'])
  dfile=BinaryFile(delay_file)
  delay=dfile.read_binary_file()

  cluster_file=os.path.join(locdir,'cluster-%s-%s'%(str(threshold),str(nbmin)))
  clfile=BinaryFile(cluster_file)
  cluster=clfile.read_binary_file()
  # ----------------------------------------------------------------------------------------
  # Input parameters
  nb_iter=2
  len_cluster_min=2

  # ----------------------------------------------------------------------------------------
  for i in cluster.keys():
    print "CLUSTER %d:"%i,cluster[i],len(cluster[i])
    iter = 0
    N = len(cluster[i])

    # Hypocentral parameters to be changed
    x,y,z,z_ph,to = coord_cluster(cluster[i],locs)

    if N > len_cluster_min:

      # Theroretical traveltimes and arrival times
      t_th, arr_times = traveltimes(x,y,z,to,stations,time_grids)
      
      x,y,z,to = do_double_diff(x,y,z,to,stations,coeff,delay,cluster[i],threshold,t_th, arr_times)

      if verbose:
        plot_events(cluster,locs,stations,x,y,z,i,threshold,nbmin,area)

    if dd_loc:
      ind=0
      for j in cluster[i]:
          locs[j-1]['x_mean']=x[ind]
          locs[j-1]['y_mean']=y[ind]
          locs[j-1]['z_mean']=z[ind]
          locs[j-1]['o_time']=to[ind]
          locs[j-1]['x_sigma']=0
          locs[j-1]['y_sigma']=0
          locs[j-1]['z_sigma']=0
          locs[j-1]['o_err_right']=0
          locs[j-1]['o_err_left']=0
          ind+=1
   
  if dd_loc:
    new_loc_filename=os.path.join(locdir,'relocations.dat')
    new_loc_file=open(new_loc_filename,'w')
    for loc in locs:
      new_loc_file.write("Max = %.2f, %s - %.2f s + %.2f s, x= %.4f pm %.4f km, y= %.4f pm %.4f km, z= %.4f pm %.4f km\n"%(loc['max_trig'],loc['o_time'].isoformat(),loc['o_err_left'], loc['o_err_right'],loc['x_mean'],loc['x_sigma'],loc['y_mean'],loc['y_sigma'],loc['z_mean'],loc['z_sigma']))
    new_loc_file.close()


###############################################################################################
if __name__ == '__main__':
  from options import WavelocOptions

  logging.basicConfig(level=logging.INFO, format="%(levelname)s : %(asctime)s : %(message)s")

  wo = WavelocOptions()
  args=wo.p.parse_args()

  wo.set_all_arguments(args)
  wo.verify_doublediff_options()

  do_double_diff_setup_and_run(wo.opdict)
