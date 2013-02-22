#!/usr/bin/env python
# encoding: utf-8

import os, h5py
import numpy as np
import logging

from obspy.core import utcdatetime
import time

from locations_trigger import read_locs_from_file, read_header_from_file, write_header_options
from correlation import BinaryFile
from NllGridLib import read_stations_file
from hdf5_grids import *

####################################################################################
# Compute theoretical traveltimes and arrival times
def traveltimes(x,y,z,t_orig,stations,time_grids):
  t_th={}
  arr_times={}
  for staname in sorted(stations):
    if not staname in time_grids.keys():
      logging.info("%s station not in time_grids"%staname)
      continue
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

  for staname in sorted(stations):
    if not staname in delay.keys():
      continue
    if not staname in t_th.keys():
      continue
    grid_id="%s.HHZ"%staname
    coord=[stations[staname]['x'],stations[staname]['y'],-stations[staname]['elev']]
    for n in range(N):
      ev1=[x[n],y[n],z[n]]
      e1=cluster[n]
      dp1 = partial_deriv(coord,ev1,t_th[staname][n]) 

      for nn in range(n+1,N):
        e2=cluster[nn]
        if delay[staname][e1-1][e2-1]!='NaN' and coeff[staname][e1-1][e2-1] >= threshold: 
          # fill G
          G.append(np.zeros(4*N))
          ev2=[x[nn],y[nn],z[nn]]
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
# Inversion: computation of the matrix m
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
def plot_events(cluster,locs,stations,x,y,z,i,threshold,nbmin,area,nbsta):
  from mayavi import mlab

  # Stations coordinates
  xsta,ysta,zsta=[],[],[]
  for sta in sorted(stations):
    xsta.append(stations[sta]['x'])
    ysta.append(stations[sta]['y'])
    zsta.append(stations[sta]['elev'])

  z_ph=[-elt for elt in z]

  # Initial hypocentral parameters
  xini, yini, zini, zini_ph, to_ini = coord_cluster(cluster[i],locs)

  s=mlab.figure(i,bgcolor=(1,1,1),fgcolor=(0,0,0),size=(1000,900))
  mlab.clf()
  s=mlab.points3d(xini,yini,zini_ph,color=(1,1,0),scale_factor=0.2) # yellow : initial locations
  s=mlab.points3d(xsta,ysta,zsta,color=(1,0,0),scale_factor=0.05,mode='cube')
  s=mlab.points3d(x,y,z_ph,color=(0,1,1),scale_factor=0.2) # cyan : new locations
  s=mlab.axes(extent=area,color=(0,0,0)) # axe des z positif vers le haut
  s=mlab.outline(extent=area,color=(0,0,0))
  s=mlab.title("cluster=%s, threshold=%s, nbmin=%s"%(i,threshold,nbmin),height=0.1,size=0.35,color=(0,0,0))


  if len(cluster[i]) < 20:
    from CZ_W_2_color import *
    for ind_I in range(len(cluster[i])):
      for ind_J in range(ind_I+1,len(cluster[i])):
        ev_I=cluster[i][ind_I]-1
        ev_J=cluster[i][ind_J]-1
        W_IJ=nbsta[ev_I,ev_J]
        if W_IJ >= nbmin:
          mlab.points3d(xini[ind_J],yini[ind_J],zini_ph[ind_J],scale_factor=0.1,color=(0,0,0))
          mlab.points3d(xini[ind_I],yini[ind_I],zini_ph[ind_I],scale_factor=0.1,color=(0,0,0))
          d=(xini[ind_J]-xini[ind_I],yini[ind_J]-yini[ind_I],zini_ph[ind_J]-zini_ph[ind_I])
          norm=np.sqrt(d[0]**2+d[1]**2+d[2]**2)
          s2=mlab.quiver3d(xini[ind_I],yini[ind_I],zini_ph[ind_I],d[0],d[1],d[2],color=tuple(CZ_W_2_color(W_IJ)),mode='2ddash',scale_factor=norm,scale_mode='scalar')

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
  opdict=read_header_from_file(loc_filename,opdict)


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
  corr_file=os.path.join(locdir,opdict['xcorr_corr'])
  cfile=BinaryFile(corr_file)
  coeff=cfile.read_binary_file()

  delay_file=os.path.join(locdir,opdict['xcorr_delay'])
  dfile=BinaryFile(delay_file)
  delay=dfile.read_binary_file()

  cluster_file=os.path.join(locdir,'cluster-%s-%s'%(str(threshold),str(nbmin)))
  clfile=BinaryFile(cluster_file)
  cluster=clfile.read_binary_file()

  # ----------------------------------------------------------------------------------------
  # Input parameters
  nb_iter=2
  len_cluster_min=2

  if dd_loc:
    new_loc_filename=os.path.join(locdir,'relocations.dat')
    new_loc_file=open(new_loc_filename,'w')
    write_header_options(new_loc_file,opdict)

  # ----------------------------------------------------------------------------------------
  for i in cluster.keys():
    print "CLUSTER %d:"%i,cluster[i],len(cluster[i])
    iter = 0
    N = len(cluster[i])

    # Hypocentral parameters to be changed
    x,y,z,z_ph,to = coord_cluster(cluster[i],locs)

    # Replace bad locations by the centroid coordinates
    centroid_x,centroid_y,centroid_z = np.mean(x), np.mean(y), np.mean(z)
    for ii in range(len(cluster[i])):
      if np.abs(x[ii]-centroid_x) > .75:
        x[ii]=centroid_x
      if np.abs(y[ii]-centroid_y) > .75:
        y[ii]=centroid_y
      if np.abs(z[ii]-centroid_z) > .75:
        z[ii]=centroid_z

    if N > len_cluster_min:

      # Theroretical traveltimes and arrival times
      t_th, arr_times = traveltimes(x,y,z,to,stations,time_grids)
      
      x,y,z,to = do_double_diff(x,y,z,to,stations,coeff,delay,cluster[i],threshold,t_th, arr_times)

      if verbose:
        from clustering import compute_nbsta
        nbsta=compute_nbsta(len(locs),coeff,threshold)
        plot_events(cluster,locs,stations,x,y,z,i,threshold,nbmin,area,nbsta)

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
          new_loc_file.write("Max = %.2f, %s - %.2f s + %.2f s, x= %.4f pm %.4f km, y= %.4f pm %.4f km, z= %.4f pm %.4f km\n"%(locs[j-1]['max_trig'],locs[j-1]['o_time'].isoformat(),locs[j-1]['o_err_left'], locs[j-1]['o_err_right'],locs[j-1]['x_mean'],locs[j-1]['x_sigma'],locs[j-1]['y_mean'],locs[j-1]['y_sigma'],locs[j-1]['z_mean'],locs[j-1]['z_sigma']))

  if dd_loc:
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
