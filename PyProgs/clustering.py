#!/usr/bin/env python
# encoding: utf-8

import os, sys, optparse, glob
import matplotlib.pyplot as plt
import numpy as np
import cProfile
from CZ_Clust_2_color import *
from CZ_W_2_color import *
import cPickle
from OP_waveforms import *
import logging
from correlation import BinaryFile
from locations_trigger import read_locs_from_file
from NllGridLib import read_stations_file
# ----------------------------------------------------------------------------------------
class Graph(object):
  def __init__(self):
    self.flag=[]
    self.cluster_index=[]
    self.voisins=[]
    self.nb_voisins=[]

  def set_flag(self,value):
    self.flag.append(value)
  def set_cluster_index(self,value):
    self.cluster_index.append(value)
  def set_voisins(self,value):
    self.voisins.append(value)
# ----------------------------------------------------------------------------------------
# Depth First Search algorithm
def CZ_DFS(GRAPH,sommet_first,cluster_ind):
  GRAPH.flag[sommet_first]=1
  GRAPH.cluster_index[sommet_first]=cluster_ind
  for ind_sommet_fils in range(len(GRAPH.voisins[sommet_first])):
    sommet_fils=GRAPH.voisins[sommet_first][ind_sommet_fils]
    if GRAPH.flag[sommet_fils]==0:
      GRAPH,cluster_ind=CZ_DFS(GRAPH,sommet_fils,cluster_ind)
  return GRAPH,cluster_ind
# -----------------------------------------------------------------------------------------
def waveval(stack_time,t_before,t_after,dt,tdeb):
  tstart=stack_time-t_before-tdeb
  tend=stack_time+t_after-tdeb
  i_start=int(round(tstart*1./dt))
  i_end=int(round(tend*1./dt))
  return i_start,i_end
# --------------------------------------------------------------------------------------------------------------------
def plot_traces(CLUSTER, delay_file, coeff, locs, stations, datadir, data_files, threshold):
  # Read the file containing the time delays
  a=BinaryFile(delay_file)
  delay=a.read_binary_file()

  t_before=0.5
  t_after=6.0

  rg_x=[362,370]
  rg_y=[7647,7653]
  rg_z=[-4.0,2.5]

  tr={}
  for data_file in data_files:
    wf=Waveform()
    wf.read_from_file(data_file)
    tr[wf.station]=wf.values
    dt=wf.delta
    tdeb=wf.starttime

  list_name=sorted(tr)

  for i in range(1,len(CLUSTER)+1): # cluster index
    for j in range(len(CLUSTER[i])): # first event index
      e1=CLUSTER[i][j]
      for k in range(j+1,len(CLUSTER[i])): # second event index
        e2=CLUSTER[i][k]
        co=0
        fig = plt.figure()
        fig.set_facecolor('white')
        for l in range(len(list_name)):
          name=list_name[l]
          if delay[name][e1-1][e2-1]!='NaN':
            stack_time_1=locs[e1-1]['o_time']
            i_start_1, i_end_1=waveval(stack_time_1,t_before,t_after,dt,tdeb)
            val1=tr[name][i_start_1-1:i_end_1]
            stack_time_2=locs[e2-1]['o_time']
            i_start_2, i_end_2=waveval(stack_time_2-delay[name][e1-1][e2-1],t_before,t_after,dt,tdeb)
            val2=tr[name][i_start_2-1:i_end_2]
            t=np.linspace(0,t_after+t_before,(t_after+t_before)/dt+1)
            ax=fig.add_subplot(len(list_name),1,l+1)
            ax.set_axis_off()
            ax.plot(t,val1/max(val1),'k')
            ax.plot(t,val2/max(val2),'y--')
            c='k'
            if coeff[name][e1-1][e2-1]>=threshold:
              co=co+1
              c='r'
            ax.text(0.2,0.5,"%s, %s, %s"%(name,str(coeff[name][e1-1][e2-1]),str(delay[name][e1-1][e2-1])),color=c)
        fig.suptitle("Cluster : %s ; Event pair : (%s,%s) ; %d"%(str(i),str(e1),str(e2),co))
        plt.show()

#        # Plot location
#        fig = plt.figure()
#        x1=locs[e1-1]['x_mean']
#        y1=locs[e1-1]['y_mean']
#        z1=-locs[e1-1]['z_mean']
#        x2=locs[e2-1]['x_mean']
#        y2=locs[e2-1]['y_mean']
#        z2=-locs[e2-1]['z_mean']
#        ax1=fig.add_subplot(221,xlabel='x',ylabel='y')
#        ax1.plot(x1,y1,'bo',x2,y2,'ro')
#        ax1.axis(rg_x+rg_y)
#        ax1.text(367.5,7652.5,"dx=%.03f"%np.abs(x1-x2))
#        ax2=fig.add_subplot(222,xlabel='x',ylabel='z')
#        ax2.plot(x1,z1,'bo',x2,z2,'ro')
#        ax2.axis(rg_x+rg_z)
#        ax2.text(367.5,2,"dz=%.03f"%np.abs(z1-z2))
#        ax3=fig.add_subplot(223,xlabel='y',ylabel='z')
#        ax3.plot(y1,z1,'bo',y2,z2,'ro')
#        ax3.axis(rg_y+rg_z)
#        ax3.text(7651,2,"dy=%.03f"%np.abs(y1-y2))
#        #plt.show()
#
#        fig = plt.figure()
#        fig.set_facecolor('white')
#        for l in range(len(list_name)):
#          name=list_name[l]
#          if delay[name][e1-1][e2-1]!='NaN' and coeff[name][e1-1][e2-1]>0.8:
#            stack_time_1=locs[e1-1]['o_time']
#            i_start_1, i_end_1=waveval(stack_time_1,t_before,t_after,dt,tdeb)
#            val1=tr[name][i_start_1-1:i_end_1]
#            stack_time_2=locs[e2-1]['o_time']
#            i_start_2, i_end_2=waveval(stack_time_2-delay[name][e1-1][e2-1],t_before,t_after,dt,tdeb)
#            val2=tr[name][i_start_2-1:i_end_2]
#            t=np.linspace(0,t_after+t_before,(t_after+t_before)/dt+1)
#            plt.plot(t,val1/max(val1))
#            plt.plot(t,val2/max(val2),'r')
#            plt.title("%s - coeff: %s"%(name,coeff[name][e1-1][e2-1]))
#            #plt.show()
# --------------------------------------------------------------
def compute_nbsta(event,coeff,threshold):
  # Compute the number of stations where the correlation value is >= threshold for every event pair
  nbsta=[]
  for i in xrange(event):
    liste=[]
    for k in xrange(i):
      liste.append(0)
    for j in xrange(i,event):
      c=0
      if i!=j:
        for name in sorted(coeff):
          if coeff[name] and coeff[name][i][j] >= threshold and coeff[name][i][j] != 'NaN':
            c=c+1
        liste.append(c)
      else:
        liste.append(0)
    nbsta.append(liste)

  nbsta=np.matrix(nbsta)
  return nbsta

# -----------------------------------------------------------------------------------------
def do_clustering(event,nbsta,nbmin):
  # CODE CHRISTOPHE - CLUSTERING : DEPTH FIRST SEARCH
  voisins_du_sommet_I__horiz, voisins_du_sommet_I__verti, voisins_du_sommet_I=[],[],[]
  GRAPH=Graph()
  sommets=[]
 
  for I in range(event):
    voisins_du_sommet_I__verti=(np.where(nbsta[:,I]>=nbmin)[0]).tolist()[0]
    voisins_du_sommet_I__horiz=(np.where(nbsta[I,:]>=nbmin)[1]).tolist()[0]
    GRAPH.set_voisins(voisins_du_sommet_I__verti+voisins_du_sommet_I__horiz)
    GRAPH.set_flag(0)
    GRAPH.set_cluster_index(0)
    if voisins_du_sommet_I__verti+voisins_du_sommet_I__horiz:
      sommets.append(I)

  if sommets:
    NB_MAX_VOISINS=0
    for ind_sommet in range(len(sommets)):
      l=len(GRAPH.voisins[sommets[ind_sommet]])
      if l > NB_MAX_VOISINS:
        sommet_first=sommets[ind_sommet]
        NB_MAX_VOISINS=l

    ind_sommet_first=0
    cluster_ind=0
    CLUSTER={}
    while 1:
      event_index_flagged=[]
      event_index_non_flagged_with_neighbours=[]
      ind_sommet_first=ind_sommet_first+1
      cluster_ind=cluster_ind+1
      GRAPH,cluster_ind=CZ_DFS(GRAPH,sommet_first,cluster_ind)
      for k in range(len(GRAPH.voisins)):
        if GRAPH.flag[k] == 1 and GRAPH.cluster_index[k] == ind_sommet_first:
          event_index_flagged.append(k)
        elif GRAPH.flag[k] == 0 and len(GRAPH.voisins[k]) != 0:
          event_index_non_flagged_with_neighbours.append(k)

      # add 1 to each event number as the first one is number one (and not zero)
      CLUSTER[cluster_ind]=list(event_index_flagged+np.ones(len(event_index_flagged),dtype=np.int))
      if len(event_index_non_flagged_with_neighbours) > 1:
        sommet_first=event_index_non_flagged_with_neighbours[0]
      else:
        break

    return CLUSTER

  else:
    return {}
# -----------------------------------------------------------------------------------------
def plot_graphs(locs,stations,nbsta,CLUSTER,nbmin,threshold):
  from mayavi import mlab

  # Event coordinates
  stack_x,stack_y,stack_z=[],[],[]
  for loc in locs:
    stack_x.append(loc['x_mean'])
    stack_y.append(loc['y_mean'])
    stack_z.append(-loc['z_mean'])

  # Extract coordinates
  xsta,ysta,zsta=[],[],[]
  for sta in sorted(stations):
    xsta.append(stations[sta]['x'])
    ysta.append(stations[sta]['y'])
    zsta.append(stations[sta]['elev'])

  # 3D PLOT USING MAYAVI
  logging.info("Plotting...")
  s1=mlab.figure(1,bgcolor=(1,1,1),fgcolor=(0,0,0),size=(1000,900))
  s1=mlab.points3d(xsta,ysta,zsta,color=(1,0,0),scale_factor=0.05,mode='cube')
  s1=mlab.axes(extent=[362,370,7647,7653,-0.5,2.5],color=(0,0,0))
  s1=mlab.outline(extent=[362,370,7647,7653,-0.5,2.5],color=(0,0,0))
  s1=mlab.points3d(stack_x,stack_y,stack_z,scale_factor=0.1,color=(0.8,0.8,0.8))
  s1=mlab.title("threshold=%s, nbmin=%s"%(threshold,nbmin),height=0.1,size=0.35,color=(0,0,0))
  for i_ev in range(len(nbsta)):
    for i_c in range(1,len(CLUSTER)+1):
      if i_ev+1 in CLUSTER[i_c]:
        s1=mlab.points3d(stack_x[i_ev],stack_y[i_ev],stack_z[i_ev],scale_factor=0.1,color=tuple(CZ_Clust_2_color(100*(len(CLUSTER)-i_c)/len(CLUSTER))))
        s1=mlab.text3d(stack_x[i_ev],stack_y[i_ev],stack_z[i_ev],str(i_c),color=(0,0,0),scale=0.1)
  logging.info("Done!")
   
  logging.info("Plotting...")
  s2=mlab.figure(2,bgcolor=(1,1,1),fgcolor=(0,0,0),size=(1000,900))
  mlab.points3d(xsta,ysta,zsta,color=(1,0,0),scale_factor=0.05,mode='cube')
  mlab.axes(extent=[362,370,7647,7653,-0.5,2.5],color=(0,0,0))
  mlab.outline(extent=[362,370,7647,7653,-0.5,2.5],color=(0,0,0))
  mlab.points3d(stack_x,stack_y,stack_z,scale_factor=0.1,color=(0.8,0.8,0.8))
  mlab.title("threshold=%s, nbmin=%s"%(threshold,nbmin),height=0.1,size=0.35,color=(0,0,0))
  print nbsta
  for ind_I in range(len(nbsta)):
    for ind_J in range(ind_I+1,len(nbsta)):
      W_IJ=nbsta[ind_I,ind_J]
      if W_IJ >= nbmin:
        mlab.points3d(stack_x[ind_J],stack_y[ind_J],stack_z[ind_J],scale_factor=0.1,color=(0,0,0))
        mlab.points3d(stack_x[ind_I],stack_y[ind_I],stack_z[ind_I],scale_factor=0.1,color=(0,0,0))
        d=(stack_x[ind_J]-stack_x[ind_I],stack_y[ind_J]-stack_y[ind_I],stack_z[ind_J]-stack_z[ind_I])
        norm=np.sqrt(d[0]**2+d[1]**2+d[2]**2)
        s2=mlab.quiver3d(stack_x[ind_I],stack_y[ind_I],stack_z[ind_I],d[0],d[1],d[2],color=tuple(CZ_W_2_color(W_IJ)),mode='2ddash',scale_factor=norm,scale_mode='scalar')
  #mlab.colorbar(s2)
  logging.info("Done!")
  mlab.show()


def do_clustering_setup_and_run(opdict):

  base_path=opdict['base_path']
  verbose=opdict['verbose']

  # stations
  stations_filename=os.path.join(base_path,'lib',opdict['stations'])

  # output directory
  output_dir=os.path.join(base_path,'out',opdict['outdir'])

  # data
  data_dir=os.path.join(base_path,'data',opdict['datadir'])
  data_glob=opdict['dataglob']
  data_files=glob.glob(os.path.join(data_dir,data_glob))
  data_files.sort()

  # location file
  locdir=os.path.join(base_path,'out',opdict['outdir'],'loc')
  loc_filename=os.path.join(locdir,'locations.dat')

  # file containing correlation values
  coeff_file=os.path.join(locdir,opdict['xcorr_corr'])
  # Read correlation values
  b=BinaryFile(coeff_file)
  coeff=b.read_binary_file()

  # file containing time delays
  delay_file=os.path.join(locdir,opdict['xcorr_delay'])

  # INPUT PARAMETERS
  nbmin=int(opdict['nbsta'])
  if nbmin > len(coeff.keys()):
    raise Error('the minimum number of stations cannot be > to the number of stations !!')
  event=len(coeff.values()[0])
  tplot=float(opdict['clus']) # threshold for which we save and plot 
  cluster_file="%s/cluster-%s-%s"%(locdir,str(tplot),str(nbmin))

  corr=[opdict['clus']]
  #corr=np.arange(0,1.1,0.1)
  for threshold in corr:
    threshold=float(threshold)
    nbsta=compute_nbsta(event,coeff,threshold)

    CLUSTER = do_clustering(event,nbsta,nbmin)

    if threshold == tplot:

      print "----------------------------------------------"
      print "THRESHOLD : ",threshold," # STATIONS : ",nbmin
      print "# CLUSTERS : ",len(CLUSTER)
      print CLUSTER

      c=BinaryFile(cluster_file)
      c.write_binary_file(CLUSTER)
      print "Written in %s"%cluster_file

      if verbose: # PLOT
        # Read location file
        locs=read_locs_from_file(loc_filename)
        # Read station file 
        stations=read_stations_file(stations_filename)

        # Look at the waveforms 
        #plot_traces(CLUSTER, delay_file, coeff, locs, stations, data_dir, data_files, threshold)

        # Plot graphs
        plot_graphs(locs,stations,nbsta,CLUSTER,nbmin,threshold)

# -----------------------------------------------------------------------------------------
if __name__ == '__main__':

  from options import WavelocOptions
  logging.basicConfig(level=logging.INFO, format='%(levelname)s : %(asctime)s : %(message)s')

  wo = WavelocOptions()
  args=wo.p.parse_args()

  wo.set_all_arguments(args)
  wo.verify_cluster_options()

  do_clustering_setup_and_run(wo.opdict)
