import os, h5py, logging
import numpy as np
import scipy.integrate as si
import matplotlib.pyplot as plt
from matplotlib import mpl
from integrate4D import *
from filters import smooth
from copy import deepcopy

def plotLocationWaveforms(loc,start_time,dt,data_dict,grad_dict,stack_wfm,fig_dir):
    """
    Creates plot for located waveforms.  Assumes data and grad are ready
    for plotting.
    """

    otime=loc['o_time']
    otime_left=-loc['o_err_left']
    otime_right=loc['o_err_right']
#    iotime_left=np.int(np.round((otime - otime_left -start_time)/dt))
#    iotime_right=np.int(np.round((otime + otime_right -start_time)/dt))

    plot_filename=os.path.join(fig_dir,'loc_%s.pdf'%(otime.isoformat()))

    t=np.arange(len(stack_wfm))*dt - (otime - start_time)

    stations=data_dict.keys()
    stations.sort()
    n_traces=len(stations)+1

    plt.clf()
    fig=plt.figure()
    ax=fig.add_subplot(n_traces,2,1,title='Data')
    ax.set_axis_off()
    ax=fig.add_subplot(n_traces,2,2,title='Kurtosis Gradient')
    ax.set_axis_off()

    i=0
    for sta in stations :
        # plot the data in the first column
        ax=fig.add_subplot(n_traces,2,2*i+1)
        ax.set_axis_off()
        ax.plot(t,data_dict[sta],'b')
        ax.axvspan(otime_left,otime_right,facecolor='r', alpha=0.2)
        # add the station name
        pos=list(ax.get_position().bounds)
        fig.text(pos[0]-0.01, pos[1]+pos[3]/2., sta, fontsize = 10, 
                horizontalalignment = 'right', 
                verticalalignment = 'center')
        # plot the kurtosis gradient in the second column
        ax=fig.add_subplot(n_traces,2,2*i+2)
        ax.set_axis_off()
        ax.plot(t,grad_dict[sta],'b')
        ax.axvspan(otime_left,otime_right,facecolor='r', alpha=0.2)
        # add the maximum kurtosis value 
        pos=list(ax.get_position().bounds)
        fig.text(pos[0]+pos[2]+0.05 , pos[1], 
                '%.1f'%np.max(grad_dict[sta]), 
                fontsize = 10, horizontalalignment = 'right') 
        i=i+1

    #plot the stack under the kurtosis gradient only
    ax=fig.add_subplot(n_traces,2,2*n_traces,xlabel='time (s)')
    ax.plot(t,stack_wfm,'r')
    # put time axis only on plot
    ax.xaxis.set_ticks_position('none')
    ax.yaxis.set_ticks([])
    ax.spines['top'].set_color('none')
    ax.spines['left'].set_color('none')
    ax.spines['right'].set_color('none')
    ax.axvspan(otime_left,otime_right,facecolor='r', alpha=0.2)
    # add the maximum kurtosis value 
    pos=list(ax.get_position().bounds)
    fig.text(pos[0]+pos[2]+0.05 , pos[1], 
            '%.1f'%np.max(stack_wfm), 
            fontsize = 10, horizontalalignment = 'right') 

    #write the origin time under the data
    ax=fig.add_subplot(n_traces,2,2*n_traces-1)
    ax.set_axis_off()
    pos=list(ax.get_position().bounds)
    fig.text(pos[0]+pos[2]/2., pos[1]+pos[3]/2., otime.isoformat(), 
            fontsize = 12, horizontalalignment = 'center',
            verticalalignment = 'top')
    plt.savefig(plot_filename)
    plt.clf()


def plotLocationGrid(loc,grid_info,fig_dir):

  # set up plot using info from grid_info
  nx,ny,nz,nt = grid_info['grid_shape']
  dx,dy,dz,dt = grid_info['grid_spacing']
  x_orig,y_orig,z_orig = grid_info['grid_orig']
  stack_starttime=grid_info['start_time']

  # Take much of the information from the grid_info
  plot_info=deepcopy(grid_info)
  plot_info['o_time']=loc['o_time']

  # get location info
  o_time=loc['o_time']
  x_mean=loc['x_mean']
  y_mean=loc['y_mean']
  z_mean=loc['z_mean']
  o_err_left=loc['o_err_left']
  o_err_right=loc['o_err_right']
  x_sigma=loc['x_sigma']
  y_sigma=loc['y_sigma']
  z_sigma=loc['z_sigma']

  #get indexes correponding to location
  it_true=np.int(np.round((o_time-stack_starttime)/dt))
  ix_true=np.int(np.round((x_mean-x_orig)/dx))
  iy_true=np.int(np.round((y_mean-y_orig)/dy))
  iz_true=np.int(np.round((z_mean-z_orig)/dz))
  plot_info['true_indexes'] = (ix_true, iy_true, iz_true, it_true)
  plot_info['true_values'] = (x_mean, y_mean, z_mean, o_time-stack_starttime)

  # get indexes corresponding to location uncertainties
  # times are wrt stack_starttime
  t_left  = o_time - o_err_left  - stack_starttime
  t_right = o_time + o_err_right - stack_starttime
  # coordinates are absolute
  x_low   = x_mean - x_sigma 
  y_low   = y_mean - y_sigma
  z_low   = z_mean - z_sigma 
  x_high  = x_mean + x_sigma
  y_high  = y_mean + y_sigma
  z_high  = z_mean + z_sigma

  plot_info['t_err'] = (t_left, t_right)
  plot_info['x_err'] = (x_low, x_high)
  plot_info['y_err'] = (y_low, y_high)
  plot_info['z_err'] = (z_low, z_high)

  plotDiracTest(plot_info,fig_dir)

def plotDiracTest(test_info,fig_dir):

  # set up plot using info from test_info
  nx,ny,nz,nt = test_info['grid_shape']
  dx,dy,dz,dt = test_info['grid_spacing']
  x_orig,y_orig,z_orig = test_info['grid_orig']
  ix_true, iy_true, iz_true, it_true= test_info['true_indexes']  
  if test_info.has_key('true_values'): 
      x_true, y_true, z_true, t_true= test_info['true_values']  
  stack_start_time=test_info['start_time']
  grid_filename=test_info['dat_file']
  stack_filename=test_info['stack_file']
  if test_info.has_key('o_time'):
    fig_filename = os.path.join(fig_dir,"grid_%s.pdf" % 
            test_info['o_time'].isoformat())
  else:
    fig_filename = os.path.join(fig_dir,"%s.pdf" % 
            os.path.splitext(os.path.basename(grid_filename))[0])

  # read the stack file
  f=h5py.File(grid_filename,'r')
  stack_grid=f['stack_grid']

  stack_3D=stack_grid[:,it_true].reshape(nx,ny,nz)


  # cut through the true location at the true time 
  xy_cut=stack_3D[:,:,iz_true] 
  xz_cut=stack_3D[:,iy_true,:]
  yz_cut=stack_3D[ix_true,:,:]

  # extract the max stacks
  f_stack=h5py.File(stack_filename,'r')
  # if have a smoothed version, use it for the plots
  if 'max_val_smooth' in f_stack:
    max_val=f_stack['max_val_smooth']
  else:
    max_val=f_stack['max_val']
  max_x=f_stack['max_x']
  max_y=f_stack['max_y']
  max_z=f_stack['max_z']
  

  # set up the 4 axes
  x=np.arange(nx)*dx
  y=np.arange(ny)*dy
  z=(np.arange(nz)*dz+z_orig)*(-1)
  # setup of t-axis depends on type of stack_start_time
  if type(stack_start_time)==float:
    t=np.arange(nt)*dt+stack_start_time
  else:
    t=np.arange(nt)*dt

  # set time origin to o_time
  #t=t-it_true*dt

  # do plot
  plt.clf()
  fig=plt.figure()


  if test_info.has_key('true_values') and test_info.has_key('o_time'): 
    fig.suptitle('%s   x = %.2fkm  y = %.2fkm  z = %.2fkm'%(test_info['o_time'].isoformat(), x_true, y_true, z_true))

  # plot xy plane
  p=plt.subplot(2,2,1)
  pos=list(p.get_position().bounds)
  fig.text(pos[0]-0.08,pos[1]+pos[3], '(a)', fontsize=12)
  plt.imshow(xy_cut.T,origin='lower',interpolation='none',extent=[np.min(x),np.max(x),np.min(y),np.max(y)])
 # if test_info.has_key('x_err'):
 #     x_low,x_high=test_info['x_err']
 #     plt.vlines(x_low,np.min(y),np.max(y),'w',linewidth=1)
 #     plt.vlines(x_high,np.min(y),np.max(y),'w',linewidth=1)
  p.tick_params(labelsize=10)
  p.xaxis.set_ticks_position('top')
  plt.xlabel('x (km wrt ref)',size=10)
  plt.ylabel('y (km wrt ref)',size=10)
  #plt.title('XY plane')

  #plot xz plane
  p=plt.subplot(4,2,5)
  pos=list(p.get_position().bounds)
  fig.text(pos[0]-0.08,pos[1]+pos[3], '(b)', fontsize=12)
  plt.imshow(xz_cut.T,origin='upper',interpolation='none',extent=[np.min(x),np.max(x),np.min(z),np.max(z)])
  p.tick_params(labelsize=10)
  p.xaxis.set_ticks_position('top')
  p.xaxis.set_ticks(())
  #plt.xlabel('x (km wrt ref)')
  plt.ylabel('z (km up)',size=10)
  #plt.title('XZ plane')

  # plot yz plane
  p=plt.subplot(4,2,7)
  pos=list(p.get_position().bounds)
  fig.text(pos[0]-0.08,pos[1]+pos[3], '(c)', fontsize=12)
  plt.imshow(yz_cut.T,origin='upper',interpolation='none',extent=[np.min(y),np.max(y),np.min(z),np.max(z)])
  p.xaxis.set_ticks_position('bottom')
  p.tick_params(labelsize=10)
  plt.xlabel('y (km wrt ref)',size=10)
  plt.ylabel('z (km up)',size=10)
  #plt.title('YZ plane')

  # choose portion of time series to plot
  if test_info.has_key('true_values'):
    llim = t_true-2.0
    rlim = t_true+2.0
  else:
    llim = t[it_true]-2.0
    rlim = t[it_true]+2.0

  illim = int((llim-t[0])/dt)
  irlim = int((rlim-t[0])/dt)

  # plot max value
  p=plt.subplot(4,2,2, frameon=False)
  pos=list(p.get_position().bounds)
  fig.text(pos[0],pos[1]+pos[3], '(d)', fontsize=12)
  p.tick_params(labelsize=10)
  plt.plot(t,max_val)
  #plt.xlabel('t (s)')
  p.xaxis.set_ticks_position('none')
  p.xaxis.set_ticks(())
  plt.ylabel('Stack max',size=10)
  p.yaxis.set_ticks_position('right')
  #plt.title('Maximum of stack')
  p.set_xlim(llim,rlim)
  p.set_ylim(0,max(max_val))
  #p.xaxis.set_ticks_position('bottom')
  #p.xaxis.set_ticks()
  plt.vlines(t[it_true],0,max(max_val),'r',linewidth=2)
  if test_info.has_key('t_err'):
      t_left,t_right=test_info['t_err']
      plt.axvspan(t_left,t_right,facecolor='r', alpha=0.2)

  # put the origin back in for the last plots
  x=np.arange(nx)*dx+x_orig
  y=np.arange(ny)*dy+y_orig
  z=np.arange(nz)*dz+z_orig

  # plot max x
  p=plt.subplot(4,2,4, frameon=False)
  pos=list(p.get_position().bounds)
  fig.text(pos[0],pos[1]+pos[3], '(e)', fontsize=12)
  p.tick_params(labelsize=10)
  #plt.plot(t,max_x,'b.',clip_on=False)
  plt.scatter(t[illim:irlim],max_x[illim:irlim],s=40, c=max_val[illim:irlim],marker='.',linewidths=(0,),clip_on=False)
  #plt.xticks([llim,t[it_true],rlim])
  #plt.xlabel('t (s)')
  p.xaxis.set_ticks_position('none')
  p.xaxis.set_ticks(())
  plt.ylabel('x (km)',size=10)
  p.yaxis.set_ticks_position('right')
  #plt.title('x at maximum')
  p.set_xlim(llim,rlim)
  if test_info.has_key('true_values'):
    if not test_info.has_key('t_err'):
      plt.hlines(x_true,llim,rlim,'r',linewidth=2)
      plt.vlines(t_true,min(max_x),max(max_x),'r',linewidth=2)
  else:
    plt.hlines(x[ix_true],llim,rlim,'r',linewidth=2)
    plt.vlines(t[it_true],min(max_x),max(max_x),'r',linewidth=2)
  if test_info.has_key('x_err'):
      x_low,x_high=test_info['x_err']
      plt.axhspan(x_low,x_high,facecolor='r', alpha=0.2)
  if test_info.has_key('t_err'):
      t_left,t_right=test_info['t_err']
      plt.axvspan(t_left,t_right,facecolor='r', alpha=0.2)

  # plot max y
  p=plt.subplot(4,2,6, frameon=False)
  pos=list(p.get_position().bounds)
  fig.text(pos[0],pos[1]+pos[3], '(f)', fontsize=12)
  p.tick_params(labelsize=10)
  #plt.plot(t,max_y,'b.')
  plt.scatter(t[illim:irlim],max_y[illim:irlim],s=40, c=max_val[illim:irlim],marker='.',linewidths=(0,),clip_on=False)
  #plt.xticks([llim,t[it_true],rlim])
  #plt.xlabel('t (s)')
  p.xaxis.set_ticks_position('none')
  p.xaxis.set_ticks(())
  plt.ylabel('y (km)',size=10)
  p.yaxis.set_ticks_position('right')
  #plt.title('y at maximum')
  p.set_xlim(llim,rlim)
  if test_info.has_key('true_values'):
    if not test_info.has_key('t_err'):
      plt.hlines(y_true,llim,rlim,'r',linewidth=2)
      plt.vlines(t_true,min(max_y),max(max_y),'r',linewidth=2)
  else:
    plt.hlines(y[iy_true],llim,rlim,'r',linewidth=2)
    plt.vlines(t[it_true],min(max_y),max(max_y),'r',linewidth=2)
  if test_info.has_key('y_err'):
      y_low,y_high=test_info['y_err']
      plt.axhspan(y_low,y_high,facecolor='r', alpha=0.2)
  if test_info.has_key('t_err'):
      t_left,t_right=test_info['t_err']
      plt.axvspan(t_left,t_right,facecolor='r', alpha=0.2)

  # plot max z
  p=plt.subplot(4,2,8, frameon=False)
  pos=list(p.get_position().bounds)
  fig.text(pos[0],pos[1]+pos[3], '(g)', fontsize=12)
  p.tick_params(labelsize=10)
  #plt.plot(t,max_z,'b.')
  plt.scatter(t[illim:irlim],max_z[illim:irlim],s=40, c=max_val[illim:irlim],marker='.',linewidths=(0,),clip_on=False)
  #plt.xticks([llim,t[it_true],rlim])
  plt.xlabel('Time (s)',size=10)
  p.xaxis.set_ticks_position('bottom')
  plt.ylabel('z (km down)',size=10)
  p.yaxis.set_ticks_position('right')
  #plt.title('z at maximum')
  p.set_xlim(llim,rlim)
  if test_info.has_key('true_values'):
    if not test_info.has_key('t_err'):
      plt.hlines(z_true,llim,rlim,'r',linewidth=2)
      plt.vlines(t_true,min(max_z),max(max_z),'r',linewidth=2)
  else:
    plt.hlines(z[iz_true],llim,rlim,'r',linewidth=2)
    plt.vlines(t[it_true],min(max_z),max(max_z),'r',linewidth=2)
  if test_info.has_key('z_err'):
      z_low,z_high=test_info['z_err']
      plt.axhspan(z_low,z_high,facecolor='r', alpha=0.2)
  if test_info.has_key('t_err'):
      t_left,t_right=test_info['t_err']
      plt.axvspan(t_left,t_right,facecolor='r', alpha=0.2)

  # add independent colorbar
  ax1 = fig.add_axes([0.40, 0.03, 0.2, 0.015])
  ax1.tick_params(labelsize=8)
  ax1.xaxis.set_ticks_position('bottom')
  cmap = mpl.cm.jet
  norm = mpl.colors.Normalize(vmin=min(max_val), vmax=max(max_val))
  cb = mpl.colorbar.ColorbarBase(ax1, cmap=cmap, norm=norm, orientation='horizontal')
  pos=list(ax1.get_position().bounds)
  fig.text(pos[0]+pos[2]/2.,pos[1]+pos[3]+0.01, 'Stack max', fontsize=8, horizontalalignment='center', verticalalignment='bottom')

  #plt.tight_layout()
  plt.savefig(fig_filename)

  f.close()
  f_stack.close()

def plot_test(curve_tuple,axes_tuple,filename_base):
  stack_x, stack_y, stack_z, stack_t = curve_tuple
  x,y,z,t = axes_tuple

  stack_x = stack_x / compute_integral1D(stack_x,x)
  stack_y = stack_y / compute_integral1D(stack_y,y)
  stack_z = stack_z / compute_integral1D(stack_z,z)
  stack_t = stack_t / compute_integral1D(stack_t,t)

  # do 1D plots
  plt.clf()
  plt.plot(x,stack_x)
  plt.xlabel('x (km)')
  plt.ylabel('p(x)')
  plt.title('Marginal probability density over x (at maximum)')
  plt.savefig(filename_base+'_test_stack_x.pdf')

  plt.clf()
  plt.plot(y,stack_y)
  plt.xlabel('y (km)')
  plt.ylabel('p(y)')
  plt.title('Marginal probability density over y (at maximum)')
  plt.savefig(filename_base+'_test_stack_y.pdf')

  plt.clf()
  plt.plot(z,stack_z)
  plt.xlabel('z (km)')
  plt.ylabel('p(z)')
  plt.title('Marginal probability density over z (at maximum)')
  plt.savefig(filename_base+'_test_stack_z.pdf')

  plt.clf()
  plt.plot(t,stack_t)
  plt.xlabel('t (km)')
  plt.ylabel('p(t)')
  plt.title('Marginal probability density over t (at maximum)')
  plt.savefig(filename_base+'_test_stack_t.pdf')

def plot_probloc_mpl3D(grid_dict,x_list,base_filename):

  #print grid_dict.keys()

  x0,x1,x2 = x_list

  # create filenames for plots
  files_dict={}
  for key in grid_dict.keys():
    files_dict[key] = base_filename + '_' + key + '.pdf'
  
  logging.debug('PLOTTING : Type x0 = %s'%x0.dtype)
  logging.debug('PLOTTING : Type prob_x0 = %s'%grid_dict['prob_x0'].dtype)
  
  # do 1D plots
  plt.clf()
  plt.plot(x0,grid_dict['prob_x0'])
  plt.xlabel('x (km)')
  plt.ylabel('p(x)')
  plt.title('Marginal probability density over x')
  plt.savefig(files_dict['prob_x0'])


  plt.clf()
  plt.plot(x1,grid_dict['prob_x1'])
  plt.xlabel('y (km)')
  plt.ylabel('p(y)')
  plt.title('Marginal probability density over y')
  plt.savefig(files_dict['prob_x1'])

  plt.clf()
  plt.plot(x2,grid_dict['prob_x2'])
  plt.xlabel('z (km)')
  plt.ylabel('p(z)')
  plt.title('Marginal probability density over z')
  plt.savefig(files_dict['prob_x2'])


  # 2D plots
  plt.clf()
  p=plt.contourf(x0,x1,grid_dict['prob_x0_x1'].T)
  plt.colorbar(p)
  plt.xlabel('x (km)')
  plt.ylabel('y (km)')
  plt.title('Marginal probability density over x and y')
  plt.savefig(files_dict['prob_x0_x1'])
 
  plt.clf()
  p=plt.contourf(x0,x2,grid_dict['prob_x0_x2'].T)
  plt.colorbar(p)
  plt.xlabel('x (km)')
  plt.ylabel('z (km)')
  plt.title('Marginal probability density over x and z')
  plt.savefig(files_dict['prob_x0_x2'])

  plt.clf()
  p=plt.contourf(x1,x2,grid_dict['prob_x1_x2'].T)
  plt.colorbar(p)
  plt.xlabel('y (km)')
  plt.ylabel('z (km)')
  plt.title('Marginal probability density over y and z')
  plt.savefig(files_dict['prob_x1_x2'])
 

def plot_probloc_mpl(grid_dict,x_list,base_filename):

  #print grid_dict.keys()

  x0,x1,x2,x3 = x_list

  # create filenames for plots
  files_dict={}
  for key in grid_dict.keys():
    files_dict[key] = base_filename + '_' + key + '.pdf'
  
  logging.debug('PLOTTING : Type x0 = %s'%x0.dtype)
  logging.debug('PLOTTING : Type prob_x0 = %s'%grid_dict['prob_x0'].dtype)
  
  # do 1D plots
  plt.clf()
  plt.plot(x0,grid_dict['prob_x0'])
  plt.xlabel('x (km)')
  plt.ylabel('p(x)')
  plt.title('Marginal probability density over x')
  plt.savefig(files_dict['prob_x0'])

  logging.debug('Integral over x0 = %.3f'%si.trapz(grid_dict['prob_x0'],x=x0))

  plt.clf()
  plt.plot(x1,grid_dict['prob_x1'])
  plt.xlabel('y (km)')
  plt.ylabel('p(y)')
  plt.title('Marginal probability density over y')
  plt.savefig(files_dict['prob_x1'])

  plt.clf()
  plt.plot(x2,grid_dict['prob_x2'])
  plt.xlabel('z (km)')
  plt.ylabel('p(z)')
  plt.title('Marginal probability density over z')
  plt.savefig(files_dict['prob_x2'])

  plt.clf()
  plt.plot(x3,grid_dict['prob_x3'])
  plt.xlabel('time (s)')
  plt.ylabel('p(t)')
  plt.title('Marginal probability density over time')
  plt.savefig(files_dict['prob_x3'])

  # 2D plots
  plt.clf()
  p=plt.contourf(x0,x3,grid_dict['prob_x0_x3'].T)
  plt.colorbar(p)
  plt.xlabel('x (km)')
  plt.ylabel('time (s)')
  plt.title('Marginal probability density over x and time')
  plt.savefig(files_dict['prob_x0_x3'])
  
  plt.clf()
  p=plt.contourf(x1,x3,grid_dict['prob_x1_x3'].T)
  plt.colorbar(p)
  plt.xlabel('y (km)')
  plt.ylabel('time (s)')
  plt.title('Marginal probability density over y and time')
  plt.savefig(files_dict['prob_x1_x3'])

  plt.clf()
  p=plt.contourf(x2,x3,grid_dict['prob_x2_x3'].T)
  plt.colorbar(p)
  plt.xlabel('z (km)')
  plt.ylabel('time (s)')
  plt.title('Marginal probability density over z and time')
  plt.savefig(files_dict['prob_x2_x3'])
 
  plt.clf()
  p=plt.contourf(x0,x1,grid_dict['prob_x0_x1'].T)
  plt.colorbar(p)
  plt.xlabel('x (km)')
  plt.ylabel('y (km)')
  plt.title('Marginal probability density over x and y')
  plt.savefig(files_dict['prob_x0_x1'])
 
  plt.clf()
  p=plt.contourf(x0,x2,grid_dict['prob_x0_x2'].T)
  plt.colorbar(p)
  plt.xlabel('x (km)')
  plt.ylabel('z (km)')
  plt.title('Marginal probability density over x and z')
  plt.savefig(files_dict['prob_x0_x2'])

  plt.clf()
  p=plt.contourf(x1,x2,grid_dict['prob_x1_x2'].T)
  plt.colorbar(p)
  plt.xlabel('y (km)')
  plt.ylabel('z (km)')
  plt.title('Marginal probability density over y and z')
  plt.savefig(files_dict['prob_x1_x2'])
 
