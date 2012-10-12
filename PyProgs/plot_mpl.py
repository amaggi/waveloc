import os, h5py, logging
import numpy as np
import scipy.integrate as si
import matplotlib.pyplot as plt
from integrate4D import *
from filters import smooth
from copy import deepcopy

def plotLocationWaveforms(loc,data_dict,grad_dict,stack_wfm,fig_dir):
    """
    Creates plot for located waveforms.  Assumes data and grad are ready
    for plotting.
    """

    otime=loc['o_time']

    plot_filename=os.path.join(fig_dir,'loc_%s.pdf'%(otime.isoformat()))

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
        ax.plot(data_dict[sta],'b')
        # add the station name
        pos=list(ax.get_position().bounds)
        fig.text(pos[0]-0.01, pos[1], sta, fontsize = 10, 
                horizontalalignment = 'right')
        # plot the kurtosis gradient in the second column
        ax=fig.add_subplot(n_traces,2,2*i+2)
        ax.set_axis_off()
        ax.plot(grad_dict[sta],'b')
        i=i+1

    #plot the stack under the kurtosis gradient only
    ax=fig.add_subplot(n_traces,2,2*n_traces)
    ax.set_axis_off()
    ax.plot(stack_wfm,'r')
    pos=list(ax.get_position().bounds)

    #write the origin time under the data
    ax=fig.add_subplot(n_traces,2,2*n_traces-1)
    ax.set_axis_off()
    fig.text(pos[0]-0.05, pos[1], otime.isoformat(), fontsize = 10, 
            horizontalalignment = 'right')
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

  # get indexes corresponding to location uncertainties
  it_left=np.int(np.round((o_time-o_err_left-stack_starttime)/dt))
  it_right=np.int(np.round((o_time+o_err_right-stack_starttime)/dt))
  ix_low=np.int(np.round((x_mean-x_sigma-x_orig)/dx))
  iy_low=np.int(np.round((y_mean-y_sigma-y_orig)/dy))
  iz_low=np.int(np.round((z_mean-z_sigma-z_orig)/dz))
  ix_high=np.int(np.round((x_mean+x_sigma-x_orig)/dx))
  iy_high=np.int(np.round((y_mean+y_sigma-y_orig)/dy))
  iz_high=np.int(np.round((z_mean+z_sigma-z_orig)/dz))

  plot_info['t_err'] = (it_left, it_left)
  plot_info['x_err'] = (ix_low, ix_high)
  plot_info['y_err'] = (iy_low, iy_high)
  plot_info['z_err'] = (iz_low, iz_high)

  plotDiracTest(plot_info,fig_dir)

def plotDiracTest(test_info,fig_dir):

  # set up plot using info from test_info
  nx,ny,nz,nt = test_info['grid_shape']
  dx,dy,dz,dt = test_info['grid_spacing']
  x_orig,y_orig,z_orig = test_info['grid_orig']
  ix_true, iy_true, iz_true, it_true= test_info['true_indexes']  
  stack_start_time=test_info['start_time']
  grid_filename=test_info['dat_file']
  stack_filename=test_info['stack_file']
  fig_filename=os.path.join(fig_dir,"%s.pdf"%os.path.basename(grid_filename))

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
  max_val=f_stack['max_val']
  max_x=f_stack['max_x']
  max_y=f_stack['max_y']
  max_z=f_stack['max_z']
  
  plt.clf()

  # set up the 4 axes
  x=np.arange(nx)*dx
  y=np.arange(ny)*dy
  z=(np.arange(nz)*dz+z_orig)*(-1)
  # setup of t-axis depends on type of stack_start_time
  if type(stack_start_time)==float:
    t=np.arange(nt)*dt+stack_start_time
  else:
    t=np.arange(nt)*dt

  # do plot
  plt.subplot(3,3,1)
  #p=plt.contourf(x,y,xy_cut.T)
  p=plt.imshow(xy_cut.T,origin='lower',interpolation='none',extent=[np.min(x),np.max(x),np.min(y),np.max(y)])
  #p=plt.imshow(x,y,xy_cut.T,origin='lower',extent=[np.min(x),np.max(x),np.min(y),np.max(y)])
  #plt.plot(x[ix_true],y[iy_true],'wo',markersize=20, alpha=0.5)
  plt.xlabel('x (km wrt origin)')
  plt.ylabel('y (km wrt to origin)')
  plt.title('XY plane')
  plt.subplot(3,3,2)
  #p=plt.contourf(x,z,xz_cut.T)
  p=plt.imshow(xz_cut.T,origin='upper',interpolation='none',extent=[np.min(x),np.max(x),np.min(z),np.max(z)])
  #plt.plot(x[ix_true],z[iz_true],'wo',markersize=20, alpha=0.5)
  plt.xlabel('x (km wrt origin)')
  plt.ylabel('z (km up)')
  plt.title('XZ plane')
  plt.subplot(3,3,3)
  #p=plt.contourf(y,z,yz_cut.T)
  p=plt.imshow(yz_cut.T,origin='upper',interpolation='none',extent=[np.min(y),np.max(y),np.min(z),np.max(z)])
  #plt.plot(y[iy_true],z[iz_true],'wo',markersize=20, alpha=0.5)
  plt.xlabel('y (km wrt origin)')
  plt.ylabel('z (km up)')
  plt.title('YZ plane')
  #plt.colorbar(p)


  llim = t[it_true]-2.0
  rlim = t[it_true]+2.0
  p=plt.subplot(3,1,2)
  plt.plot(t,max_val)
  #plt.xticks([llim,t[it_true],rlim])
  plt.xlabel('t (s)')
  plt.ylabel('Stack max ')
  plt.title('Maximum of stack')
  p.set_xlim(llim,rlim)
  p.set_ylim(0,max(max_val))
#  plt.hlines(loclevel,llim,rlim,'r',linewidth=2)
  plt.vlines(t[it_true],0,max(max_val),'r',linewidth=2)

  # put the origin back in for the last plots
  x=np.arange(nx)*dx+x_orig
  y=np.arange(ny)*dy+y_orig
  z=np.arange(nz)*dz+z_orig
  # plot max x
  p=plt.subplot(3,3,7)
  plt.plot(t,max_x)
  plt.xticks([llim,t[it_true],rlim])
  plt.xlabel('t (s)')
  plt.ylabel('x (km) ')
  plt.title('x at maximum')
  p.set_xlim(llim,rlim)
  plt.hlines(x[ix_true],llim,rlim,'r',linewidth=2)
  plt.vlines(t[it_true],min(max_x),max(max_x),'r',linewidth=2)
  # plot max y
  p=plt.subplot(3,3,8)
  plt.plot(t,max_y)
  plt.xticks([llim,t[it_true],rlim])
  plt.xlabel('t (s)')
  plt.ylabel('y (km) ')
  plt.title('y at maximum')
  p.set_xlim(llim,rlim)
  plt.hlines(y[iy_true],llim,rlim,'r',linewidth=2)
  plt.vlines(t[it_true],min(max_y),max(max_y),'r',linewidth=2)
  # plot max z
  p=plt.subplot(3,3,9)
  plt.plot(t,max_z)
  plt.xticks([llim,t[it_true],rlim])
  plt.xlabel('t (s)')
  plt.ylabel('z (km down) ')
  plt.title('z at maximum')
  p.set_xlim(llim,rlim)
  plt.hlines(z[iz_true],llim,rlim,'r',linewidth=2)
  plt.vlines(t[it_true],min(max_z),max(max_z),'r',linewidth=2)

  plt.tight_layout()
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
 
