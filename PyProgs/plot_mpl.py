import os, h5py, logging
import numpy as np
import scipy.integrate as si
import matplotlib.pyplot as plt
from integrate4D import *
from filters import smooth

def plotLocationGrid(loc,grid_info,stack_grid,fig_dir):

  # set up plot using info from grid_info
  nx,ny,nz,nt = grid_info['grid_shape']
  dx,dy,dz,dt = grid_info['grid_spacing']
  x_orig,y_orig,z_orig = grid_info['grid_orig']
  stack_starttime=grid_info['start_time']
  stack_otime=grid_info['stack_otime']
  grid_filename=grid_info['dat_file']

  fig_filename=os.path.join(fig_dir,"%s.pdf"%os.path.basename(grid_filename))

  # read the stack file
#  stack_grid=np.fromfile(grid_filename).reshape(nx,ny,nz,nt)

  # set up the 4 axes
  x=np.arange(nx)*dx+x_orig
  y=np.arange(ny)*dy+y_orig
  z=np.arange(nz)*dz+z_orig
  t=np.arange(nt)*dt

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

  #get index of origin time 
  it_true=np.int(np.round((o_time-stack_starttime)/dt))
  ix_true=np.int(np.round((x_mean-x_orig)/dx))
  iy_true=np.int(np.round((y_mean-y_orig)/dy))
  iz_true=np.int(np.round((z_mean-z_orig)/dz))

  it_left=np.int(np.round((o_time-o_err_left-stack_starttime)/dt))
  it_right=np.int(np.round((o_time+o_err_right-stack_starttime)/dt))
  ix_low=np.int(np.round((x_mean-x_sigma-x_orig)/dx))
  iy_low=np.int(np.round((y_mean-y_sigma-y_orig)/dy))
  iz_low=np.int(np.round((z_mean-z_sigma-z_orig)/dz))
  ix_high=np.int(np.round((x_mean+x_sigma-x_orig)/dx))
  iy_high=np.int(np.round((y_mean+y_sigma-y_orig)/dy))
  iz_high=np.int(np.round((z_mean+z_sigma-z_orig)/dz))

  # cut through the true location at the true time 
  xy_cut=stack_grid[:,:,iz_true,it_true]
  xz_cut=stack_grid[:,iy_true,:,it_true]
  yz_cut=stack_grid[ix_true,:,:,it_true]

  # extract the max stacks
  max_val=smooth(stack_grid.max(0).max(0).max(0))
  max_x=stack_grid.max(2).max(1).argmax(0)*dx + x_orig
  max_y=stack_grid.max(2).max(0).argmax(0)*dy + y_orig
  max_z=stack_grid.max(1).max(0).argmax(0)*dz + z_orig

  print it_true, np.argmax(max_val)

  # do plot
  plt.subplot(3,3,1)
  p=plt.contourf(x,y,xy_cut.T)
  #plt.plot(x[ix_true],y[iy_true],'wo',markersize=20, alpha=0.5)
  plt.xlabel('x (km)')
  plt.ylabel('y (km)')
  plt.title('XY plane')
  plt.subplot(3,3,2)
  p=plt.contourf(x,z,xz_cut.T)
  #plt.plot(x[ix_true],z[iz_true],'wo',markersize=20, alpha=0.5)
  plt.xlabel('x (km)')
  plt.ylabel('z (km)')
  plt.title('XZ plane')
  plt.subplot(3,3,3)
  p=plt.contourf(y,z,yz_cut.T)
  #plt.plot(y[iy_true],z[iz_true],'wo',markersize=20, alpha=0.5)
  plt.xlabel('y (km)')
  plt.ylabel('z (km)')
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
  plt.vlines(t[it_left],0,max(max_val),'r',linewidth=1)
  plt.vlines(t[it_right],0,max(max_val),'r',linewidth=1)

  # plot max x
  p=plt.subplot(3,3,7)
  plt.plot(t,max_x)
  plt.xticks([llim,t[it_true],rlim])
  plt.xlabel('t (s)')
  plt.ylabel('x (km) ')
  plt.title('x at maximum')
  p.set_xlim(llim,rlim)
  plt.hlines(x[ix_true],llim,rlim,'r',linewidth=2)
  plt.hlines(x[ix_low],llim,rlim,'r',linewidth=1)
  plt.hlines(x[ix_high],llim,rlim,'r',linewidth=1)
  plt.vlines(t[it_true],min(max_x),max(max_x),'r',linewidth=2)
  plt.vlines(t[it_left],min(max_x),max(max_x),'r',linewidth=1)
  plt.vlines(t[it_right],min(max_x),max(max_x),'r',linewidth=1)
  # plot max y
  p=plt.subplot(3,3,8)
  plt.plot(t,max_y)
  plt.xticks([llim,t[it_true],rlim])
  plt.xlabel('t (s)')
  plt.ylabel('y (km) ')
  plt.title('y at maximum')
  p.set_xlim(llim,rlim)
  plt.hlines(y[iy_true],llim,rlim,'r',linewidth=2)
  plt.hlines(y[iy_low],llim,rlim,'r',linewidth=1)
  plt.hlines(y[iy_high],llim,rlim,'r',linewidth=1)
  plt.vlines(t[it_true],min(max_y),max(max_y),'r',linewidth=2)
  plt.vlines(t[it_left],min(max_y),max(max_y),'r',linewidth=1)
  plt.vlines(t[it_right],min(max_y),max(max_y),'r',linewidth=1)
  # plot max z
  p=plt.subplot(3,3,9)
  plt.plot(t,max_z)
  plt.xticks([llim,t[it_true],rlim])
  plt.xlabel('t (s)')
  plt.ylabel('z (km) ')
  plt.title('z at maximum')
  p.set_xlim(llim,rlim)
  plt.hlines(z[iz_true],llim,rlim,'r',linewidth=2)
  plt.hlines(z[iz_low],llim,rlim,'r',linewidth=1)
  plt.hlines(z[iz_high],llim,rlim,'r',linewidth=1)
  plt.vlines(t[it_true],min(max_z),max(max_z),'r',linewidth=2)
  plt.vlines(t[it_left],min(max_z),max(max_z),'r',linewidth=1)
  plt.vlines(t[it_right],min(max_z),max(max_z),'r',linewidth=1)

  plt.tight_layout()
  plt.savefig(fig_filename)



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
  z=np.arange(nz)*dz
  t=np.arange(nt)*dt+stack_start_time

  # do plot
  plt.subplot(3,3,1)
  p=plt.contourf(x,y,xy_cut.T)
  #plt.plot(x[ix_true],y[iy_true],'wo',markersize=20, alpha=0.5)
  plt.xlabel('x (km)')
  plt.ylabel('y (km)')
  plt.title('XY plane')
  plt.subplot(3,3,2)
  p=plt.contourf(x,z,xz_cut.T)
  #plt.plot(x[ix_true],z[iz_true],'wo',markersize=20, alpha=0.5)
  plt.xlabel('x (km)')
  plt.ylabel('z (km)')
  plt.title('XZ plane')
  plt.subplot(3,3,3)
  p=plt.contourf(y,z,yz_cut.T)
  #plt.plot(y[iy_true],z[iz_true],'wo',markersize=20, alpha=0.5)
  plt.xlabel('y (km)')
  plt.ylabel('z (km)')
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
  plt.ylabel('z (km) ')
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
 
