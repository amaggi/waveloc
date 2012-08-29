import os
import numpy as np
import matplotlib.pyplot as plt
from visualization import setup_test_grid

def plot_locations_static_matplotlib(grid4D,grid_dict,x_list,base_filename):

  print grid_dict.keys()

  x0,x1,x2,x3 = x_list

  # create filenames for plots
  files_dict={}
  for key in grid_dict.keys():
    files_dict[key] = base_filename + '_' + key + '.pdf'
  
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

  plt.clf()
  plt.plot(x3,grid_dict['prob_x3'])
  plt.xlabel('time (s)')
  plt.ylabel('p(t)')
  plt.title('Marginal probability density over time')
  plt.savefig(files_dict['prob_x3'])

  # 2D plots
  plt.clf()
  plt.contourf(x0,x3,grid_dict['prob_x0_x3'].T)
  plt.xlabel('x (km)')
  plt.ylabel('time (s)')
  plt.title('Marginal probability density over x and time')
  plt.savefig(files_dict['prob_x0_x3'])
  
  plt.clf()
  plt.contourf(x1,x3,grid_dict['prob_x1_x3'].T)
  plt.xlabel('y (km)')
  plt.ylabel('time (s)')
  plt.title('Marginal probability density over y and time')
  plt.savefig(files_dict['prob_x1_x3'])

  plt.clf()
  plt.contourf(x2,x3,grid_dict['prob_x2_x3'].T)
  plt.xlabel('z (km)')
  plt.ylabel('time (s)')
  plt.title('Marginal probability density over z and time')
  plt.savefig(files_dict['prob_x2_x3'])
 
  plt.clf()
  plt.contourf(x0,x1,grid_dict['prob_x0_x1'].T)
  plt.xlabel('x (km)')
  plt.ylabel('y (km)')
  plt.title('Marginal probability density over x and y')
  plt.savefig(files_dict['prob_x0_x1'])
 
  plt.clf()
  plt.contourf(x0,x2,grid_dict['prob_x0_x2'].T)
  plt.xlabel('x (km)')
  plt.ylabel('z (km)')
  plt.title('Marginal probability density over x and z')
  plt.savefig(files_dict['prob_x0_x2'])

  plt.clf()
  plt.contourf(x1,x2,grid_dict['prob_x1_x2'].T)
  plt.xlabel('y (km)')
  plt.ylabel('z (km)')
  plt.title('Marginal probability density over y and z')
  plt.savefig(files_dict['prob_x1_x2'])
 
if __name__ == '__main__':

  base_name=os.getenv('WAVELOC_PATH')
  fig_dir = os.path.join(base_name,'test_figures')
  if not os.path.exists(fig_dir): os.makedirs(fig_dir)

  base_filename = os.path.join(fig_dir,'testplot_st_mpl')
  grid4D, grid_dict,x_list = setup_test_grid()
  plot_locations_static_matplotlib(grid4D,grid_dict,x_list,base_filename)
  

