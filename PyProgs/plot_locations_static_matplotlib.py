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
  
  print files_dict

  plt.plot(x0,grid_dict['prob_x0'])
  plt.xlabel('x (km)')
  plt.ylabel('p(x)')
  plt.title('Marginal probability density over x')
  plt.savefig(files_dict['prob_x0'])
  


if __name__ == '__main__':

  base_name=os.getenv('WAVELOC_PATH')
  fig_dir = os.path.join(base_name,'test_figures')
  if not os.path.exists(fig_dir): os.makedirs(fig_dir)

  base_filename = os.path.join(fig_dir,'testplot_st_mpl')
  grid4D, grid_dict,x_list = setup_test_grid()
  plot_locations_static_matplotlib(grid4D,grid_dict,x_list,base_filename)
  

