import numpy as np
from integrate4D import *


def setup_test_grid():

  dims=(16,12,6,200)
  grid4D=np.zeros(dims)
  x0=np.linspace(0,8,dims[0])
  x1=np.linspace(0,6,dims[1])
  x2=np.linspace(0,3,dims[2])
  x3=np.linspace(0,20,dims[3])
  grid4D[8,4,3,60]=60  # add something interesting to find in the matrix

  #grid4D = grid4D / compute_integral4D(grid4D,x0,x1,x2,x3)

  exp0,exp1,exp2,exp3,cov_matrix,grid_dict = compute_expected_coordinates4D(grid4D,x0,x1,x2,x3,return_2Dgrids=True)

  
  return grid4D,grid_dict,[x0,x1,x2,x3]

