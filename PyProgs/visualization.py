import numpy as np
from integrate4D import *


def setup_test_grid():

  dims=(20,30,40,50)
  grid4D=np.zeros(dims)
  x0=np.linspace(0,4,dims[0])
  x1=np.linspace(0,6,dims[1])
  x2=np.linspace(0,8,dims[2])
  x3=np.linspace(0,10,dims[3])
  grid4D[1,2,3,7]=2
  grid4D[1,2,4,7]=4  # add something interesting to find in the matrix
  grid4D[1,2,5,7]=2

  grid4D = grid4D / compute_integral4D(grid4D,x0,x1,x2,x3)

  exp0,exp1,exp2,exp3,cov_matrix,grid_dict = compute_expected_coordinates4D(grid4D,x0,x1,x2,x3,return_2Dgrids=True)

  
  return grid4D,grid_dict

