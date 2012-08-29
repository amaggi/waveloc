import numpy as np
import scipy.integrate as si

def compute_integral4D(grid4D,x0,x1,x2,x3):
  grid_norm=si.trapz(si.trapz(si.trapz(si.trapz(grid4D,x=x0,axis=0),x=x1,axis=0),x=x2,axis=0),x=x3,axis=0)
  return grid_norm

def compute_expected_coordinates4D(grid4D,x0,x1,x2,x3,return_2Dgrids=False):
  # normalize grid
  grid4D = grid4D / compute_integral4D(grid4D,x0,x1,x2,x3)

  # get 1D marginals, expected values and variances
  prob_x0=si.trapz(si.trapz(si.trapz(grid4D,x=x1,axis=1),x=x2,axis=1),x=x3,axis=1)
  exp_x0=si.trapz(x0*prob_x0,x=x0)
  var_x0=si.trapz((x0-exp_x0)*(x0-exp_x0)*prob_x0,x=x0)

  prob_x1=si.trapz(si.trapz(si.trapz(grid4D,x=x0,axis=0),x=x2,axis=1),x=x3,axis=1)
  exp_x1=si.trapz(x1*prob_x1,x=x1)
  var_x1=si.trapz((x1-exp_x1)*(x1-exp_x1)*prob_x1,x=x1)
   
  prob_x2=si.trapz(si.trapz(si.trapz(grid4D,x=x0,axis=0),x=x1,axis=0),x=x3,axis=1)
  exp_x2=si.trapz(x2*prob_x2,x=x2)
  var_x2=si.trapz((x2-exp_x2)*(x2-exp_x2)*prob_x2,x=x2)
 
  prob_x3=si.trapz(si.trapz(si.trapz(grid4D,x=x0,axis=0),x=x1,axis=0),x=x2,axis=0)
  exp_x3=si.trapz(x3*prob_x3,x=x3)
  var_x3=si.trapz((x3-exp_x3)*(x3-exp_x3)*prob_x3,x=x3)

  # get 2D marginals and covariances
  prob_x0_x1=si.trapz(si.trapz(grid4D,x=x2,axis=2),x=x3,axis=2)
  cov_x0_x1=si.trapz((x0-exp_x0)*si.trapz((x1-exp_x1)*prob_x0_x1,x=x1,axis=1),x=x0,axis=0)

  prob_x0_x2=si.trapz(si.trapz(grid4D,x=x1,axis=1),x=x3,axis=2)
  cov_x0_x2=si.trapz((x0-exp_x0)*si.trapz((x2-exp_x2)*prob_x0_x2,x=x2,axis=1),x=x0,axis=0)

  prob_x0_x3=si.trapz(si.trapz(grid4D,x=x1,axis=1),x=x2,axis=1)
  cov_x0_x3=si.trapz((x0-exp_x0)*si.trapz((x3-exp_x3)*prob_x0_x3,x=x3,axis=1),x=x0,axis=0)

  prob_x1_x2=si.trapz(si.trapz(grid4D,x=x0,axis=0),x=x3,axis=2)
  cov_x1_x2=si.trapz((x1-exp_x1)*si.trapz((x2-exp_x2)*prob_x1_x2,x=x2,axis=1),x=x1,axis=0)

  prob_x1_x3=si.trapz(si.trapz(grid4D,x=x0,axis=0),x=x2,axis=1)
  cov_x1_x3=si.trapz((x1-exp_x1)*si.trapz((x3-exp_x3)*prob_x1_x3,x=x3,axis=1),x=x1,axis=0)

  prob_x2_x3=si.trapz(si.trapz(grid4D,x=x0,axis=0),x=x1,axis=0)
  cov_x2_x3=si.trapz((x2-exp_x2)*si.trapz((x3-exp_x3)*prob_x2_x3,x=x3,axis=1),x=x2,axis=0)
  

  cov_matrix=np.eye(4)
  cov_matrix[0,0]=var_x0
  cov_matrix[1,1]=var_x1
  cov_matrix[2,2]=var_x2
  cov_matrix[3,3]=var_x3

  cov_matrix[0,1]=cov_x0_x1
  cov_matrix[1,0]=cov_x0_x1

  cov_matrix[0,2]=cov_x0_x2
  cov_matrix[2,0]=cov_x0_x2

  cov_matrix[0,3]=cov_x0_x3
  cov_matrix[3,0]=cov_x0_x3

  cov_matrix[1,2]=cov_x1_x2
  cov_matrix[2,1]=cov_x1_x2

  cov_matrix[1,3]=cov_x1_x3
  cov_matrix[3,1]=cov_x1_x3

  cov_matrix[2,3]=cov_x2_x3
  cov_matrix[3,2]=cov_x2_x3

  if return_2Dgrids:
    grid_dict={}
    grid_dict['prob_x0']=prob_x0
    grid_dict['prob_x1']=prob_x1
    grid_dict['prob_x2']=prob_x2
    grid_dict['prob_x3']=prob_x3
    grid_dict['prob_x0_x1']=prob_x0_x1
    grid_dict['prob_x0_x2']=prob_x0_x2
    grid_dict['prob_x0_x3']=prob_x0_x3
    grid_dict['prob_x1_x2']=prob_x1_x2
    grid_dict['prob_x1_x3']=prob_x1_x3
    grid_dict['prob_x2_x3']=prob_x2_x3

    return exp_x0,exp_x1,exp_x2,exp_x3,cov_matrix,grid_dict
  else:
    return exp_x0,exp_x1,exp_x2,exp_x3,cov_matrix
