import numpy as np
import scipy.integrate as si

def compute_integral4D(grid4D,x0,x1,x2,x3):
  grid_norm=si.trapz(si.trapz(si.trapz(si.trapz(grid4D,x=x0,axis=0),x=x1,axis=0),x=x2,axis=0),x=x3,axis=0)
  return grid_norm

def compute_expected_coordinates4D(grid4D,x0,x1,x2,x3):
  # expects a normalized grid  

  prob_x0=si.trapz(si.trapz(si.trapz(grid4D_norm,x=x1,axis=1),x=x2,axis=1),x=x3,axis=1)
  exp_x0=si.trapz(x0*prob_x0,x=x0)

  prob_x1=si.trapz(si.trapz(si.trapz(grid4D_norm,x=x0,axis=0),x=x2,axis=1),x=x3,axis=1)
  exp_x1=si.trapz(x1*prob_x1,x=x1)
   
  prob_x2=si.trapz(si.trapz(si.trapz(grid4D_norm,x=x0,axis=0),x=x1,axis=0),x=x3,axis=1)
  exp_x2=si.trapz(x2*prob_x2,x=x2)
 
  prob_x3=si.trapz(si.trapz(si.trapz(grid4D_norm,x=x0,axis=0),x=x1,axis=0),x=x2,axis=0)
  exp_x3=si.trapz(x3*prob_x3,x=x3)

  return x0,x1,x2,x3
