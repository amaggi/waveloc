import os, glob, unittest
from options import WavelocOptions
from OP_waveforms import Waveform
from locations_trigger import do_locations_trigger_setup_and_run
from locations_prob import do_locations_prob_setup_and_run
from integrate4D import *

def suite():
  suite = unittest.TestSuite()
  suite.addTest(IntegrationTests('test_integration'))
  suite.addTest(IntegrationTests('test_expected_values'))
  suite.addTest(LocationTests('test_locations_trigger'))
  suite.addTest(LocationTests('test_locations_prob'))
  return suite

class IntegrationTests(unittest.TestCase):

  def test_integration(self):
    dims=(40,60,80,100)
    grid4D=np.ones(dims)
    x0=np.linspace(0,4,dims[0])
    x1=np.linspace(0,6,dims[1])
    x2=np.linspace(0,8,dims[2])
    x3=np.linspace(0,10,dims[3])
    grid_area=4.0*6*8*10

    grid_integral=compute_integral4D(grid4D,x0,x1,x2,x3)

    self.assertAlmostEqual(grid_area,grid_integral,7)
    
  def test_expected_values(self):
    dims=(20,30,40,50)
    grid4D=np.zeros(dims)
    x0=np.linspace(0,4,dims[0])
    x1=np.linspace(0,6,dims[1])
    x2=np.linspace(0,8,dims[2])
    x3=np.linspace(0,10,dims[3])
    grid4D[1,2,3,7]=2
    grid4D[1,2,4,7]=4  # add something interesting to find in the matrix
    grid4D[1,2,5,7]=2
    my_exp0=x0[1]
    my_exp1=x1[2]
    my_exp2=x2[4]
    my_exp3=x3[7]


    #grid4D = grid4D / compute_integral4D(grid4D,x0,x1,x2,x3)
    exp0,exp1,exp2,exp3,cov_matrix = compute_expected_coordinates4D(grid4D,x0,x1,x2,x3)
    var_x0=cov_matrix[0,0]
    var_x1=cov_matrix[1,1]
    var_x2=cov_matrix[2,2]
    var_x3=cov_matrix[3,3]

    self.assertAlmostEqual(my_exp0,exp0,7)
    self.assertAlmostEqual(my_exp1,exp1,7)
    self.assertAlmostEqual(my_exp2,exp2,7)
    self.assertAlmostEqual(my_exp3,exp3,7)

    self.assertAlmostEqual(var_x0,0.0,7)
    self.assertAlmostEqual(var_x1,0.0,7)
    self.assertAlmostEqual(var_x2,0.0210,4)
    self.assertAlmostEqual(var_x3,0.0,7)

class LocationTests(unittest.TestCase):

  def setUp(self):

    self.wo=WavelocOptions()
    self.wo.set_test_options()
    self.wo.verify_location_options()


  def test_locations_trigger(self):

    do_locations_trigger_setup_and_run(self.wo.opdict)
    self.assertTrue(True)

  def test_locations_prob(self):

    #self.wo.opdict['loclevel']=50
    do_locations_prob_setup_and_run(self.wo.opdict)
    self.assertTrue(True)

if __name__ == '__main__':

  import logging
  logging.basicConfig(level=logging.INFO, format='%(levelname)s : %(asctime)s : %(message)s')
 
  unittest.TextTestRunner(verbosity=2).run(suite())
 
