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
  suite.addTest(LocationTests('test_locations_trigger_fullRes'))
  suite.addTest(LocationTests('test_locations_prob'))
  suite.addTest(LocationTests('test_locations_prob_fullRes'))
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


#  @unittest.skip('Not bothering with trigger test')
  def test_locations_trigger(self):

    self.wo.opdict['outdir']='TEST'
    self.wo.verify_location_options()

    base_path=self.wo.opdict['base_path']
    test_datadir=self.wo.opdict['test_datadir']
    outdir=self.wo.opdict['outdir']

    exp_loc_fname = os.path.join(base_path,test_datadir,'TEST_locations.dat')
    exp_loc_file = open(exp_loc_fname,'r') 
    exp_lines=exp_loc_file.readlines()

    do_locations_trigger_setup_and_run(self.wo.opdict)

    loc_fname = os.path.join(base_path,'out',outdir,'loc','locations.dat')
    loc_file = open(loc_fname,'r') 
    lines=loc_file.readlines()

    self.assertEquals(lines,exp_lines)

  def test_locations_trigger_fullRes(self):

    self.wo.opdict['outdir']='TEST_fullRes'
    self.wo.verify_location_options()

    base_path=self.wo.opdict['base_path']
    test_datadir=self.wo.opdict['test_datadir']
    outdir=self.wo.opdict['outdir']

    exp_loc_fname = os.path.join(base_path,test_datadir,'TEST_fullRes_locations.dat')
    exp_loc_file = open(exp_loc_fname,'r') 
    exp_lines=exp_loc_file.readlines()

    do_locations_trigger_setup_and_run(self.wo.opdict)

    loc_fname = os.path.join(base_path,'out',outdir,'loc','locations.dat')
    loc_file = open(loc_fname,'r') 
    lines=loc_file.readlines()

    self.assertEquals(lines,exp_lines)

#  @unittest.skip('Not bothering with low res test')
  def test_locations_prob(self):

    self.wo.opdict['search_grid'] = 'test_grid.search.hdr'
    self.wo.opdict['outdir']='TEST'

    self.wo.verify_location_options()
    self.wo.verify_migration_options()

    base_path=self.wo.opdict['base_path']
    test_datadir=self.wo.opdict['test_datadir']
    outdir=self.wo.opdict['outdir']

    exp_loc_fname = os.path.join(base_path,test_datadir,'TEST_locations_prob.dat')
    exp_loc_file = open(exp_loc_fname,'r') 
    exp_lines=exp_loc_file.readlines()

    do_locations_prob_setup_and_run(self.wo.opdict)

    loc_fname = os.path.join(base_path,'out',outdir,'loc','locations_prob.dat')
    loc_file = open(loc_fname,'r') 
    lines=loc_file.readlines()

    self.assertEquals(lines,exp_lines)


#  @unittest.skip('Not bothering with high res test')
  def test_locations_prob_fullRes(self):

    self.wo.opdict['search_grid'] = 'grid.Taisne.search.hdr'
    self.wo.opdict['outdir']='TEST_fullRes'

    self.wo.verify_location_options()
    self.wo.verify_migration_options()

    base_path=self.wo.opdict['base_path']
    test_datadir=self.wo.opdict['test_datadir']
    outdir=self.wo.opdict['outdir']

    exp_loc_fname = os.path.join(base_path,test_datadir,'TEST_fullRes_locations_prob.dat')
    exp_loc_file = open(exp_loc_fname,'r') 
    exp_lines=exp_loc_file.readlines()

    do_locations_prob_setup_and_run(self.wo.opdict)

    loc_fname = os.path.join(base_path,'out',outdir,'loc','locations_prob.dat')
    loc_file = open(loc_fname,'r') 
    lines=loc_file.readlines()

    self.assertEquals(lines,exp_lines)

if __name__ == '__main__':

  import logging
  logging.basicConfig(level=logging.INFO, format='%(levelname)s : %(asctime)s : %(message)s')
 
  unittest.TextTestRunner(verbosity=2).run(suite())
 
