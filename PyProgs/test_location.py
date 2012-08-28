import unittest
import os, glob
from locations_trigger import do_locations_trigger_setup_and_run
from locations_prob import do_locations_prob_setup_and_run
from OP_waveforms import Waveform
from integrate4D import *

def suite():
  suite = unittest.TestSuite()
  suite.addTest(IntegrationTests('test_integration'))
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
    dims=(40,60,80,100)
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

    grid_norm = compute_integral4D(grid4D,x0,x1,x2,x3)
    grid4D=grid4D/grid_norm

    exp0,epx1,exp2,exp3 = compute_expected_coordintes4D(grid4D,x0,x1,x2,x3)

    self.assertAlmostEqual(my_exp0,exp0,7)
    self.assertAlmostEqual(my_exp1,exp1,7)
    self.assertAlmostEqual(my_exp2,exp2,7)
    self.assertAlmostEqual(my_exp3,exp3,7)


class LocationTests(unittest.TestCase):

  def setUp(self):

    self.base_path=os.getenv('WAVELOC_PATH')

    self.test_datadir='test_data'
    self.datadir='TEST'
    self.net_list='YA'
    self.sta_list="FJS,FLR,FOR,HDL,RVL,SNE,UV01,UV02,UV03,UV04,UV05,UV06,UV07,UV08,UV09,UV10,UV11,UV12,UV13,UV14,UV15"
    self.comp_list="HHZ"
    self.starttime="2010-10-14T00:14:00.0Z"
    self.endtime="2010-10-14T00:18:00.0Z"
    self.time_grid='Slow_len.100m.P'
    self.search_grid='test_grid.search.hdr'
    self.coord_stations='coord_stations_test'
    self.kwin=4
    self.loclevel=50
    self.snr_limit=10.0
    self.sn_time=10.0
    self.dataglob="*kurt.mseed"



  def test_locations_trigger(self):

    do_locations_trigger_setup_and_run(base_path=self.base_path, outdir=self.datadir, datadir=self.datadir, dataglob=self.dataglob, reloc=False, loclevel=self.loclevel,snr_limit=self.snr_limit,sn_time=self.sn_time,n_kurt_min=self.kwin)
    self.assertTrue(True)

  def test_locations_prob(self):

    do_locations_prob_setup_and_run(base_path=self.base_path, outdir=self.datadir, loclevel=self.loclevel, datadir=self.datadir, dataglob=self.dataglob, snr_limit=self.snr_limit, sn_time=self.sn_time)
    self.assertTrue(True)

if __name__ == '__main__':

  import logging
  logging.basicConfig(level=logging.DEBUG, format='%(levelname)s : %(asctime)s : %(message)s')
 
  unittest.TextTestRunner(verbosity=2).run(suite())
 
