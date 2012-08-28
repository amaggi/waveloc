import unittest
import os, glob
from locations_trigger import do_locations_trigger_setup_and_run
from locations_prob import do_locations_prob_setup_and_run
from OP_waveforms import Waveform

def suite():
  suite = unittest.TestSuite()
  suite.addTest(LocationTests('test_locations_trigger'))
  suite.addTest(LocationTests('test_locations_prob'))
  return suite

    
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
 
