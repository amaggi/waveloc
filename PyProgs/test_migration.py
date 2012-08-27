import unittest
import os, glob
from migration import do_migration_setup_and_run
from OP_waveforms import Waveform
from test_processing import waveforms_to_signature

def suite():
  suite = unittest.TestSuite()
  suite.addTest(MigrationTests('test_migration'))
  return suite

    
class MigrationTests(unittest.TestCase):

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
    self.data_length=600
    self.data_overlap=20


  @unittest.expectedFailure
  def test_migration(self):

    expected_signature_filename = os.path.join(self.base_path,self.test_datadir,'test_stack_signature.dat')
    expected_signature_file = open(expected_signature_filename,'r') 
    expected_lines=expected_signature_file.readlines()

    do_migration_setup_and_run(base_path=self.base_path,runtime=True,verbose=True, twoD=False, time_grid=self.time_grid, search_grid=self.search_grid, stations=self.coord_stations, outdir=self.datadir, datadir=self.datadir, dataglob='*kurt_grad.mseed', starttime=self.starttime, endtime=self.endtime, data_length=self.data_length, data_overlap=self.data_overlap,load_ttimes_buf=True)

    waveforms_to_signature(self.base_path,os.path.join('out',self.datadir,'stack'),'*mseed','stack_signature.dat')
    signature_filename=os.path.join(self.base_path,'out',self.datadir,'stack','stack_signature.dat')
    signature_file = open(signature_filename,'r') 
    lines=signature_file.readlines()

    self.assertSequenceEqual(lines,expected_lines)


if __name__ == '__main__':

  import logging
  logging.basicConfig(level=logging.DEBUG, format='%(levelname)s : %(asctime)s : %(message)s')
 
  unittest.TextTestRunner(verbosity=2).run(suite)
 
