import unittest
import os
import logging

def suite():
  suite = unittest.TestSuite()
  suite.addTest(SetupTests('test_setup'))
  return suite

def setUpModule():

  from make_SDS_data_links import make_SDS_data_links

  base_path=os.getenv('WAVELOC_PATH')
  test_data_dir=os.path.join(base_path,'test_data','raw_data')
  data_dir=os.path.join(base_path,'data','TEST')
  make_SDS_data_links(test_data_dir,'*MSEED',data_dir)
 
  # make link for test grid file
  try:
    os.symlink(os.path.join(base_path,'test_data','test_grid.search.hdr'),os.path.join(base_path,'aux','test_grid.search.hdr'))
  except OSError:
    logging.debug("File %s already linked"%'test_grid_search.hdr')
  try:
    os.symlink(os.path.join(base_path,'test_data','coord_stations_test'),os.path.join(base_path,'aux','coord_stations_test'))
  except OSError:
    logging.debug("File %s already linked"%'coord_stations_test')
  

class SetupTests(unittest.TestCase):

  def test_setup(self):
    self.assertTrue(True)

if __name__ == '__main__':

  import test_processing, test_migration
  import logging
  logging.basicConfig(level=logging.DEBUG, format='%(levelname)s : %(asctime)s : %(message)s')
 

  suite01 = suite()
  suite02 = test_processing.suite()
  suite03 = test_migration.suite()

  alltests=unittest.TestSuite([suite01, suite02, suite03])

  unittest.TextTestRunner(verbosity=2).run(alltests)
 
