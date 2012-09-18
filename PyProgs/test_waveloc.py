import unittest
import os, glob
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
 
  # make link for test grid file etc
  test_files=['test_grid.search.hdr', 'coord_stations_test']
  for tfile in test_files:
    try:
      os.symlink(os.path.join(base_path,'test_data',tfile),os.path.join(base_path,'lib',tfile))
      logging.info("Linked %s"%tfile)
    except OSError:
      logging.info("File %s already linked"%'tfile')
      logging.info("Removing old %s"%tfile)
      os.remove(os.path.join(base_path,'lib',tfile))
      os.symlink(os.path.join(base_path,'test_data',tfile),os.path.join(base_path,'lib',tfile))
      logging.info("Linked %s"%tfile)

  # make links for PDF time grids
  test_files=glob.glob(os.path.join(base_path,'test_data', 'time_grids', 'Slow*'))
  if test_files==[]: 
    logging.error('Dowload https://github.com/downloads/amaggi/waveloc/TEST_time_grids.tgz and unpack it in the %s directory, then re-run'%os.path.join(base_path,'test_data'))
  for tfile in test_files:
    try:
      os.symlink(os.path.join(base_path,'test_data','time_grids',os.path.basename(tfile)),os.path.join(base_path,'lib',os.path.basename(tfile)))
      logging.info("Linked %s"%tfile)
     except OSError:
       pass

   


class SetupTests(unittest.TestCase):

  def test_setup(self):
    self.assertTrue(True)

if __name__ == '__main__':

  import test_processing, test_migration, test_location
  import logging
  logging.basicConfig(level=logging.INFO, format='%(levelname)s : %(asctime)s : %(message)s')
 

  suite01 = suite()
  suite02 = test_processing.suite()
  suite03 = test_migration.suite()
  suite04 = test_location.suite()

  alltests=unittest.TestSuite([suite01, suite02, suite03, suite04])

  unittest.TextTestRunner(verbosity=2).run(alltests)
 
