import unittest
import os

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
 

class SetupTests(unittest.TestCase):

  def test_setup(self):
    self.assertTrue(True)

if __name__ == '__main__':

  import test_processing

  suite01 = suite()
  suite02 = test_processing.suite()

  alltests=unittest.TestSuite([suite01, suite02])

  unittest.TextTestRunner(verbosity=2).run(alltests)
 
