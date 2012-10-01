import h5py, os, unittest
import numpy as np
from NllGridLib import *

def suite():
  suite = unittest.TestSuite()
  suite.addTest(NLLTests('test_NllReadHdr'))
  suite.addTest(NLLTests('test_ReadStations'))
  suite.addTest(ProjectionTests('test_TransUnknown'))
  suite.addTest(ProjectionTests('test_TransGlobal'))
  suite.addTest(ProjectionTests('test_TransSimple'))
  return suite

class ProjectionTests(unittest.TestCase):

  def setUp(self):
    self.true_lat=-90+np.random.ranf()*180
    self.true_lon=-180+np.random.ranf()*360

    self.proj_info={}
    self.proj_info['orig_lat'] = self.true_lat + np.random.ranf()*10.0 - 5.0
    self.proj_info['orig_lon'] = self.true_lon + np.random.ranf()*10.0 - 5.0
    self.proj_info['map_rot'] =  np.random.ranf()*36.0 - 180.0

  def test_TransGlobal(self):
    x,y = latlon2rect('TRANS_GLOBAL',self.true_lat,self.true_lon)
    lat,lon = rect2latlon('TRANS_GLOBAL',x,y)
    self.assertAlmostEqual(self.true_lat,lat,6)
    self.assertAlmostEqual(self.true_lon,lon,6)

  def test_TransUnknown(self):
    self.assertRaises(UserWarning,latlon2rect,'TRANS_BUG',self.true_lat,self.true_lon)
    self.assertRaises(UserWarning,rect2latlon,'TRANS_BUG',self.true_lat,self.true_lon)

  def test_TransSimple(self):
    x,y = latlon2rect('TRANS_SIMPLE',self.true_lat,self.true_lon,self.proj_info)
    lat,lon = rect2latlon('TRANS_SIMPLE',x,y,self.proj_info)
    self.assertAlmostEqual(self.true_lat,lat,6)
    self.assertAlmostEqual(self.true_lon,lon,6)


class NLLTests(unittest.TestCase):

  def test_NllReadHdr(self):

    base_path=os.getenv('WAVELOC_PATH')
    nll_hdr_file=os.path.join(base_path,'test_data','test_grid.search.geo.hdr')
    info=read_hdr_file(nll_hdr_file)
    self.assertEqual(info['nx'],51)
    self.assertEqual(info['y_orig'],-10.)
    self.assertEqual(info['dz'],2.)
    self.assertEqual(info['proj_name'],'TRANS_SIMPLE')
    self.assertAlmostEqual(info['orig_lat'],44.727000)
    self.assertAlmostEqual(info['orig_lon'],11.086000)
    self.assertAlmostEqual(info['map_rot'],0.)

  def test_ReadStations(self):
    base_path=os.getenv('WAVELOC_PATH')
    stations_file=os.path.join(base_path,'test_data','coord_stations_test')
    stations=read_stations_file(stations_file)
    self.assertTrue(stations.has_key('RVL'))
    self.assertAlmostEqual(stations['SNE']['x'],366.958)


if __name__ == '__main__':

  import logging
  logging.basicConfig(level=logging.INFO, format='%(levelname)s : %(asctime)s : %(message)s')
 
  unittest.TextTestRunner(verbosity=2).run(suite())
 
