import h5py, os, unittest
import numpy as np
from hdf5_grids import *

def suite():
  suite = unittest.TestSuite()
  suite.addTest(H5Tests('test_RandomRead'))
  suite.addTest(H5Tests('test_compression'))
  suite.addTest(H5SingleGridTests('test_init_del'))
  suite.addTest(H5SingleGridTests('test_NllReadHdr'))
  suite.addTest(H5SingleGridTests('test_NllRead'))
  return suite

class H5Tests(unittest.TestCase):

  def test_RandomRead(self):

    # set up some random data
    data=np.random.rand(100,200,50)

    # put this data in a hdf5 file
    filename='randomtest.hdf5'
    f=h5py.File(filename,'w')
    dset = f.create_dataset('random',data=data)
    f.close()

    # re-read the information from the file
    f=h5py.File(filename,'r')
    dset=f['random']

    # assert statements
    self.assertEqual(data.shape,dset.shape)
    np.testing.assert_allclose(data,dset)

    # close and remove the hdf5 file
    f.close()
    os.remove(filename)

  def test_compression(self):

    # set up some random data
    #data=np.random.rand(100,200,50)
    data=np.ones((100,200,50))

    # put this data in a hdf5 file
    filename='randomtest.hdf5'
    f=h5py.File(filename,'w')
    dset = f.create_dataset('random',data=data)
    f.close()
    size_none=os.stat(filename).st_size
    os.remove(filename)

    # put this data in a hdf5 file
    f=h5py.File(filename,'w')
    dset = f.create_dataset('random',data=data,compression='lzf')
    f.close()
    size_lzf=os.stat(filename).st_size
    os.remove(filename)

    # assert statements
    self.assertLess(size_lzf,size_none)

class H5SingleGridTests(unittest.TestCase):

  def test_init_del(self):

    data=np.ones((100,200,50))
    info={}
    info['nx']=100
    info['ny']=200
    info['nz']=50

    # put this data in a hdf5 file
    filename='init_del.hdf5'
    sg=H5SingleGrid(filename,data,info)

    # assert statements
    self.assertTrue(os.path.isfile(filename))
    self.assertEqual(data.shape,sg.grid_data.shape)
    np.testing.assert_allclose(data,sg.grid_data)
    self.assertEqual(info['ny'],sg.grid_info['ny'])

    # now get rid of the grid
    del sg
    #self.assertRaisesRegexp(UnboundLocalError,'sg',type,sg)

    sg1=H5SingleGrid(filename)
    self.assertEqual(data.shape,sg1.grid_data.shape)
    np.testing.assert_allclose(data,sg1.grid_data)
    self.assertEqual(info['ny'],sg1.grid_info['ny'])

    del sg1
    # clean up file
    os.remove(filename)

  def test_NllReadHdr(self):

    from NllGridLib import read_hdr_file

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

  def test_NllRead(self):
    from array import array

    base_path=os.getenv('WAVELOC_PATH')
    nll_name=os.path.join(base_path,'test_data','test.time')
    nll_hdr_file="%s.hdr"%nll_name 
    nll_buf_file="%s.buf"%nll_name 

    info=read_hdr_file(nll_hdr_file)
    nx=info['nx']
    ny=info['ny']
    nz=info['nz']

    f=open(nll_buf_file,'rb')
    buf=array('f')
    buf.fromfile(f,nx*ny*nz)
    f.close()
    np_buf=np.array(buf)

    b_index=np.random.randint(0,nx)*np.random.randint(0,ny)*np.random.randint(0,nz)

    sg=H5NllSingleGrid('test_nll.hdf5',nll_name)
    self.assertEqual(info['ny'],sg.grid_info['ny'])
    self.assertAlmostEqual(info['orig_lat'],sg.grid_info['orig_lat'])
    self.assertEqual(np_buf.shape,sg.grid_data.shape)
    self.assertAlmostEqual(np_buf[b_index],sg.grid_data[b_index])

    del sg

if __name__ == '__main__':

  import logging
  logging.basicConfig(level=logging.INFO, format='%(levelname)s : %(asctime)s : %(message)s')
 
  unittest.TextTestRunner(verbosity=2).run(suite())
 
