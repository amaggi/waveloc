import h5py, os, unittest
import numpy as np

def suite():
  suite = unittest.TestSuite()
  suite.addTest(H5Tests('test_RandomRead'))
  suite.addTest(H5Tests('test_compression'))
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

if __name__ == '__main__':

  import logging
  logging.basicConfig(level=logging.INFO, format='%(levelname)s : %(asctime)s : %(message)s')
 
  unittest.TextTestRunner(verbosity=2).run(suite())
 
