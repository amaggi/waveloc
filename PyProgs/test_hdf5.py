import h5py, os, unittest
import numpy as np
from hdf5_grids import *

def suite():
  suite = unittest.TestSuite()
  suite.addTest(H5Tests('test_RandomRead'))
  suite.addTest(H5Tests('test_Persistance'))
  suite.addTest(H5Tests('test_MemoryPerformance'))
  suite.addTest(H5Tests('test_compression'))
  suite.addTest(H5SingleGridTests('test_init_del'))
  suite.addTest(H5SingleGridTests('test_NllRead'))
  suite.addTest(H5SingleGridTests('test_nll2hdf5'))
  suite.addTest(H5SingleGridTests('test_interpolation_ones'))
  suite.addTest(H5SingleGridTests('test_interpolation_sinc'))
  suite.addTest(H5SingleGridTests('test_interpolation_newgrid'))
  return suite

class H5Tests(unittest.TestCase):

  def set_random(self,f,npts):
    dset = f.create_dataset('random',(npts,),'f')
    dset[...]=1.
    dset[:] = np.random.rand(npts)
    dsum=np.sum(dset)
    return dsum

  def set_random2d(self,f,n0,n1):
    dset = f.create_dataset('random',(n0,n1),'f',chunks=(1,n1))
    dset[:,:] = np.random.rand(n0,n1)


  def test_Persistance(self):
    npts=100
    filename='randomtest2.hdf5'
    f=h5py.File(filename,'w')
    dsum = self.set_random(f,npts)
    dset=f['random']
    self.assertAlmostEqual(dsum,np.sum(dset))
    f.close()
    os.remove(filename)

  def test_MemoryPerformance(self):

    nb=32*24*12
    nt=1000
    filename='memorytest.hdf5'
    f=h5py.File(filename,'w')
    self.set_random2d(f,nb,nt)
    bigdata=f['random']
    maxdata=f.create_dataset('maxdata',(nt,),'f')
    imax=f.create_dataset('imax',(nt,),'i')
    t_slice=int(5e5/nb)

    n_slices=nt/t_slice
    for i in range(n_slices):
      imax[i*t_slice:(i+1)*t_slice]=np.argmax(bigdata[:,i*t_slice:(i+1)*t_slice],0)
      maxdata[i*t_slice:(i+1)*t_slice]=np.max(bigdata[:,i*t_slice:(i+1)*t_slice],0)
    imax[n_slices*t_slice:nt]=np.argmax(bigdata[:,n_slices*t_slice:nt],0)
    maxdata[n_slices*t_slice:nt]=np.max(bigdata[:,n_slices*t_slice:nt],0)

    f.close()
    os.remove(filename)
    

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

    filename='test_nll.hdf5'
    sg=H5NllSingleGrid(filename,nll_name)
    self.assertEqual(info['ny'],sg.grid_info['ny'])
    self.assertAlmostEqual(info['orig_lat'],sg.grid_info['orig_lat'])
    self.assertEqual(info['station'],'FIU')
    self.assertEqual(info['station'],sg.grid_info['station'])
    self.assertEqual(np_buf.shape,sg.grid_data.shape)
    self.assertAlmostEqual(np_buf[b_index],sg.grid_data[b_index])

    del sg
    del buf
    os.remove(filename)

  def test_nll2hdf5(self):
    base_path=os.getenv('WAVELOC_PATH')
    nll_name=os.path.join(base_path,'test_data','test.time')
    h5_name="%s.hdf5"%nll_name
    
    if os.path.isfile(h5_name) : os.remove(h5_name)
    nll2hdf5(nll_name,h5_name)
    self.assertTrue(os.path.isfile(h5_name))

  def test_interpolation_ones(self):

    data=np.ones((100,200,50)).flatten()
    info={}
    info['nx']=100
    info['ny']=200
    info['nz']=50
    info['dx']=1.
    info['dy']=0.1
    info['dz']=0.5
    info['x_orig']=10.
    info['y_orig']=30.
    info['z_orig']=15.
    
    x=np.random.rand()*info['nx']*info['dx']+info['x_orig']
    y=np.random.rand()*info['ny']*info['dy']+info['y_orig']
    z=np.random.rand()*info['nz']*info['dz']+info['z_orig']

    # put this data in a hdf5 file
    filename='interpolate.hdf5'
    if os.path.isfile(filename): os.remove(filename)
    sg=H5SingleGrid(filename,data,info)
    interp=sg.value_at_point(x,y,z)

    self.assertAlmostEqual(interp,1.)

    del sg
    os.remove(filename)
    
  def test_interpolation_sinc(self):

    data=np.zeros((100,200,50))
    info={}
    info['nx']=100
    info['ny']=200
    info['nz']=50
    info['dx']=1.
    info['dy']=0.1
    info['dz']=0.5
    info['x_orig']=10.
    info['y_orig']=30.
    info['z_orig']=15.
    
    X=np.arange(info['nx'])*info['dx']+info['x_orig']
    Y=np.arange(info['ny'])*info['dy']+info['y_orig']
    Z=np.arange(info['nz'])*info['dz']+info['z_orig']

    xx,yy=np.meshgrid(X,Y)
    data1=np.sin(xx**2 + yy**2) / (xx**2 + yy**2)
    data[:,:,0]=data1.T
    

    x=np.random.rand()*info['nx']*info['dx']+info['x_orig']
    y=np.random.rand()*info['ny']*info['dy']+info['y_orig']
    z=info['z_orig']

    true_answer=np.sin(x**2 + y**2) / (x**2 + y**2)

    # put this data in a hdf5 file
    filename='interpolate.hdf5'
    if os.path.isfile(filename): os.remove(filename)
    sg=H5SingleGrid(filename,data.flatten(),info)
    interp=sg.value_at_point(x,y,z)

    self.assertAlmostEqual(interp,true_answer,2)

    del sg
    del data
    os.remove(filename)

  def test_interpolation_newgrid(self):

    data=np.zeros((100,200,50))
    info={}
    info['nx']=100
    info['ny']=200
    info['nz']=50
    info['dx']=1.
    info['dy']=0.1
    info['dz']=0.5
    info['x_orig']=10.
    info['y_orig']=30.
    info['z_orig']=15.
    
    new_info={}
    new_info['nx']=5
    new_info['ny']=10
    new_info['nz']=2
    new_info['dx']=1.5
    new_info['dy']=0.15
    new_info['dz']=0.55
    new_info['x_orig']=10.5
    new_info['y_orig']=30.5
    new_info['z_orig']=15

    X=np.arange(info['nx'])*info['dx']+info['x_orig']
    Y=np.arange(info['ny'])*info['dy']+info['y_orig']
    Z=np.arange(info['nz'])*info['dz']+info['z_orig']

    xx,yy=np.meshgrid(X,Y)
    data1=np.sin(xx**2 + yy**2) / (xx**2 + yy**2)
    data[:,:,0]=data1.T
    

    ix=np.random.randint(new_info['nx'])
    iy=np.random.randint(new_info['ny'])
    x=ix*new_info['dx']+new_info['x_orig']
    y=iy*new_info['dy']+new_info['y_orig']
    z=new_info['z_orig']

    true_answer=np.sin(x**2 + y**2) / (x**2 + y**2)

    # put this data in a hdf5 file
    filename='interpolate.hdf5'
    new_filename='interpolate_new.hdf5'
    if os.path.isfile(filename): os.remove(filename)
    sg=H5SingleGrid(filename,data.flatten(),info)
    new_sg=sg.interp_to_newgrid(new_filename,new_info)


    interp=sg.value_at_point(x,y,z)
    new_interp=new_sg.value_at_point(x,y,z)

    self.assertAlmostEqual(interp,true_answer,2)
    self.assertAlmostEqual(new_interp,true_answer,2)

    del sg
    del new_sg
    del data
    os.remove(filename)
    os.remove(new_filename)
 
if __name__ == '__main__':

  import logging
  logging.basicConfig(level=logging.INFO, format='%(levelname)s : %(asctime)s : %(message)s')
 
  unittest.TextTestRunner(verbosity=2).run(suite())
 
