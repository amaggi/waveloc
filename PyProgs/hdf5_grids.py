import h5py,os
import numpy as np
from NllGridLib import read_hdr_file

class H5SingleGrid(object):

  def __init__(self, filename=None, grid_data=None, grid_info=None):

    if os.path.isfile(filename):
      self._f=h5py.File(filename,'r')
      self.grid_data=self._f['grid_data']
      self.grid_info=self.grid_data.attrs

    else:
      self._f=h5py.File(filename,'w')

      if not grid_data==None:
        self.grid_data=self._f.create_dataset('grid_data',data=grid_data,compression='lzf')
      else : self.grid_data=None

      if not grid_info==None:
        self.grid_info=self.grid_data.attrs
        for key,value in grid_info.iteritems():
          self.grid_info[key]=value
      else : self.grid_info = None
      
  def __del__(self):
    self._f.close()

      
class H5NllSingleGrid(H5SingleGrid):

  def __init__(self,filename,nll_filename):

    from array import array

    hdr="%s.hdr"%nll_filename
    buf="%s.buf"%nll_filename

    info=read_hdr_file(hdr)
    nx=info['nx']
    ny=info['ny']
    nz=info['nz']

    f=open(buf,'rb')
    buf=array('f')
    buf.fromfile(f,nx*ny*nz)
    f.close() 

    H5SingleGrid.__init__(self,filename,buf,info)
    

if __name__=='__main__' : 
  
  pass
