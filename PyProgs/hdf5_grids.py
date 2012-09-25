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

  def value_at_point(self,x,y,z,epsilon=0.001):
    """
    Performs n-linear interpolation on the regular grid
    """
    nx=self.grid_info['nx']
    ny=self.grid_info['ny']
    nz=self.grid_info['nz']
    x_orig=self.grid_info['x_orig']
    y_orig=self.grid_info['y_orig']
    z_orig=self.grid_info['z_orig']
    dx=self.grid_info['dx']
    dy=self.grid_info['dy']
    dz=self.grid_info['dz']

    min_x=x_orig
    max_x=x_orig+nx*dx
    min_y=y_orig
    max_y=y_orig+ny*dy
    min_z=z_orig
    max_z=z_orig+nz*dz

    grid_size=nx*nz*nz

    # sanity check for point being within grid
    # use a tollerance value to avoid problems with numerical errors
    if x < min_x-epsilon or x > max_x+epsilon \
    or y < min_y-epsilon or y > max_y+epsilon \
    or z < min_z-epsilon or z > max_z+epsilon :
      raise UserWarning('Point (%f, %f, %f) is outside the grid (tollerance=%f).'%(x,y,z,epsilon))

    # fix up lower and upper bounds if they are still (just) outside the grid
    if x < min_x : x=min_x
    if y < min_y : y=min_y
    if z < min_z : z=min_z

    if x > max_x : x=max_x
    if y > max_y : y=max_y
    if z > max_z : z=max_z
    
    # make arrays of X, Y and Z ranges 
    X=np.arange(nx)*dx+x_orig
    Y=np.arange(ny)*dy+y_orig
    Z=np.arange(nz)*dz+z_orig

    # get the position this point would have in the X,Y,Z arrays if they were extended by 1
    ix=X.searchsorted(x)
    iy=Y.searchsorted(y)
    iz=Z.searchsorted(z)

    # set the interpolation "box" for extreme cases
    if nx==1 : # special case of 2D grid
      ix1=0
      ix2=0
    elif ix==0: # lower bound
      ix1=0
      ix2=1
    elif ix==nx: # upper bound
      ix1=nx-2
      ix2=nx-1
    else :	# general case
      ix1=ix-1
      ix2=ix

    if iy==0:	# lower bound
      iy1=0
      iy2=1
    elif iy==ny: # upper bound
      iy1=ny-2
      iy2=ny-1
    else :	# general case
      iy1=iy-1
      iy2=iy

    if iz==0:	# lower bound
      iz1=0
      iz2=1
    elif iz==nz: # upper bound
      iz1=nz-2
      iz2=nz-1
    else :	# general case
      iz1=iz-1
      iz2=iz

    # set up the values
    # bottom four values counterclockwise from x1y1
    v_x1y1z1=self.grid_data[np.ravel_multi_index((ix1,iy1,iz1),(nx,ny,nz))]
    v_x2y1z1=self.grid_data[np.ravel_multi_index((ix2,iy1,iz1),(nx,ny,nz))]
    v_x2y2z1=self.grid_data[np.ravel_multi_index((ix2,iy2,iz1),(nx,ny,nz))]
    v_x1y2z1=self.grid_data[np.ravel_multi_index((ix1,iy2,iz1),(nx,ny,nz))]
    # top four values counterclockwise from x1y1
    v_x1y1z2=self.grid_data[np.ravel_multi_index((ix1,iy1,iz2),(nx,ny,nz))]
    v_x2y1z2=self.grid_data[np.ravel_multi_index((ix2,iy1,iz2),(nx,ny,nz))]
    v_x2y2z2=self.grid_data[np.ravel_multi_index((ix2,iy2,iz2),(nx,ny,nz))]
    v_x1y2z2=self.grid_data[np.ravel_multi_index((ix1,iy2,iz2),(nx,ny,nz))]

    # set up interpolators
    # take extra care over the X interpolator in case of 2D grid
    if ix2==ix1:
      tx=0
    else:
      tx=(x-X[ix1])/(X[ix2]-X[ix1])
    ty=(y-Y[iy1])/(Y[iy2]-Y[iy1])
    tz=(z-Z[iz1])/(Z[iz2]-Z[iz1])

    # do bilinear interpolation
    result = (1-tx) * (1-ty) * (1-tz) * v_x1y1z1 + \
                tx  * (1-ty) * (1-tz) * v_x2y1z1 + \
                tx  *    ty  * (1-tz) * v_x2y2z1 + \
             (1-tx) *    ty  * (1-tz) * v_x1y2z1 + \
             (1-tx) * (1-ty) *    tz  * v_x1y1z2 + \
                tx  * (1-ty) *    tz  * v_x2y1z2 + \
                tx  *    ty  *    tz  * v_x2y2z2 + \
             (1-tx) *    ty  *    tz  * v_x1y2z2 

    return result

  #@profile
  def interp_to_newgrid(self,new_filename,new_grid_info):

    nx=new_grid_info['nx']
    ny=new_grid_info['ny']
    nz=new_grid_info['nz']

    dx=new_grid_info['dx']
    dy=new_grid_info['dy']
    dz=new_grid_info['dz']

    x_orig=new_grid_info['x_orig']
    y_orig=new_grid_info['y_orig']
    z_orig=new_grid_info['z_orig']
    
    # if you're calling this function, you want any existing file overwritten
    f=h5py.File(new_filename,'w')
    buf=f.create_dataset('grid_data',(nx*ny*nz,),'f')
    for key,value in new_grid_info.iteritems():
      buf.attrs[key]=value

    #initialize new buffer

    new_x=np.arange(nx)*dx+x_orig
    new_y=np.arange(ny)*dy+y_orig
    new_z=np.arange(nz)*dz+z_orig

    # loop doing interpolation
    for ix in range(nx):
     x=new_x[ix]
     for iy in range(ny):
       y=new_y[iy]
       for iz in range(nz):
         z=new_z[iz]
         buf[np.ravel_multi_index((ix,iy,iz),(nx,ny,nz))]=self.value_at_point(x,y,z)

    f.close()

    # create the new Grid file and object
    new_grid=H5SingleGrid(new_filename)
    return new_grid


      
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
    
####
# helper functions
####

def nll2hdf5(nll_name,h5_name):
  h5=H5NllSingleGrid(h5_name,nll_name)
  del h5

if __name__=='__main__' : 
  
  pass
