import h5py,os,logging, tempfile
import numpy as np
from time import time
from NllGridLib import read_hdr_file

class H5SingleGrid(object):
  """
  Class that wraps several methods for regular grids in HDF5 format.

  **Attributes**

  .. attribute:: grid_data

      An HDF5 dataset.

  .. attribute:: grid_info

      Attributes of an HDF5 dataset.

  **Methods**
  """

  def __init__(self, filename=None, grid_data=None, grid_info=None):
    """
    Initialises an instance of HDF5SingleGrid.
      
    :param filename: If a file of this name exists, it is opened for
    reading.  Attributes grid_data and grid_info will then refer  to the HDF5
    dataset named 'grid_data' and to its HDF5  attributes.  If the file does
    not already exist, then a new  HDF5 file is opened for writing, the
    grid_data paramter is  used to initialise the 'grid_data' HDF5 dataset,
    and the  grid_info parameter is used to initialise its attributes.
    :param grid_data: An array containing the grid data to be written to the
    file.
    :param grid_info: Attributes to be used for the 'grid_data' dataset.
    """

    if os.path.isfile(filename):
      self._f=h5py.File(filename,'r')
      self.grid_data=self._f['grid_data']
      self.grid_info=self.grid_data.attrs

    else:
      self._f=h5py.File(filename,'w')

      if not grid_data==None:
        self.grid_data=self._f.create_dataset('grid_data',data=grid_data,\
		compression='lzf')
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
    Performs n-linear interpolation on the regular grid.
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
    if nx==1 : # special case of 2D grid in x
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

    if ny==1 : # special case of 2D grid in y
      iy1=0
      iy2=0
    elif iy==0:	# lower bound
      iy1=0
      iy2=1
    elif iy==ny: # upper bound
      iy1=ny-2
      iy2=ny-1
    else :	# general case
      iy1=iy-1
      iy2=iy

    if nz==1 : # special case of 2D grid in y
      iz1=0
      iz2=0
    elif iz==0:	# lower bound
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
    # take extra care over the interpolators in case of 2D grids
    if ix2==ix1:
      tx=0
    else:
      tx=(x-X[ix1])/(X[ix2]-X[ix1])
    if iy2==iy1:
      ty=0
    else:
      ty=(y-Y[iy1])/(Y[iy2]-Y[iy1])
    if iz2==iz1:
      tz=0
    else:
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
    for ix in xrange(nx):
     x=new_x[ix]
     for iy in xrange(ny):
       y=new_y[iy]
       for iz in xrange(nz):
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
# HDF5 enabled functions
####


def nll2hdf5(nll_name,h5_name):
  h5=H5NllSingleGrid(h5_name,nll_name)
  del h5

def get_interpolated_time_grids(opdict):
  import glob
  from NllGridLib import read_hdr_file

  base_path=opdict['base_path']
  full_time_grids=glob.glob(os.path.join(base_path,'lib',opdict['time_grid']+'*.hdf5'))
  full_time_grids.sort()
  if len(full_time_grids)==0 : raise UserWarning('No .hdf5 time grids found in directory %s'%(os.path.join(base_path,'lib')))

  # read the search grid
  search_grid=os.path.join(base_path,'lib',opdict['search_grid'])
  tgrid_dir=os.path.join(base_path,'out',opdict['outdir'],'time_grids')
  if not os.path.exists(tgrid_dir) : os.makedirs(tgrid_dir)
  search_info=read_hdr_file(search_grid)
  
  time_grids={}

  # for each of the full-length time grids
  logging.info('Loading time grids ... ')
  for f_timegrid in full_time_grids:
    f_basename=os.path.basename(f_timegrid)
    # get the filename of the corresponding short-length grid (the one for the search grid in particular)
    tgrid_filename=os.path.join(tgrid_dir,f_basename)

    # if file exists and we want to load it, then open the file and give it to the dictionary
    if os.path.isfile(tgrid_filename) and opdict['load_ttimes_buf']:    
      logging.debug('Loading %s'%tgrid_filename)
      grid=H5SingleGrid(tgrid_filename)
      name=grid.grid_info['station']
      time_grids[name]=grid

    # if the file does not exist, or want to force re-creation, then create it 
    if not os.path.isfile(tgrid_filename) or not opdict['load_ttimes_buf']:    
      logging.info('Creating %s - Please be patient'%tgrid_filename)
      full_grid=H5SingleGrid(f_timegrid)
      # copy the common part of the grid info
      new_info={}
      for name,value in full_grid.grid_info.iteritems():
        new_info[name]=value
      # set the new part of the grid info to correspond to the search grid
      new_info['x_orig']=search_info['x_orig']
      new_info['y_orig']=search_info['y_orig']
      new_info['z_orig']=search_info['z_orig']
      new_info['nx']=search_info['nx']
      new_info['ny']=search_info['ny']
      new_info['nz']=search_info['nz']
      new_info['dx']=search_info['dx']
      new_info['dy']=search_info['dy']
      new_info['dz']=search_info['dz']
      # do interpolation
      grid=full_grid.interp_to_newgrid(tgrid_filename,new_info)
      # add to dictionary
      name=grid.grid_info['station']
      time_grids[name]=grid
      # close full grid safely
      del full_grid

  return time_grids


if __name__=='__main__' : 
  
  pass
