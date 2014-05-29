import h5py
import os
import logging
import numpy as np
from scipy import ndimage
from NllGridLib import read_hdr_file
from ugrids import read_ugrid, nll2reg_ugrid

"""
Contains wrapper classes to easily manipulate hdf5 files.
"""


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

    grid_data = None
    grid_info = None
    _f = None

    def __init__(self, filename=None, grid_data=None, grid_info=None):
        """
        Initialises an instance of HDF5SingleGrid.

        :param filename: If a file of this name exists, it is opened for
            reading.  Attributes grid_data and grid_info will then refer  to
            the HDF5 dataset named 'grid_data' and to its HDF5  attributes.  If
            the file does not already exist, then a new  HDF5 file is opened
            for writing, the grid_data paramter is  used to initialise the
            'grid_data' HDF5 dataset, and the  grid_info parameter is used to
            initialise its attributes.
        :param grid_data: An array containing the grid data to be written to
            the file.
        :param grid_info: Attributes to be used for the 'grid_data' dataset.
        """

        if os.path.isfile(filename):
            self._f = h5py.File(filename, 'r')
            self.grid_data = self._f['grid_data']
            self.grid_info = self.grid_data.attrs

        else:
            self._f = h5py.File(filename, 'w')

            if not grid_data is None:
                self.grid_data = self._f.create_dataset('grid_data',
                                                        data=grid_data,
                                                        compression='lzf')

            if not grid_info is None:
                self.grid_info = self.grid_data.attrs
                for key, value in grid_info.iteritems():
                    self.grid_info[key] = value

    def __del__(self):
        self._f.close()

    def value_at_point(self, x, y, z):
        """
        Performs 3D interpolation for a single point. Calls value at points.

        :param x: x-coordinate of the point of interest
        :param y: y-coordinate of the point of interest
        :param z: z-coordinate of the point of interest

        :type x: float
        :type y: float
        :type z: float

        :rtype: float
        :returns: value at point

        """

        x_array = np.array([x, x])
        y_array = np.array([y, y])
        z_array = np.array([z, z])

        result = self.value_at_points(x_array, y_array, z_array)
        return result[0]

    def value_at_points(self, x, y, z):
        """
        Performs 3D interpolation on the regular grid.
        Uses scipy.ndimage. Works on numpy arrays of points

        :param x: x-coordinates of points
        :param y: y-coordinates of points
        :param z: z-coordinates of points

        :type x: numpy array
        :type y: numpy array
        :type z: numpy array

        :returns: values at points

        """

        grid_data = np.empty(self.grid_data.shape, dtype=float)
        grid_data[:] = self.grid_data[:]
        grid_data = grid_data.reshape(self.grid_info['nx'],
                                      self.grid_info['ny'],
                                      self.grid_info['nz'])

        ix = (x - self.grid_info['x_orig']) / self.grid_info['dx']
        iy = (y - self.grid_info['y_orig']) / self.grid_info['dy']
        iz = (z - self.grid_info['z_orig']) / self.grid_info['dz']

        coords = np.array([ix, iy, iz])

        result = ndimage.map_coordinates(grid_data, coords, order=3,
                                         mode='nearest')

        return result

    def interp_to_newgrid(self, new_filename, new_grid_info):
        """
        Intrepolates to a new grid.

        :param new_filename: Filename for new grid
        :param new_grid_info: Grid info dictionary for new grid

        """

        nx = new_grid_info['nx']
        ny = new_grid_info['ny']
        nz = new_grid_info['nz']

        dx = new_grid_info['dx']
        dy = new_grid_info['dy']
        dz = new_grid_info['dz']

        x_orig = new_grid_info['x_orig']
        y_orig = new_grid_info['y_orig']
        z_orig = new_grid_info['z_orig']

        # if you're calling this function, you want any existing file
        # overwritten
        f = h5py.File(new_filename, 'w')
        buf = f.create_dataset('grid_data', (nx*ny*nz, ), 'f')
        for key, value in new_grid_info.iteritems():
            buf.attrs[key] = value

        #set coordinates for interpolation
        npts = nx*ny*nz
        new_shape = (nx, ny, nz)
        x = np.empty(npts, dtype=float)
        y = np.empty(npts, dtype=float)
        z = np.empty(npts, dtype=float)

        for i in xrange(npts):
            ix, iy, iz = np.unravel_index(i, (new_shape))
            x[i] = x_orig+ix*dx
            y[i] = y_orig+iy*dy
            z[i] = z_orig+iz*dz

        # do interpolation
        buf[:] = self.value_at_points(x, y, z)

        # close the old grid file
        f.close()

        # create the new Grid file and object
        new_grid = H5SingleGrid(new_filename)
        return new_grid


class H5NllSingleGrid(H5SingleGrid):
    """
    Class derived from H5SingleGrid that upon initialization reads a NLL pair
    of files (.hdr and .buf) and writes an equivalent hdf5 file.

    :param nll_name: NLL filename without extension
    :param h5_name: hdf5 filename
    """

    def __init__(self, filename, nll_filename):

        from array import array

        hdr = "%s.hdr" % nll_filename
        buf = "%s.buf" % nll_filename

        info = read_hdr_file(hdr)
        nx = info['nx']
        ny = info['ny']
        nz = info['nz']

        f = open(buf, 'rb')
        buf = array('f')
        buf.fromfile(f, nx*ny*nz)
        f.close()

        H5SingleGrid.__init__(self, filename, buf, info)


####
# HDF5 enabled functions
####


def nll2hdf5(nll_name, h5_name):
    """
    Translates NLL files to hdf5 format.

    :param nll_name: NLL filename without extension
    :param h5_name: hdf5 filename
    """

    h5 = H5NllSingleGrid(h5_name, nll_name)
    del h5


def interpolateTimeGrid(tgrid_file, x, y, z):
    """
    Interpolates a time_grid.hdf5 file to another file containing only the
    travel-times for the points in x, y, z.

    :param tgrid_file: HDF5 file containing the time-grid
    :param x: x-coordinates for points of interest
    :param y: y-coordinates for points of interest
    :param z: z-coordinates for points of interest

    """

    # sanity check for length of coordinate arrays
    if (len(x) != len(y)) or (len(y) != len(z)):
        msg = 'Coordinate arrays of different lengths.'
        raise ValueError(msg)

    # read the file to be interpolated
    time_grid = H5SingleGrid(tgrid_file)

    # do the interpolation
    ttimes = time_grid.value_at_points(x, y, z)
    return ttimes


def write_interpolated_time_ugrid(ttimes, sta, out_file):
    """
    Write time grid in an unstructured-grid format.

    :param ttimes: time ugrid to write
    :param sta: station name
    :param out_file: HDF5 file for output
    """
    # create the file for output
    f = h5py.File(out_file, 'w')
    buf = f.create_dataset('ttimes', data=ttimes, compression='lzf')
    buf.attrs['station'] = sta
    f.close()


def read_interpolated_time_ugrid(filename):
    """
    Read a time grid in an unstructured-grid format. The grid geometry itself
    must be read from the unstructured-grid x, y, z file.

    :param filename: HDF5 file to read

    :returns: sta, ttimes
    """
    # create the file for output
    f = h5py.File(filename, 'r')
    buf = f['ttimes']
    ttimes = buf[:]
    sta = buf.attrs['station']
    f.close()

    return sta, ttimes


def get_interpolated_time_ugrids(opdict):
    """
    Interpolates the NLL time grids onto an given ugrid. Uses options
    contained in a WavelocOptions.opdict.

    :param opdict: Dictionary of options in WavelocOptions.opdict format
    """
    import glob

    base_path = opdict['base_path']
    ugrid_type = opdict['ugrid_type']

    # read or create the x, y, z coordinates of search points
    if ugrid_type == 'USER':
        ugrid_file = opdict['ugrid_file']
        filename = os.path.join(base_path, 'lib', ugrid_file)
        x, y, z = read_ugrid(filename)
    elif ugrid_type == 'FULL':
        search_grid = opdict['search_grid']
        filename = os.path.join(base_path, 'lib', search_grid)
        x, y, z = nll2reg_ugrid(filename)

    # read full time grids

    full_time_grids = glob.glob(os.path.join(base_path, 'lib',
                                opdict['time_grid']+'*.hdf5'))
    full_time_grids.sort()

    # create the dictionary of interpolated time_grids
    time_grids = {}
    tgrid_dir = os.path.join(base_path, 'out', opdict['outdir'], 'time_grids')

    # for each of the full-length time grids
    logging.info('Loading time grids ... ')
    for f_timegrid in full_time_grids:
        # get the filename of the corresponding short-length grid (the one for
        # the search grid in particular)
        f_basename = os.path.basename(f_timegrid)
        tgrid_filename = os.path.join(tgrid_dir, f_basename)

        # if file exists and we want to load it, then open the file and give it
        # to the dictionary
        if os.path.isfile(tgrid_filename) and opdict['load_ttimes_buf']:
            logging.debug('Loading %s' % tgrid_filename)
            sta, ttimes = read_interpolated_time_ugrid(tgrid_filename)

        # otherwise create / overwrite it
        else:
            logging.info('Creating %s - Please be patient' % tgrid_filename)
            ttimes = interpolateTimeGrid(f_timegrid, x, y, z)
            # get station name
            full_grid = H5SingleGrid(f_timegrid)
            sta = full_grid.grid_info['station']
            del full_grid
            # write new file
            write_interpolated_time_ugrid(ttimes, sta, tgrid_filename)

        # add to dictionary
        time_grids[sta] = ttimes

    return x, y, z, time_grids


def get_interpolated_time_grids(opdict):
    """
    Interpolates the NLL time grids onto the search grids. Uses options
    contained in a WavelocOptions.opdict.

    :param opdict: Dictionary of options in WavelocOptions.opdict format
    """
    import glob
    from NllGridLib import read_hdr_file

    base_path = opdict['base_path']
    full_time_grids = glob.glob(os.path.join(base_path, 'lib',
                                opdict['time_grid']+'*.hdf5'))
    full_time_grids.sort()
    if len(full_time_grids) == 0:
        raise UserWarning('No .hdf5 time grids found in directory %s' %
                          (os.path.join(base_path, 'lib')))

    # read the search grid
    search_grid = os.path.join(base_path, 'lib', opdict['search_grid'])
    tgrid_dir = os.path.join(base_path, 'out', opdict['outdir'], 'time_grids')
    if not os.path.exists(tgrid_dir):
        os.makedirs(tgrid_dir)
    search_info = read_hdr_file(search_grid)

    time_grids = {}

    # for each of the full-length time grids
    logging.info('Loading time grids ... ')
    for f_timegrid in full_time_grids:
        f_basename = os.path.basename(f_timegrid)
        # get the filename of the corresponding short-length grid (the one for
        # the search grid in particular)
        tgrid_filename = os.path.join(tgrid_dir, f_basename)

        # if file exists and we want to load it, then open the file and give it
        # to the dictionary
        if os.path.isfile(tgrid_filename) and opdict['load_ttimes_buf']:
            logging.debug('Loading %s' % tgrid_filename)
            grid = H5SingleGrid(tgrid_filename)
            name = grid.grid_info['station']
            time_grids[name] = grid.grid_data[:]
            del grid

        # if the file does not exist, or want to force re-creation, then create
        # it
        if not os.path.isfile(tgrid_filename) or not opdict['load_ttimes_buf']:
            logging.info('Creating %s - Please be patient' % tgrid_filename)
            full_grid = H5SingleGrid(f_timegrid)
            # copy the common part of the grid info
            new_info = {}
            for name, value in full_grid.grid_info.iteritems():
                new_info[name] = value
            # set the new part of the grid info to correspond to the search
            # grid
            new_info['x_orig'] = search_info['x_orig']
            new_info['y_orig'] = search_info['y_orig']
            new_info['z_orig'] = search_info['z_orig']
            new_info['nx'] = search_info['nx']
            new_info['ny'] = search_info['ny']
            new_info['nz'] = search_info['nz']
            new_info['dx'] = search_info['dx']
            new_info['dy'] = search_info['dy']
            new_info['dz'] = search_info['dz']
            # do interpolation
            grid = full_grid.interp_to_newgrid(tgrid_filename, new_info)
            # add to dictionary
            name = grid.grid_info['station']
            time_grids[name] = grid.grid_data[:]
            # close full grid safely
            del grid
            del full_grid

    # x, y, z ranges
    x_range = np.arange(0, search_info['nx']*search_info['dx'],
                        search_info['dx']) + search_info['x_orig']
    y_range = np.arange(0, search_info['ny']*search_info['dy'],
                        search_info['dy']) + search_info['y_orig']
    z_range = np.arange(0, search_info['nz']*search_info['dz'],
                        search_info['dz']) + search_info['z_orig']

    npts = search_info['nx']*search_info['ny']*search_info['nz']
    grid_shape = (search_info['nx'], search_info['ny'], search_info['nz'])

    # create space for the points
    x = np.empty(npts, dtype='float')
    y = np.empty(npts, dtype='float')
    z = np.empty(npts, dtype='float')

    # fill in the points
    for ib in xrange(npts):
        ix, iy, iz = np.unravel_index(ib, grid_shape)
        x[ib] = x_range[ix]
        y[ib] = y_range[iy]
        z[ib] = z_range[iz]

    return x, y, z, time_grids


if __name__ == '__main__':

    pass
