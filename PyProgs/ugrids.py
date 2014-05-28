"""
Helper functions for creating and manipulating unstructured grids.
"""

import os
import h5py
import numpy as np
from NllGridLib import read_hdr_file


def create_random_ugrid(xmin, xmax, ymin, ymax, zmin, zmax, npts):
    """
    Creates a random unstructured grid for use with waveloc, by sampling with a
    uniform probability distribution.

    :param xmin: Minimum x-coordinate.
    :param xmax: Maximum x-coordinate.
    :param ymin: Minimum y-coordinate.
    :param ymax: Maximum y-coordinate.
    :param zmin: Minimum z-coordinate.
    :param zmax: Maximum z-coordinate.
    :param npts: Number of points
    :param filename: Filename for HDF5 output

    :type xmin: float
    :type xmax: float
    :type ymin: float
    :type ymax: float
    :type zmin: float
    :type zmax: float
    :type npts: integer
    :type filename: string
    """

    # create random points
    x = np.random.uniform(low=xmin, high=xmax, size=npts)
    y = np.random.uniform(low=ymin, high=ymax, size=npts)
    z = np.random.uniform(low=zmin, high=zmax, size=npts)

    return x, y, z


def write_ugrid(x, y, z, filename):
    """
    Writes a ugrid HDF5 file given x, y, z vectors of points

    :param x: x-coordinates of points
    :param y: y-coordinates of points
    :param z: z-coordinates of points
    :param filname: File to write
    """

    # sanity check
    if not (len(x) == len(y) and len(y) == len(z)):
        msg = 'x, y, z arrays of different lengths'
        raise UserWarning(msg)

    npts = len(x)

    # write to HDF5 file
    f = h5py.File(filename, 'w')
    xd = f.create_dataset('x', data=x, compression='lzf')
    yd = f.create_dataset('y', data=y, compression='lzf')
    zd = f.create_dataset('z', data=z, compression='lzf')
    xd.attrs['npts'] = npts
    yd.attrs['npts'] = npts
    zd.attrs['npts'] = npts

    f.close()


def read_ugrid(filename):
    """
    Reads an unstructured grid

    :param filename: File to read
    :type filename: string

    :rtype: numpy array
    :returns: x, y, z
    """

    if not os.path.isfile(filename):
        msg = 'File %s does not exist' % filename
        raise UserWarning(msg)

    f = h5py.File(filename, 'r')
    xd = f['x']
    yd = f['y']
    zd = f['z']
    npts = xd.attrs['npts']

    x = np.empty(npts, dtype='float')
    y = np.empty(npts, dtype='float')
    z = np.empty(npts, dtype='float')

    x = xd[:]
    y = yd[:]
    z = zd[:]

    f.close()

    return x, y, z


def nll2random_ugrid(nll_filename, npts):
    """
    Reads a NLL .hdr file, and creates a random unstructured-grid by sampling
    the extent of the NonLinLoc grid unifomly using npts points.

    :param nll_filename: NonLinLoc filename
    :param ugrid_filename: Filename for output unstructured-grid

    :type nll_filename: string
    :type ugrid_filename: string

    """

    info = read_hdr_file(nll_filename)

    xmin = info['x_orig']
    ymin = info['y_orig']
    zmin = info['z_orig']

    xmax = xmin+info['nx']*info['dx']
    ymax = ymin+info['ny']*info['dy']
    zmax = zmin+info['nz']*info['dz']

    x, y, z = create_random_ugrid(xmin, xmax, ymin, ymax, zmin, zmax, npts)
    return x, y, z


def nll2reg_ugrid(nll_filename):
    """
    Reads a NLL .hdr file, and creates an unstructured grid that is equivalent
    to the regular NLL grid.

    :param nll_filename: NonLinLoc filename
    :param ugrid_filename: Filename for output unstructured-grid

    :type nll_filename: string
    :type ugrid_filename: string

    """

    info = read_hdr_file(nll_filename)

    npts = info['nx']*info['ny']*info['nz']

    xmin = info['x_orig']
    ymin = info['y_orig']
    zmin = info['z_orig']

    xmax = xmin+info['nx']*info['dx']
    ymax = ymin+info['ny']*info['dy']
    zmax = zmin+info['nz']*info['dz']

    grid_shape = (info['nx'], info['ny'], info['nz'])

    # get the ranges
    x_range = np.arange(xmin, xmax, info['dx'])
    y_range = np.arange(ymin, ymax, info['dy'])
    z_range = np.arange(zmin, zmax, info['dz'])

    # create space for the points
    x = np.empty(npts, dtype='float')
    y = np.empty(npts, dtype='float')
    z = np.empty(npts, dtype='float')

    # fill in the points
    for ib in xrange(npts):
        ix, iy, iz = np.unravel_index(ib, grid_shape )
        x[ib] = x_range[ix]
        y[ib] = y_range[iy]
        z[ib] = z_range[iz]

    return x, y, z
