"""
Helper functions for creating and manipulating unstructured grids.
"""

import os
import h5py
import numpy as np


def create_random_ugrid(xmin, xmax, ymin, ymax, zmin, zmax, npts, filename):
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

    # write to HDF5 file
    f = h5py.File(filename, 'w')
    xd = f.create_dataset('x', data=x)
    yd = f.create_dataset('y', data=y)
    zd = f.create_dataset('z', data=z)
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


def nll2ugrid(nll_filename, ugrid_filename, npts):
    """
    Reads a NLL .hdr file, and creates a random unstructured-grid by sampling
    the extent of the NonLinLoc grid unifomly using npts points.

    :param nll_filename: NonLinLoc filename
    :param ugrid_filename: Filename for output unstructured-grid

    :type nll_filename: string
    :type ugrid_filename: string

    """
    from NllGridLib import read_hdr_file

    info = read_hdr_file(nll_filename)

    xmin = info['x_orig']
    ymin = info['y_orig']
    zmin = info['z_orig']

    xmax = xmin+info['nx']*info['dx']
    ymax = ymin+info['ny']*info['dy']
    zmax = zmin+info['nz']*info['dz']

    create_random_ugrid(xmin, xmax, ymin, ymax, zmin, zmax, npts,
                        ugrid_filename)
