"""
Helper functions for creating and manipulating unstructured grids.
"""

import os
import h5py
import numpy as np
from NLL_grid_lib import read_hdr_file
from sklearn import preprocessing
from sklearn.svm import SVR


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


def ugrid_closest_point_index(x, y, z, xi, yi, zi):
    """
    Returns the index of the closest point to xi, yi, zi.

    :param x: x-coordinates of the unstructured grid
    :param y: y-coordinates of the unstructured grid
    :param z: z-coordinates of the unstructured grid
    :param xi: x-corrdinate of point of interest
    :param yi: y-corrdinate of point of interest
    :param zi: z-corrdinate of point of interest
    """

    dist = (x-xi)**2 + (y-yi)**2 + (z-zi)**2
    ic = np.argmin(dist)    # index of closest point

    return ic


def select_points_closeto_plane(x, y, z, p_string, p_coord, p_dist):
    """
    Selects points within a certain distance of an x, y, or z plane.

    """
    if p_string == 'x':
        result = np.abs(x - p_coord) < p_dist
    elif p_string == 'y':
        result = np.abs(y - p_coord) < p_dist
    elif p_string == 'z':
        result = np.abs(z - p_coord) < p_dist
    else:
        raise UserWarning('Unknown plane %s' % p_string)

    return result


def ugrid_svr(x, y, z, values, xi, yi, zi, C=1.0, epsilon=0.1):
    """
    Uses support vector regression to estimate values at points xi, yi, zi when
    values at irregular points x, y, z are known.

    :param x: x-coordinates of the unstructured grid
    :param y: y-coordinates of the unstructured grid
    :param z: z-coordinates of the unstructured grid
    :param values: values at irregular points x, y, z
    :param xi: x-corrdinate of point of interest
    :param yi: y-corrdinate of point of interest
    :param zi: z-corrdinate of point of interest
    """

    X = np.array([x, y, z]).T
    Xi = np.array([xi, yi, zi]).T

    scaler = preprocessing.StandardScaler().fit(X)

    svr = SVR(C=C, epsilon=epsilon)
    yi = svr.fit(scaler.transform(X), values).predict(scaler.transform(Xi))

    return yi


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
        ix, iy, iz = np.unravel_index(ib, grid_shape)
        x[ib] = x_range[ix]
        y[ib] = y_range[iy]
        z[ib] = z_range[iz]

    return x, y, z
