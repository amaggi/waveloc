import os
import unittest
import numpy as np
from NllGridLib import read_hdr_file
from ugrids import create_random_ugrid, read_ugrid, nll2random_ugrid,\
    nll2reg_ugrid


def suite():

    suite = unittest.TestSuite()

    suite.addTest(UgridTests('test_create_read_random'))
    suite.addTest(UgridTests('test_nll2random_ugrid'))
    suite.addTest(UgridTests('test_nll2reg_ugrid'))

    return suite


class UgridTests(unittest.TestCase):

    def test_create_read_random(self):

        filename = 'random_ugrid.hdf5'
        npts = 1000

        xmin, ymin, zmin = np.random.uniform(low=0., high=10., size=3)
        xmax, ymax, zmax = np.random.uniform(low=20., high=50., size=3)

        create_random_ugrid(xmin, xmax, ymin, ymax, zmin, zmax, npts, filename)
        x, y, z = read_ugrid(filename)
        os.remove(filename)

        self.assertEqual(len(x), npts)
        self.assertTrue(np.min(x) >= xmin)
        self.assertTrue(np.max(x) < xmax)
        self.assertTrue(np.min(y) >= ymin)
        self.assertTrue(np.max(y) < ymax)
        self.assertTrue(np.min(z) >= zmin)
        self.assertTrue(np.max(z) < zmax)

    def test_nll2random_ugrid(self):

        npts = 1000

        base_path = os.getenv('WAVELOC_PATH')
        nll_filename = os.path.join(base_path, 'test_data',
                                    'test_grid.search.hdr')
        ugrid_filename = 'test_ugrid.hdf5'

        nll2random_ugrid(nll_filename, ugrid_filename, npts)
        x, y, z = read_ugrid(ugrid_filename)
        os.remove(ugrid_filename)

        info = read_hdr_file(nll_filename)
        xmin = info['x_orig']
        ymin = info['y_orig']
        zmin = info['z_orig']

        xmax = xmin+info['nx']*info['dx']
        ymax = ymin+info['ny']*info['dy']
        zmax = zmin+info['nz']*info['dz']

        self.assertEqual(len(x), npts)
        self.assertTrue(np.min(x) >= xmin)
        self.assertTrue(np.max(x) < xmax)
        self.assertTrue(np.min(y) >= ymin)
        self.assertTrue(np.max(y) < ymax)
        self.assertTrue(np.min(z) >= zmin)
        self.assertTrue(np.max(z) < zmax)

    def test_nll2reg_ugrid(self):

        npts = 1000

        base_path = os.getenv('WAVELOC_PATH')
        nll_filename = os.path.join(base_path, 'test_data',
                                    'test_grid.search.hdr')
        ugrid_filename = 'test_ugrid.hdf5'

        nll2reg_ugrid(nll_filename, ugrid_filename)
        x, y, z = read_ugrid(ugrid_filename)
        os.remove(ugrid_filename)

        info = read_hdr_file(nll_filename)
        nx = info['nx']
        ny = info['ny']
        nz = info['nz']
        grid_shape = (nx, ny, nz)

        ib = np.random.randint(0, nx*ny*nz)
        ix, iy, iz = np.unravel_index(ib, grid_shape)

        xval_exp = info['x_orig']+ix*info['dx']
        yval_exp = info['y_orig']+iy*info['dy']
        zval_exp = info['z_orig']+iz*info['dz']

        self.assertAlmostEqual(x[ib], xval_exp)
        self.assertAlmostEqual(y[ib], yval_exp)
        self.assertAlmostEqual(z[ib], zval_exp)


if __name__ == '__main__':

    unittest.TextTestRunner(verbosity=2).run(suite())
