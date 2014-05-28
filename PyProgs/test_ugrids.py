import os
import unittest
import numpy as np
from NllGridLib import read_hdr_file
from ugrids import create_random_ugrid, nll2random_ugrid, nll2reg_ugrid,\
    ugrid_closest_point_index


def suite():

    suite = unittest.TestSuite()

    suite.addTest(UgridTests('test_create_read_random'))
    suite.addTest(UgridTests('test_nll2random_ugrid'))
    suite.addTest(UgridTests('test_nll2reg_ugrid'))
    suite.addTest(UgridTests('test_closest_point'))

    return suite


class UgridTests(unittest.TestCase):

    def test_create_read_random(self):

        npts = 1000

        xmin, ymin, zmin = np.random.uniform(low=0., high=10., size=3)
        xmax, ymax, zmax = np.random.uniform(low=20., high=50., size=3)

        x, y, z = create_random_ugrid(xmin, xmax, ymin, ymax, zmin, zmax, npts)

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

        x, y, z = nll2random_ugrid(nll_filename, npts)

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

        base_path = os.getenv('WAVELOC_PATH')
        nll_filename = os.path.join(base_path, 'test_data',
                                    'test_grid.search.hdr')

        x, y, z = nll2reg_ugrid(nll_filename)

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

    def test_closest_point(self):

        npts = 1000
        x, y, z = create_random_ugrid(0, 100, 200, 300, 400, 500, npts)

        i_chosen = np.random.randint(0, npts)
        xi = x[i_chosen]
        yi = y[i_chosen]
        zi = z[i_chosen]

        ic, xc, yc, zc = ugrid_closest_point_index(x, y, z, xi, yi, zi)

        self.assertEqual(ic, i_chosen)
        self.assertAlmostEqual(xc, xi)
        self.assertAlmostEqual(yc, yi)
        self.assertAlmostEqual(zc, zi)

if __name__ == '__main__':

    unittest.TextTestRunner(verbosity=2).run(suite())
