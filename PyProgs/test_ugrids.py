import os
import unittest
import numpy as np
from NllGridLib import read_hdr_file
from ugrids import create_random_ugrid, nll2random_ugrid, nll2reg_ugrid,\
    ugrid_closest_point_index, ugrid_svr, select_points_closeto_plane


def suite():

    suite = unittest.TestSuite()

    suite.addTest(UgridTests('test_create_read_random'))
    suite.addTest(UgridTests('test_nll2random_ugrid'))
    suite.addTest(UgridTests('test_nll2reg_ugrid'))
    suite.addTest(UgridTests('test_closest_point'))
    suite.addTest(UgridTests('test_ugrid_svr'))
    suite.addTest(UgridTests('test_select_close_points'))

    return suite


class UgridTests(unittest.TestCase):

    def setUp(self):

        self.npts = 1000
        self.xmin, self.ymin, self.zmin = np.random.uniform(low=0., high=10.,
                                                            size=3)
        self.xmax, self.ymax, self.zmax = np.random.uniform(low=20., high=50.,
                                                            size=3)
        self.x, self.y, self.z = create_random_ugrid(self.xmin, self.xmax,
                                                     self.ymin, self.ymax,
                                                     self.zmin, self.zmax,
                                                     self.npts)

    def test_create_read_random(self):

        self.assertEqual(len(self.x), self.npts)
        self.assertTrue(np.min(self.x) >= self.xmin)
        self.assertTrue(np.max(self.x) < self.xmax)
        self.assertTrue(np.min(self.y) >= self.ymin)
        self.assertTrue(np.max(self.y) < self.ymax)
        self.assertTrue(np.min(self.z) >= self.zmin)
        self.assertTrue(np.max(self.z) < self.zmax)

    def test_nll2random_ugrid(self):

        npts = self.npts

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

        i_chosen = np.random.randint(0, self.npts)
        xi = self.x[i_chosen]
        yi = self.y[i_chosen]
        zi = self.z[i_chosen]

        ic, xc, yc, zc = ugrid_closest_point_index(self.x, self.y, self.z,
                                                   xi, yi, zi)

        self.assertEqual(ic, i_chosen)
        self.assertAlmostEqual(xc, xi)
        self.assertAlmostEqual(yc, yi)
        self.assertAlmostEqual(zc, zi)

    def test_select_close_points(self):

        x_coord = self.xmin+(self.xmax-self.xmin)/2.
        x_dist = (self.xmax-self.xmin)/10.
        iclose_x = select_points_closeto_plane(self.x, self.y, self.z, 'x',
                                               x_coord, x_dist)
        y_coord = self.ymin+(self.ymax-self.ymin)/2.
        y_dist = (self.ymax-self.ymin)/10.
        iclose_y = select_points_closeto_plane(self.x, self.y, self.z, 'y',
                                               y_coord, y_dist)
        z_coord = self.zmin+(self.zmax-self.zmin)/2.
        z_dist = (self.zmax-self.zmin)/10.
        iclose_z = select_points_closeto_plane(self.x, self.y, self.z, 'z',
                                               z_coord, z_dist)

        x_close = self.x[iclose_x]
        self.assertLess(len(x_close), self.npts)
        self.assertLess(np.max(x_close), x_coord+x_dist)
        self.assertGreater(np.min(x_close), x_coord-x_dist)

        y_close = self.y[iclose_y]
        self.assertLess(len(y_close), self.npts)
        self.assertLess(np.max(y_close), y_coord+y_dist)
        self.assertGreater(np.min(y_close), y_coord-y_dist)

        z_close = self.z[iclose_z]
        self.assertLess(len(z_close), self.npts)
        self.assertLess(np.max(z_close), z_coord+z_dist)
        self.assertGreater(np.min(z_close), z_coord-z_dist)

    @unittest.expectedFailure
    def test_ugrid_svr(self):

        npts = 1000
        x0 = 50
        y0 = 250
        z0 = 450
        sigma = 10
        x, y, z = create_random_ugrid(0, 100, 200, 300, 400, 500, npts)
        values = np.exp(-((x-x0)**2/(2*sigma**2) +
                          (y-y0)**2/(2*sigma**2) +
                          (z-z0)**2/(2*sigma**2)))

        values_i = ugrid_svr(x, y, z, values, x, y, z)

        np.testing.assert_allclose(values, values_i, rtol=1e-3)


if __name__ == '__main__':

    unittest.TextTestRunner(verbosity=2).run(suite())
