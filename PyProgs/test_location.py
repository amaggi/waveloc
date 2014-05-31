import os
import unittest
import h5py
import numpy as np
from options import WavelocOptions
from locations_trigger import do_locations_trigger_setup_and_run, \
    trigger_locations_inner, read_locs_from_file
from locations_prob import do_locations_prob_setup_and_run, \
    read_prob_locs_from_file
from integrate4D import compute_integral4D, compute_expected_coordinates4D


def suite():
    suite = unittest.TestSuite()
    suite.addTest(IntegrationTests('test_integration'))
    suite.addTest(IntegrationTests('test_expected_values'))
    suite.addTest(IntegrationTests('test_reshaping'))
    suite.addTest(LocationTests('test_locations_trigger'))
    suite.addTest(LocationTests('test_locations_prob'))
    suite.addTest(TriggeringTests('test_simple_trigger'))
    suite.addTest(TriggeringTests('test_smoothing'))
    suite.addTest(TriggeringTests('test_gaussian_trigger'))
    return suite


class TriggeringTests(unittest.TestCase):

    def test_simple_trigger(self):

        max_val = np.random.rand(100)
        max_x = np.random.rand(100)
        max_y = np.random.rand(100)
        max_z = np.random.rand(100)

        max_val[10] = 5
        max_val[20] = 7
        max_val[45] = 2
        max_val[80] = 10

        left_trig = right_trig = 3
        locs = trigger_locations_inner(max_val, max_x, max_y, max_z, left_trig,
                                       right_trig, 0.0, 1.0)
        self.assertEqual(len(locs), 3)
        self.assertAlmostEqual(locs[0]['max_trig'], 5)
        self.assertAlmostEqual(locs[1]['max_trig'], 7)
        self.assertAlmostEqual(locs[2]['max_trig'], 10)
        self.assertAlmostEqual(locs[0]['o_time'], 10)
        self.assertAlmostEqual(locs[1]['o_time'], 20)
        self.assertAlmostEqual(locs[2]['o_time'], 80)

    def test_smoothing(self):

        from filters import smooth

        x = np.arange(100)
        max_val = 100.*np.exp(-(x-50.)*(x-50.)/(10.*10.))+np.random.rand(100)
        max_x = np.random.rand(100)
        max_y = np.random.rand(100)
        max_z = np.random.rand(100)

        max_val_smooth = smooth(max_val)

        left_trig = right_trig = 3
        locs_smooth = trigger_locations_inner(max_val_smooth, max_x, max_y,
                                              max_z, left_trig, right_trig,
                                              0.0, 1.0)
        self.assertAlmostEqual(locs_smooth[0]['o_time'], 50., 2)

    def test_gaussian_trigger(self):

        x = np.arange(100)
        max_val = 10.*np.exp(-(x-50.)*(x-50.)/(10.*10.))
        max_x = np.random.rand(100)
        max_y = np.random.rand(100)
        max_z = np.random.rand(100)

        left_trig = right_trig = 3
        locs = trigger_locations_inner(max_val, max_x, max_y, max_z, left_trig,
                                       right_trig, 0.0, 1.0)
        self.assertAlmostEqual(locs[0]['max_trig'], 10)
        self.assertAlmostEqual(locs[0]['o_time'], 50)


class IntegrationTests(unittest.TestCase):

    def test_integration(self):

        dims = (40, 60, 80, 100)
        grid4D = np.ones(dims)
        x0 = np.linspace(0, 4, dims[0])
        x1 = np.linspace(0, 6, dims[1])
        x2 = np.linspace(0, 8, dims[2])
        x3 = np.linspace(0, 10, dims[3])
        grid_area = 4.0*6*8*10

        grid_integral = compute_integral4D(grid4D, x0, x1, x2, x3)

        self.assertAlmostEqual(grid_area, grid_integral, 7)

    def test_expected_values(self):

        dims = (20, 30, 40, 50)
        grid4D = np.zeros(dims)
        x0 = np.linspace(0, 4, dims[0])
        x1 = np.linspace(0, 6, dims[1])
        x2 = np.linspace(0, 8, dims[2])
        x3 = np.linspace(0, 10, dims[3])
        grid4D[1, 2, 3, 7] = 2
        grid4D[1, 2, 4, 7] = 4   # add something interesting to find
        grid4D[1, 2, 5, 7] = 2
        my_exp0 = x0[1]
        my_exp1 = x1[2]
        my_exp2 = x2[4]
        my_exp3 = x3[7]

        exp0, exp1, exp2, exp3, cov_matrix = \
            compute_expected_coordinates4D(grid4D, x0, x1, x2, x3)
        var_x0 = cov_matrix[0, 0]
        var_x1 = cov_matrix[1, 1]
        var_x2 = cov_matrix[2, 2]
        var_x3 = cov_matrix[3, 3]

        self.assertAlmostEqual(my_exp0, exp0, 7)
        self.assertAlmostEqual(my_exp1, exp1, 7)
        self.assertAlmostEqual(my_exp2, exp2, 7)
        self.assertAlmostEqual(my_exp3, exp3, 7)

        self.assertAlmostEqual(var_x0, 0.0, 7)
        self.assertAlmostEqual(var_x1, 0.0, 7)
        self.assertAlmostEqual(var_x2, 0.0210, 4)
        self.assertAlmostEqual(var_x3, 0.0, 7)

    def test_reshaping(self):
        true_dims = (20, 30, 40, 50)
        dims_2D = (20*30*40, 50)
        array_2D = np.zeros(dims_2D)
        array_2D[:, 20] = 1

        array_4D = array_2D.reshape(true_dims)
        np.testing.assert_allclose(array_4D.shape,  true_dims)
        np.testing.assert_allclose(array_4D[:, :, :, 20],
                                   np.ones((20, 30, 40)))
        np.testing.assert_allclose(array_4D[:, :, :, 21],
                                   np.zeros((20, 30, 40)))


class LocationTests(unittest.TestCase):

    def setUp(self):

        self.wo = WavelocOptions()
        self.wo.set_test_options()
        self.wo.verify_location_options()

    def test_locations_trigger(self):

        self.wo.opdict['outdir'] = 'TEST'
        self.wo.verify_location_options()

        base_path = self.wo.opdict['base_path']
        test_datadir = self.wo.opdict['test_datadir']
        outdir = self.wo.opdict['outdir']

        exp_loc_fname = os.path.join(base_path, test_datadir,
                                     'TEST_locations.dat')
        exp_locs = read_locs_from_file(exp_loc_fname)

        locs = do_locations_trigger_setup_and_run(self.wo.opdict)

        loc_fname = os.path.join(base_path, 'out', outdir, 'loc',
                                 'locations.dat')
        locs = read_locs_from_file(loc_fname)

        self.assertEqual(len(locs), len(exp_locs))
        for i in xrange(len(locs)):
            loc = locs[i]
            exp_loc = exp_locs[i]
            self.assertGreater(loc['o_time'],
                               exp_loc['o_time']-exp_loc['o_err_left'])
            self.assertLess(loc['o_time'],
                            exp_loc['o_time']+exp_loc['o_err_right'])
            self.assertLess(np.abs(loc['x_mean']-exp_loc['x_mean']),
                            exp_loc['x_sigma'])
            self.assertLess(np.abs(loc['x_mean']-exp_loc['x_mean']),
                            loc['x_sigma'])
            self.assertLess(np.abs(loc['y_mean']-exp_loc['y_mean']),
                            exp_loc['y_sigma'])
            self.assertLess(np.abs(loc['y_mean']-exp_loc['y_mean']),
                            loc['y_sigma'])
            self.assertLess(np.abs(loc['z_mean']-exp_loc['z_mean']),
                            exp_loc['z_sigma'])
            self.assertLess(np.abs(loc['z_mean']-exp_loc['z_mean']),
                            loc['z_sigma'])

    @unittest.skip('Skipping prob loc until location stuff is ok')
    def test_locations_prob(self):

        self.wo.opdict['outdir'] = 'TEST'
        self.wo.opdict['probloc_spaceonly'] = True
        self.wo.verify_location_options()

        base_path = self.wo.opdict['base_path']
        outdir = self.wo.opdict['outdir']

        loc_fname = os.path.join(base_path, 'out', outdir, 'loc',
                                 'locations.dat')
        prob_fname = os.path.join(base_path, 'out', outdir, 'loc',
                                  'locations_prob.dat')
        hdf5_fname = os.path.join(base_path, 'out', outdir, 'loc',
                                  'locations_prob.hdf5')

        do_locations_prob_setup_and_run(self.wo.opdict)

        locs = read_locs_from_file(loc_fname)
        prob_locs = read_prob_locs_from_file(prob_fname)
        f_marginals = h5py.File(hdf5_fname, 'r')
        self.assertEqual(len(locs), len(prob_locs))

        for i in xrange(len(locs)):
            loc = locs[i]
            prob_loc = prob_locs[i]
            self.assertGreater(prob_loc['o_time'],
                               loc['o_time']-loc['o_err_left'])
            self.assertLess(prob_loc['o_time'],
                            loc['o_time']+loc['o_err_right'])
            self.assertLess(np.abs(loc['o_time']-prob_loc['o_time']),
                            prob_loc['o_err'])
            self.assertLess(np.abs(loc['x_mean']-prob_loc['x_mean']),
                            prob_loc['x_sigma'])
            self.assertLess(np.abs(loc['y_mean']-prob_loc['y_mean']),
                            prob_loc['y_sigma'])
            self.assertLess(np.abs(loc['z_mean']-prob_loc['z_mean']),
                            prob_loc['z_sigma'])

            grp = f_marginals[prob_loc['o_time'].isoformat()]
            nx = grp['x'].shape[0]
            ny = grp['y'].shape[0]
            nz = grp['z'].shape[0]
            self.assertEqual(grp['prob_x'].shape, (nx, ))
            self.assertEqual(grp['prob_y'].shape, (ny, ))
            self.assertEqual(grp['prob_z'].shape, (nz, ))
            self.assertEqual(grp['prob_xy'].shape, (nx, ny))
            self.assertEqual(grp['prob_xz'].shape, (nx, nz))
            self.assertEqual(grp['prob_yz'].shape, (ny, nz))
            # if is a 4D grid
            if 't' in grp:
                nt = grp['t'].shape[0]
                self.assertEqual(grp['prob_t'].shape, (nt, ))
                self.assertEqual(grp['prob_xt'].shape, (nx, nt))
                self.assertEqual(grp['prob_yt'].shape, (ny, nt))
                self.assertEqual(grp['prob_zt'].shape, (nz, nt))

        f_marginals.close()


if __name__ == '__main__':

    import logging
    logging.basicConfig(level=logging.INFO,
                        format='%(levelname)s : %(asctime)s : %(message)s')

    unittest.TextTestRunner(verbosity=2).run(suite())
