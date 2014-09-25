import os
import glob
import h5py
import logging
import unittest
import numpy as np
from options import WavelocOptions
from migration import do_migration_setup_and_run
from synth_migration import generateSyntheticDirac
from plotting import plotWavelocResults


def suite():

    suite = unittest.TestSuite()

    suite.addTest(SyntheticMigrationTests('test_dirac_migration'))
    suite.addTest(MigrationTests('test_take'))
    suite.addTest(MigrationTests('test_migration'))
    suite.addTest(MigrationTests('test_migration_fullRes'))
    suite.addTest(MigrationTests('test_migration_use_ram'))
    suite.addTest(UgridMigrationTests('test_time_grid_ugrids'))
    suite.addTest(UgridMigrationTests('test_syn_ugrid_migration'))

    return suite


def hdf5_to_signature(base_path, datadir, dataglob, output_filename):

    sig_file = open(os.path.join(base_path, datadir, output_filename), 'w')
    allfiles = glob.glob(os.path.join(base_path, datadir, dataglob))
    for filename in allfiles:
        basename = os.path.basename(filename)
        f = h5py.File(filename, 'r')
        for name in f:
            logging.debug('Signature for %s %s : ' % (basename, name))
            dset = f[name]
            maximum = np.max(dset)
            datasum = np.sum(dset)
            datalen = len(dset)
            sig_file.write("%s \t %s \t %.6f \t %.6f \t %d\n" %
                          (basename, name, maximum, datasum, datalen))
        f.close()


def read_signature_file(filename):
    f = open(filename, 'r')
    lines = f.readlines()
    f.close()

    nlines = len(lines)
    maximum = np.empty(nlines, dtype='float')
    datasum = np.empty(nlines, dtype='float')
    datalen = np.empty(nlines, dtype='int')

    for i in xrange(nlines):
        line = lines[i]
        maximum[i] = np.float(line.split()[2])
        datasum[i] = np.float(line.split()[3])
        datalen[i] = np.float(line.split()[4])

    return maximum, datasum, datalen


def signature_comparison(fname1, fname2):

    m1, d1, l1 = read_signature_file(fname1)
    m2, d2, l2 = read_signature_file(fname2)

    r_m = np.sum(np.abs(m1-m2))/np.sum(m1)
    r_d = np.sum(np.abs(d1-d2))/np.sum(d1)
    r_l = np.sum(np.abs(l1-l2))

    return r_m, r_d, r_l


def hdf5_to_sig_values(filename):

    f = h5py.File(filename, 'r')
    sig_values = []
    for name in f:
        logging.debug('Signature for %s : ' % filename)
        dset = f[name]
        maximum = np.max(dset)
        datasum = np.sum(dset)
        datalen = len(dset)
        sig_values.append([maximum, datasum, datalen])
    f.close()

    return sig_values


def hdf5_max_values(filename):

    f = h5py.File(filename, 'r')
    dset = f['max_val']
    np_dset = np.empty(len(dset), dtype='float32')
    np_dset[:] = dset
    f.close()

    return np_dset


class SyntheticMigrationTests(unittest.TestCase):

    def test_dirac_migration(self):
        from locations_trigger import trigger_locations_inner

        outdir = 'TEST_Dirac'

        # run the synthetic test (call external function to avoid copying code)
        wo, plotopt = run_synthetic_test(outdir, ugrid=True)

        # retrieve info
        stack_filename = plotopt.getStackFilename()
        n_buf = plotopt.opdict['n_buf']
        nt = plotopt.opdict['nt']
        dt = plotopt.opdict['dt']
        x_true = plotopt.opdict['x_loc']
        y_true = plotopt.opdict['y_loc']
        z_true = plotopt.opdict['z_loc']
        stack_start_time = plotopt.opdict['start_time']

        # loclevel for triggers
        loclevel = wo.opdict['loclevel']

        # extract the max stacks
        f_stack = h5py.File(stack_filename, 'r')
        max_val = f_stack['max_val']
        max_x = f_stack['max_x']
        max_y = f_stack['max_y']
        max_z = f_stack['max_z']

        locs = trigger_locations_inner(max_val, max_x, max_y, max_z,
                                       loclevel, loclevel,
                                       stack_start_time, dt)
        self.assertTrue(len(locs) > 0)

        # This is a dirac test, so only have one element in locs
        imax = np.argmax([loc['max_trig'] for loc in locs])
        trig_loc = locs[imax]
        self.assertAlmostEqual(wo.opdict['syn_otime'], trig_loc['o_time'], 2)
        self.assertAlmostEqual(wo.opdict['syn_x'], trig_loc['x_mean'])
        self.assertAlmostEqual(wo.opdict['syn_y'], trig_loc['y_mean'])
        self.assertAlmostEqual(wo.opdict['syn_z'], trig_loc['z_mean'])

        f_stack.close()

        # run the plotting
        plotWavelocResults(plotopt)


#@unittest.skip('Skip for rapidity')
class MigrationTests(unittest.TestCase):

    def setUp(self):

        self.wo = WavelocOptions()
        self.wo.set_test_options()
        self.wo.verify_migration_options()

    def test_migration(self):

        self.wo.opdict['load_ttimes_buf'] = False
        self.wo.opdict['data_length'] = 300

        base_path = self.wo.opdict['base_path']
        test_datadir = self.wo.opdict['test_datadir']
        outdir = self.wo.opdict['outdir']

        # do migration
        do_migration_setup_and_run(self.wo.opdict)

        # write signature
        hdf5_to_signature(base_path, os.path.join('out', outdir, 'stack'),
                          'stack_all_2010-10-14T00:14:00.000000Z.hdf5',
                          'stack_signature.dat')
        signature_filename = os.path.join(base_path, 'out', outdir, 'stack',
                                          'stack_signature.dat')

        # read expected signature
        expected_signature_filename =\
            os.path.join(base_path, test_datadir, 'test_stack_signature.dat')

        m, d, l = signature_comparison(signature_filename,
                                       expected_signature_filename)
        self.assertAlmostEqual(m, 0.0)
        self.assertAlmostEqual(d, 0.0)
        self.assertAlmostEqual(l, 0.0)

    def test_take(self):

        nbuf = 9000
        nt = 3600

        matrix = np.random.randn(nbuf, nt)

        max_ib = np.argmax(matrix, 0)
        max_val = np.max(matrix, 0)
        max_val_take = np.diag(matrix.take(max_ib, 0))

        self.assertEqual(len(max_ib), nt)
        self.assertEqual(len(max_val), nt)
        self.assertEqual(len(max_val_take), nt)
        self.assertEqual(max_val.shape, max_val_take.shape)
        np.testing.assert_allclose(max_val_take, max_val)

    def test_migration_use_ram(self):

        self.wo.opdict['load_ttimes_buf'] = True
        self.wo.opdict['data_length'] = 300

        base_path = self.wo.opdict['base_path']
        outdir = self.wo.opdict['outdir']

        # do migration without use_ram
        self.wo.opdict['use_ram'] = False
        do_migration_setup_and_run(self.wo.opdict)
        lines_no_ram =\
            hdf5_max_values(
                os.path.join(base_path, 'out', outdir, 'stack',
                             'stack_all_2010-10-14T00:14:00.000000Z.hdf5'))

        # do migration with use_ram
        self.wo.opdict['use_ram'] = True
        do_migration_setup_and_run(self.wo.opdict)
        lines_use_ram =\
            hdf5_max_values(
                os.path.join(base_path, 'out', outdir, 'stack',
                             'stack_all_2010-10-14T00:14:00.000000Z.hdf5'))

        # verify that the two give the same result
        np.testing.assert_allclose(lines_use_ram, lines_no_ram)

    @unittest.skip('Time-consuming and uninteresting for code checking')
    def test_migration_fullRes(self):

        self.wo.opdict['search_grid'] = 'grid.Taisne.search.hdr'
        self.wo.opdict['outdir'] = 'TEST_fullRes'
        self.wo.opdict['load_ttimes_buf'] = False
        self.wo.opdict['data_length'] = 300
        self.wo.opdict['use_ram'] = True
        self.wo.verify_migration_options()

        base_path = self.wo.opdict['base_path']
        test_datadir = self.wo.opdict['test_datadir']
        outdir = self.wo.opdict['outdir']

        expected_signature_filename =\
            os.path.join(base_path, test_datadir,
                         'TEST_fullRes_stack_signature.dat')

        do_migration_setup_and_run(self.wo.opdict)

        hdf5_to_signature(base_path, os.path.join('out', outdir, 'stack'),
                          'stack*hdf5', 'stack_signature.dat')
        signature_filename = os.path.join(base_path, 'out', outdir, 'stack',
                                          'stack_signature.dat')

        m, d, l = signature_comparison(signature_filename,
                                       expected_signature_filename)
        self.assertAlmostEqual(m, 0.0)
        self.assertAlmostEqual(d, 0.0)
        self.assertAlmostEqual(l, 0.0)

#@unittest.skip('Skip for rapidity')
class UgridMigrationTests(unittest.TestCase):

    def test_time_grid_ugrids(self):

        from hdf5_grids import get_interpolated_time_grids,\
            get_interpolated_time_ugrids

        wo = WavelocOptions()
        wo.set_test_options()
        wo.verify_base_path()

        # force creation of time grids
        wo.opdict['load_ttimes_buf'] = False

        # do old-style interpolation
        wo.opdict['outdir'] = 'TEST_Dirac'
        wo.verify_migration_options()
        x, y, z, time_grids = get_interpolated_time_grids(wo.opdict)

        # do new-style interpolation
        wo.opdict['outdir'] = 'TEST_Dirac_ugrid'
        wo.opdict['ugrid_type'] = 'FULL'
        wo.verify_migration_options()
        ux, uy, yz, time_ugrids = get_interpolated_time_ugrids(wo.opdict)

        # check that we have the same time grids
        for sta in time_grids.keys():
            np.testing.assert_allclose(time_grids[sta], time_ugrids[sta])

    def test_syn_ugrid_migration(self):

        # run the same synthetic test the old way and with ugrid
        wo, test_info = run_synthetic_test('TEST_Dirac', ugrid=False)
        wo_ugrid, test_info_ugrid = run_synthetic_test('TEST_Dirac_ugrid',
                                                       ugrid=True)

        # get old-style max_val
        stack_filename = test_info.getStackFilename()
        f = h5py.File(stack_filename, 'r')
        mv = f['max_val']
        max_val = mv[:]
        f.close()

        # get new-style max_val
        stack_filename = test_info_ugrid.getStackFilename()
        f = h5py.File(stack_filename, 'r')
        mv = f['max_val']
        max_val_ugrid = mv[:]
        f.close()

        # compare
        np.testing.assert_allclose(max_val, max_val_ugrid)


def run_synthetic_test(outdir, ugrid=False):

    wo = WavelocOptions()
    wo.set_test_options()

    wo.opdict['outdir'] = outdir
    wo.opdict['search_grid'] = 'grid.Taisne.search.hdr'
    wo.opdict['loclevel'] = 10
    wo.opdict['load_ttimes_buf'] = False
    wo.opdict['syn_addnoise'] = False
    wo.opdict['syn_amplitude'] = 1.0
    wo.opdict['syn_datalength'] = 20.0
    wo.opdict['syn_samplefreq'] = 100.0
    wo.opdict['syn_kwidth'] = 0.1
    wo.opdict['syn_otime'] = 6.0
    wo.opdict['syn_x'] = 366.
    wo.opdict['syn_y'] = 7650.
    wo.opdict['syn_z'] = -1.
    wo.opdict['syn_filename'] = 'test_grid4D_hires.hdf5'
    wo.opdict['ugrid_type'] = 'FULL'
    wo.opdict['otime_window'] = 5.0

    wo.verify_base_path
    wo.verify_synthetic_options()

    ##########################
    # generate the test case and retrieve necessary information
    ##########################

    if ugrid:
        plotopt = generateSyntheticDirac(wo.opdict, ugrid=True)
    else:
        plotopt = generateSyntheticDirac(wo.opdict, ugrid=False)
    return wo, plotopt


if __name__ == '__main__':

    logging.basicConfig(level=logging.INFO,
                        format='%(levelname)s : %(asctime)s : %(message)s')

    unittest.TextTestRunner(verbosity=2).run(suite())
