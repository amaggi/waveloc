import os
import glob
import h5py
import logging
import unittest
import numpy as np
from options import WavelocOptions
from migration import do_migration_setup_and_run
from synth_migration import generateSyntheticDirac


def suite():

    suite = unittest.TestSuite()

    suite.addTest(SyntheticMigrationTests('test_dirac_migration'))
    suite.addTest(MigrationTests('test_take'))
    suite.addTest(MigrationTests('test_migration'))
    suite.addTest(MigrationTests('test_migration_fullRes'))
    suite.addTest(MigrationTests('test_migration_use_ram'))

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

        wo = WavelocOptions()
        wo.set_test_options()

        wo.opdict['outdir'] = 'TEST_Dirac'
        wo.opdict['search_grid'] = 'grid.Taisne.search.hdr'
        wo.opdict['loclevel'] = 10
        wo.opdict['load_ttimes_buf'] = True
        wo.opdict['syn_addnoise'] = False
        wo.opdict['syn_amplitude'] = 1.0
        wo.opdict['syn_datalength'] = 20.0
        wo.opdict['syn_samplefreq'] = 100.0
        wo.opdict['syn_kwidth'] = 0.1
        wo.opdict['syn_otime'] = 6.0
        wo.opdict['syn_ix'] = 16
        wo.opdict['syn_iy'] = 8
        wo.opdict['syn_iz'] = 6
        wo.opdict['syn_filename'] = 'test_grid4D_hires.hdf5'

        wo.verify_synthetic_options()

        ##########################
        # generate the test case and retrieve necessary information
        ##########################

        logging.info('Running synthetic test case generation...')
        test_info = generateSyntheticDirac(wo.opdict)
        logging.debug(test_info)

        # retrieve info
        stack_filename = test_info['stack_file']
        nx, ny, nz, nt = test_info['grid_shape']
        dx, dy, dz, dt = test_info['grid_spacing']
        x_orig, y_orig, z_orig = test_info['grid_orig']
        ix_true, iy_true, iz_true, it_true = test_info['true_indexes']
        stack_start_time = test_info['start_time']

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
        self.assertAlmostEqual(wo.opdict['syn_ix']*dx+x_orig,
                               trig_loc['x_mean'])
        self.assertAlmostEqual(wo.opdict['syn_iy']*dy+y_orig,
                               trig_loc['y_mean'])
        self.assertAlmostEqual(wo.opdict['syn_iz']*dz+z_orig,
                               trig_loc['z_mean'])

        f_stack.close()


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

        expected_signature_filename = os.path.join(base_path,
                                                   test_datadir,
                                                   'test_stack_signature.dat')
        expected_signature_file = open(expected_signature_filename, 'r')
        expected_lines = expected_signature_file.readlines()

        do_migration_setup_and_run(self.wo.opdict)

        hdf5_to_signature(base_path, os.path.join('out', outdir, 'stack'),
                          'stack_all_2010-10-14T00:14:00.000000Z.hdf5',
                          'stack_signature.dat')
        signature_filename = os.path.join(base_path, 'out', outdir, 'stack',
                                          'stack_signature.dat')
        signature_file = open(signature_filename, 'r')
        lines = signature_file.readlines()

        self.assertSequenceEqual(lines, expected_lines)

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
        expected_signature_file = open(expected_signature_filename, 'r')
        expected_lines = expected_signature_file.readlines()

        do_migration_setup_and_run(self.wo.opdict)

        hdf5_to_signature(base_path, os.path.join('out', outdir, 'stack'),
                          'stack*hdf5', 'stack_signature.dat')
        signature_filename = os.path.join(base_path, 'out', outdir, 'stack',
                                          'stack_signature.dat')
        signature_file = open(signature_filename, 'r')
        lines = signature_file.readlines()

        self.assertSequenceEqual(lines, expected_lines)


if __name__ == '__main__':

    logging.basicConfig(level=logging.INFO,
                        format='%(levelname)s : %(asctime)s : %(message)s')

    unittest.TextTestRunner(verbosity=2).run(suite())
