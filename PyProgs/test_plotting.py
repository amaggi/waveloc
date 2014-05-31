import os
import h5py
import unittest
import logging
import numpy as np
from options import WavelocOptions
from plot_options import PlotOptions
from plotting import plotWavelocResults


def suite():

    suite = unittest.TestSuite()
    suite.addTest(PlottingTests('test_plotOptions'))
    suite.addTest(PlottingTests('test_plotWavelocResults'))
    suite.addTest(FullPlottingTests('test_waveloc_example'))

    return suite


class PlottingTests(unittest.TestCase):

    def setUp(self):

        wo = WavelocOptions()
        wo.set_test_options()
        wo.verify_base_path()

        self.plotopt = PlotOptions(wo.opdict)
        self.plotopt.opdict['grid_filename'] = 'GRID_FNAME.hdf5'
        self.plotopt.opdict['stack_filename'] = 'STACK_FNAME.hdf5'

        self._create_dummy_grid()
        self._create_dummy_stack()

        self.plotopt.opdict['plot_otime_window'] = 5.
        self.plotopt.opdict['t_err'] = (0.5, 0.5)
        self.plotopt.opdict['x_err'] = (0.5, 0.5)
        self.plotopt.opdict['y_err'] = (0.5, 0.5)
        self.plotopt.opdict['z_err'] = (0.1, 0.1)

    def _create_dummy_grid(self):
        
        x, y, z = self.plotopt.getXYZ()
        dt = 0.5
        tlen = 50.

        n_buf = len(x)
        nt = int(tlen/dt)
        grid = np.empty((n_buf, nt), dtype='float')
        t = np.arange(0, tlen, dt)

        x_range = np.max(x)-np.min(x)
        y_range = np.max(y)-np.min(y)
        z_range = np.max(z)-np.min(z)

        xc = np.min(x)+x_range/2.
        yc = np.min(y)+y_range/3.
        zc = np.min(z)+z_range/5.
        tc = tlen/2.0

        a = 0.75*x_range/2
        b = 0.75*y_range/2
        c = 0.75*z_range/2

        for it in xrange(nt):
            for ib in xrange(n_buf):
                val = ((x[ib]-xc)/a)**2 + ((y[ib]-yc)/b)**2 + ((z[ib]-zc)/c)**2
                if val > 1:
                    grid[ib, it] = 0.
                else:
                    grid[ib, it] = 1-val

        # write file
        grid_filename = self.plotopt.getGridFilename()
        f = h5py.File(grid_filename, 'w')
        mg = f.create_dataset('migrated_grid', data=grid)
        mg.attrs['n_buf'] = n_buf
        mg.attrs['nt'] = nt
        mg.attrs['dt'] = dt
        f.create_dataset('x', data=x)
        f.create_dataset('y', data=y)
        f.create_dataset('z', data=z)
        f.close()

        # save stuff that needs to be saved in plot options
        self.plotopt.opdict['dt'] = dt
        self.plotopt.opdict['nt'] = nt
        self.plotopt.opdict['x_loc'] = xc
        self.plotopt.opdict['y_loc'] = yc
        self.plotopt.opdict['z_loc'] = zc
        self.plotopt.opdict['t_loc_rel'] = tc

    def _create_dummy_stack(self):

        x, y, z = self.plotopt.getXYZ()
        dt = 0.5
        tlen = 50
        nt = int(tlen/dt)
        t = np.arange(0, tlen, dt)

        x_range = np.max(x)-np.min(x)
        y_range = np.max(y)-np.min(y)
        z_range = np.max(z)-np.min(z)

        xc = np.min(x)+x_range/2.
        yc = np.min(y)+y_range/3.
        zc = np.min(z)+z_range/5.
        tc = tlen/2.0

        max_val = np.empty(nt, dtype='float')
        max_x = np.empty(nt, dtype='float')
        max_y = np.empty(nt, dtype='float')
        max_z = np.empty(nt, dtype='float')

        sig=tlen/30.
        max_val = np.exp(-(t-tc)**2/(2*sig**2))
        max_x = np.random.normal(loc=xc, scale=0.5, size=nt)
        max_y = np.random.normal(loc=yc, scale=0.5, size=nt)
        max_z = np.random.normal(loc=zc, scale=0.1, size=nt)

        stack_filename = self.plotopt.getStackFilename()
        f = h5py.File(stack_filename, 'w')
        f.create_dataset('max_val', data=max_val)
        mv = f.create_dataset('max_val_smooth', data=max_val)
        f.create_dataset('max_x', data=max_x)
        f.create_dataset('max_y', data=max_y)
        f.create_dataset('max_z', data=max_z)
        f.create_dataset('t', data=t)
        f.close()

    def test_plotWavelocResults(self):

        plotWavelocResults(self.plotopt)

    def test_plotOptions(self):

        # check that the grid filenames are sensible
        exp_grid_filename = os.path.join(self.plotopt.opdict['base_path'],
                                         'out', 'TEST', 'grid',
                                         'GRID_FNAME.hdf5')
        exp_stack_filename = os.path.join(self.plotopt.opdict['base_path'],
                                         'out', 'TEST', 'stack',
                                         'STACK_FNAME.hdf5')
        exp_fig_filename = os.path.join(self.plotopt.opdict['base_path'],
                                         'out', 'TEST', 'fig',
                                         'GRID_FNAME_grid.pdf')
        exp_wfm_fig_filename = os.path.join(self.plotopt.opdict['base_path'],
                                         'out', 'TEST', 'fig',
                                         'GRID_FNAME_wfm.pdf')

        grid_filename = self.plotopt.getGridFilename()
        stack_filename = self.plotopt.getStackFilename()
        fig_filename = self.plotopt.getFigFilename()
        wfm_fig_filename = self.plotopt.getWfmFigFilename()

        self.assertEqual(grid_filename, exp_grid_filename)
        self.assertEqual(stack_filename, exp_stack_filename)
        self.assertEqual(fig_filename, exp_fig_filename)
        self.assertEqual(wfm_fig_filename, exp_wfm_fig_filename)

class FullPlottingTests(unittest.TestCase):

    def test_waveloc_example(self):
        from options import WavelocOptions
        from SDS_processing import do_SDS_processing_setup_and_run
        from migration import do_migration_setup_and_run
        from locations_trigger import do_locations_trigger_setup_and_run
        from plot_locations2 import do_plotting_setup_and_run

        # set up default parameters
        wo = WavelocOptions()
        wo.verify_base_path()

        wo.opdict['time'] = True
        wo.opdict['verbose'] = False

        wo.opdict['test_datadir'] = 'test_data'
        wo.opdict['datadir'] = 'TEST'
        wo.opdict['outdir'] = 'TEST_fullRes'

        wo.opdict['net_list'] = 'YA'
        wo.opdict['sta_list'] = "FJS,FLR,FOR,HDL,RVL,SNE,UV01,UV02,UV03,UV04,UV05,\
                                 UV06,UV07,UV08,UV09,UV10,UV11,UV12,UV13,UV14,UV15"
        wo.opdict['comp_list'] = "HHZ"

        wo.opdict['starttime'] = "2010-10-14T00:14:00.0Z"
        wo.opdict['endtime'] = "2010-10-14T00:18:00.0Z"

        wo.opdict['time_grid'] = 'Slow_len.100m.P'
        wo.opdict['search_grid'] = 'grid.Taisne.search.hdr'
        wo.opdict['stations'] = 'coord_stations_test'

        wo.opdict['resample'] = False
        wo.opdict['fs'] = None

        wo.opdict['c1'] = 4.0
        wo.opdict['c2'] = 10.0

        wo.opdict['kwin'] = 4
        wo.opdict['krec'] = False
        wo.opdict['kderiv'] = True

        wo.opdict['data_length'] = 300
        wo.opdict['data_overlap'] = 20

        wo.opdict['dataglob'] = '*filt.mseed'
        wo.opdict['kurtglob'] = '*kurt.mseed'
        wo.opdict['gradglob'] = '*grad.mseed'

        wo.opdict['load_ttimes_buf'] = True

        wo.opdict['loclevel'] = 5000.0
        wo.opdict['snr_limit'] = 10.0
        wo.opdict['sn_time'] = 10.0
        wo.opdict['n_kurt_min'] = 4

        wo.opdict['plot_tbefore'] = 4
        wo.opdict['plot_tafter'] = 6
        wo.opdict['plot_otime_window'] = 2

        ##########################################
        # end of option setting - start processing
        ##########################################

        wo.verify_SDS_processing_options()
        #do_SDS_processing_setup_and_run(wo.opdict)

        wo.verify_migration_options()
        #do_migration_setup_and_run(wo.opdict)

        # do trigger location
        wo.verify_location_options()
        #do_locations_trigger_setup_and_run(wo.opdict)

        # This will do plotting of grids and stacks for locations
        wo.verify_plotting_options()
        do_plotting_setup_and_run(wo.opdict, plot_wfm=True, plot_grid=True)

if __name__ == '__main__':

    logging.basicConfig(level=logging.INFO,
                        format='%(levelname)s : %(asctime)s : %(message)s')

    unittest.TextTestRunner(verbosity=2).run(suite())
