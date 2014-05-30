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

        self.plotopt.opdict['otime_window'] = 5.
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


if __name__ == '__main__':

    logging.basicConfig(level=logging.INFO,
                        format='%(levelname)s : %(asctime)s : %(message)s')

    unittest.TextTestRunner(verbosity=2).run(suite())
