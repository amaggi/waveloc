import os
import unittest
import numpy as np
from waveloc.kurtogram import Fast_Kurtogram, getBandwidthAndFrequency, \
    get_h_parameters, Find_wav_kurt, getFTSquaredEnvelope


def suite():
    suite = unittest.TestSuite()
    suite.addTest(KurtogramTests('test_MatlabGrid'))
    suite.addTest(KurtogramTests('test_BandwidthFrequency'))
    suite.addTest(KurtogramTests('test_FindWavKurt'))
    return suite


class KurtogramTests(unittest.TestCase):

    def setUp(self):
        # read the input data for the test - this is the original data from
        # Antoni
        # input data come in the form of a Matlab file
        from scipy.io.matlab import loadmat

        base_path = os.getenv('WAVELOC_PATH')
        matlab_file = os.path.join(base_path, 'test_data', 'VOIE1.mat')
        v1 = loadmat(matlab_file)
        self.x = v1['v1']
        self.Fs = 100
        self.nlevel = 8
        # calculate the kurtogram
        self.grid, self.level_w, self.freq_w, self.c, self.flow, self.fup = \
            Fast_Kurtogram(self.x, self.nlevel, Fs=self.Fs)

    def test_MatlabGrid(self):

        # Test you get the same output for the kurtogram grid as the original
        # Matlab code does
        base_path = os.getenv('WAVELOC_PATH')
        np_file = os.path.join(base_path, 'test_data', 'matlab_grid.np')
        matlab_grid = np.fromfile(np_file).reshape(16, 768)

        # do test
        np.testing.assert_allclose(self.grid, matlab_grid, atol=0.003)

    def test_BandwidthFrequency(self):

        level_index = 11
        freq_index = 24
        level_exp = 6
        bw_exp = 0.78125
        fc_exp = 1.953125

        # get the bandwith and central frequency
        bw, fc, fi, l1 = \
            getBandwidthAndFrequency(self.nlevel, self.Fs, self.level_w,
                                     self.freq_w, level_index, freq_index)
        # do test
        self.assertAlmostEqual(bw, bw_exp)
        self.assertAlmostEqual(fc, fc_exp)
        self.assertEqual(self.level_w[level_index], level_exp)

    def test_FindWavKurt(self):
        from scipy.io.matlab import loadmat

        N = 16
        fcut = 0.4
        level_index = 11
        freq_index = 24
        lev = self.level_w[level_index]

        base_path = os.getenv('WAVELOC_PATH')
        matlab_file = os.path.join(base_path, 'test_data', 'c.mat')
        c_dict = loadmat(matlab_file)
        c_exp = c_dict['c']

        matlab_file = os.path.join(base_path, 'test_data', 'S.mat')
        S_dict = loadmat(matlab_file)
        S_exp = S_dict['S']

        # get bw and frequency (Hz)
        bw_hz, fc_hz, fi, l1 = \
            getBandwidthAndFrequency(self.nlevel, self.Fs, self.level_w,
                                     self.freq_w, level_index, freq_index)

        # get basic filter parameters
        h, g, h1, h2, h3 = get_h_parameters(N, fcut)
        c, s, threshold, Bw, fc = Find_wav_kurt(self.x, h, g, h1, h2, h3,
                                                self.nlevel, lev, fi,
                                                Fs=self.Fs)
        S = getFTSquaredEnvelope(c)

        # do tests
        self.assertAlmostEqual(Bw*self.Fs, bw_hz)
        self.assertAlmostEqual(fc*self.Fs, fc_hz)
        np.testing.assert_allclose(c.flatten(), c_exp.flatten(), atol=1e-3)
        np.testing.assert_allclose(S.flatten(), S_exp.flatten(), atol=1e-6)


if __name__ == '__main__':

    import logging
    logging.basicConfig(level=logging.INFO,
                        format='%(levelname)s : %(asctime)s : %(message)s')

    unittest.TextTestRunner(verbosity=2).run(suite())
