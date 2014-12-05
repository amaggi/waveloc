import unittest
from waveloc.correlation import correlate
import numpy as np


def suite():
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(CorrelationTests))
    return suite


class CorrelationTests(unittest.TestCase):

    def setUp(self):

        # Sweep
        from math import pi
        t = np.arange(0, 10, 0.01)
        a = 2
        b = 0.5
        shift = 100
        self.x = np.cos(2*pi*(a+b*t)*t)
        self.y = np.zeros(len(self.x))
        self.y[shift:] = 10*self.x[:-shift]

        self.dt = 0.01
        self.expected_tau = -shift*self.dt

    def test_auto_corr_time(self):
        self.v = 0
        dtime, cval = correlate(self.x, self.x, self.dt, self.v, 't')
        self.assertAlmostEqual(cval, 1)
        self.assertAlmostEqual(dtime, 0)

    def test_auto_corr_fr(self):
        self.v = 0
        dtime = correlate(self.x, self.x, self.dt, self.v, 'f')
        self.assertEqual(dtime, 0)

    def test_correlation(self):
        self.v = 0
        dtime_t, cval = correlate(self.x, self.y, self.dt, self.v, 't')
        ndtime = int(round(dtime_t*1./self.dt))
        yy = np.zeros(len(self.y))
        yy[:ndtime] = self.y[-ndtime:]
        dtime_f = correlate(self.x, yy, self.dt, self.v, 'f')
        dtime = dtime_t+dtime_f

        self.assertAlmostEqual(dtime, self.expected_tau, 2)


if __name__ == '__main__':
    unittest.TextTestRunner(verbosity=2).run(suite())
