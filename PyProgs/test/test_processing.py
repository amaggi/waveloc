import unittest, os, glob
import numpy as np
from SDS_processing import do_SDS_processing_setup_and_run
from OP_waveforms import Waveform
from options import WavelocOptions


def suite():
  suite = unittest.TestSuite()
  suite.addTest(KurtosisTests('test_ss_kurtosis'))
  suite.addTest(ProcessingTests('test_positive_gradient'))
  suite.addTest(ProcessingTests('test_processing'))
  return suite

    
def waveforms_to_signature(base_path,datadir,dataglob,output_filename):

  sig_file=open(os.path.join(base_path,datadir,output_filename),'w')
  allfiles=glob.glob(os.path.join(base_path,datadir, dataglob))
  for filename in sorted(allfiles) :
    basename=os.path.basename(filename)
    wf=Waveform()
    wf.read_from_file(filename,format='MSEED')
    (maximum, datasum) = wf.compute_signature()
    sig_file.write("%s \t\t %.6f \t %.6f\n"%(basename,maximum,datasum))
 
class KurtosisTests(unittest.TestCase):

  def test_ss_kurtosis(self):
    import numpy as np
    from filters import sw_kurtosis1, sw_kurtosis2

    npts=1000
    nkurt=7
    r=np.random.randn(npts)
    s=np.zeros(npts)
    s[npts/2:]=10.0
    sig=s+r

    k1=sw_kurtosis1(sig,nkurt)
    k2=sw_kurtosis2(sig,nkurt)

    self.assertEquals(k1.shape,k2.shape)
    np.testing.assert_allclose(k1,k2,5)
    self.assertAlmostEquals(np.max(k1), np.max(k2))
    self.assertEquals(np.argmax(k1), np.argmax(k2))
   
    
   

class ProcessingTests(unittest.TestCase):

  def setUp(self):

    self.wo=WavelocOptions()
    self.wo.set_test_options()
    self.wo.verify_SDS_processing_options()

  def test_positive_gradient(self):
    from OP_waveforms import stream_positive_derivative
    from obspy.core import read

    base_path=self.wo.opdict['base_path']
    test_datadir=self.wo.opdict['test_datadir']

    st=read(os.path.join(base_path,test_datadir,'raw_data','YA.UV15.00.HHZ.MSEED'))
    tr=st[0]
    npts= len(tr.data)
    dt  = tr.stats.delta
    x = np.arange(npts)*dt

    # set up a polynomial function
    y = (3 + 2*x +4*x*x +5*x*x*x)
    dy_exp = (2 + 8*x +15*x*x)

    tr.data=y
    st=stream_positive_derivative(st,dt)
    np.testing.assert_almost_equal(tr.data[20:100], dy_exp[20:100],2)


  @unittest.expectedFailure
  @unittest.skip('Skip for now')
  def test_processing(self):

    base_path=self.wo.opdict['base_path']
    datadir=self.wo.opdict['datadir']
    test_datadir=self.wo.opdict['test_datadir']

    expected_signature_filename = os.path.join(base_path,test_datadir,'test_data_signature.dat')
    expected_signature_file = open(expected_signature_filename,'r') 
    expected_lines=expected_signature_file.readlines()

    self.wo.opdict['load_ttimes_buf']=False
    do_SDS_processing_setup_and_run(self.wo.opdict)
   
    waveforms_to_signature(base_path,os.path.join('data',datadir),'*mseed','data_signature.dat')
    signature_filename=os.path.join(base_path,'data',datadir,'data_signature.dat')
    signature_file = open(signature_filename,'r') 
    lines=signature_file.readlines()

    self.assertSequenceEqual(lines,expected_lines)

if __name__ == '__main__':

  import logging
  logging.basicConfig(level=logging.INFO, format='%(levelname)s : %(asctime)s : %(message)s')
 
  unittest.TextTestRunner(verbosity=2).run(suite())
 
