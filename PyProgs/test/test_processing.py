import unittest, os, glob
from waveloc.SDS_processing import do_SDS_processing_setup_and_run
from waveloc.OP_waveforms import Waveform
from waveloc.options import WavelocOptions

def suite():
  suite = unittest.TestSuite()
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
 

class ProcessingTests(unittest.TestCase):

  def setUp(self):

    self.wo=WavelocOptions()
    self.wo.set_test_options()
    self.wo.verify_SDS_processing_options()


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
 
