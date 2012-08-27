import unittest
import os, glob
from SDS_processing import do_SDS_processing_setup_and_run
from OP_waveforms import Waveform

def suite():
  suite = unittest.TestSuite()
  suite.addTest(ProcessingTests('test_processing'))
  return suite

    
def waveforms_to_signature(base_path,datadir,dataglob,output_filename):

  sig_file=open(os.path.join(base_path,datadir,output_filename),'w')
  allfiles=glob.glob(os.path.join(base_path,datadir, dataglob))
  for filename in allfiles :
    basename=os.path.basename(filename)
    wf=Waveform()
    wf.read_from_file(filename,format='MSEED')
    (maximum, datasum) = wf.compute_signature()
    sig_file.write("%s \t\t %.6f \t %.6f\n"%(basename,maximum,datasum))
 

class ProcessingTests(unittest.TestCase):

  def setUp(self):

    self.base_path=os.getenv('WAVELOC_PATH')

    self.test_datadir='test_data'
    self.datadir='TEST'
    self.net_list='YA'
    self.sta_list="FJS,FLR,FOR,HDL,RVL,SNE,UV01,UV02,UV03,UV04,UV05,UV06,UV07,UV08,UV09,UV10,UV11,UV12,UV13,UV14,UV15"
    self.comp_list="HHZ"
    self.starttime="2010-10-14T00:14:00.0Z"
    self.endtime="2010-10-14T00:18:00.0Z"
    self.resample=False
    self.fs=None
    self.c1=4.0
    self.c2=10.0
    self.kwin=4
    self.krec=False
    self.kderiv=True


  def test_processing(self):

    expected_signature_filename = os.path.join(self.base_path,self.test_datadir,'test_data_signature.dat')
    expected_signature_file = open(expected_signature_filename,'r') 
    expected_lines=expected_signature_file.readlines()

#    do_SDS_processing_setup_and_run( datadir=self.datadir, net_list=self.net_list, sta_list=self.sta_list, comp_list=self.comp_list, starttime=self.starttime, endtime=self.endtime, resample=self.resample, fs=self.fs, c1=self.c1, c2=self.c2, kwin=self.kwin, krec=self.krec, kderiv=self.kderiv)
   
    waveforms_to_signature(self.base_path,os.path.join('data',self.datadir),'*mseed','data_signature.dat')
    signature_filename=os.path.join(self.base_path,'data',self.datadir,'data_signature.dat')
    signature_file = open(signature_filename,'r') 
    lines=signature_file.readlines()

    self.assertSequenceEqual(lines,expected_lines)

if __name__ == '__main__':

#  import logging
#  logging.basicConfig(level=logging.DEBUG, format='%(levelname)s : %(asctime)s : %(message)s')
 
  unittest.TextTestRunner(verbosity=2).run(suite)
 
