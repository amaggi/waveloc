import unittest
import os, glob
from SDS_processing import do_SDS_processing_setup_and_run
from OP_waveforms import Waveform

def suite():
  suite = unittest.TestSuite()
  suite.addTest(ProcessingTests('test_processing'))
  return suite

def setUpModule():

  datadir='TEST'
  net_list='YA'
  sta_list="FJS,FLR,FOR,HDL,RVL,SNE,UV01,UV02,UV03,UV04,UV05,UV06,UV07,UV08,UV09,UV10,UV11,UV12,UV13,UV14,UV15"
  comp_list="HHZ"
  starttime="2010-10-14T00:14:00.0Z"
  endtime="2010-10-14T00:18:00.0Z"
  resample=False
  fs=None
  c1=4.0
  c2=10.0
  kwin=4
  krec=False
  kderiv=True

#  do_SDS_processing_setup_and_run(
#     datadir=datadir,
#     net_list=net_list,
#     sta_list=sta_list,
#     comp_list=comp_list,
#     starttime=starttime,
#     endtime=endtime,
#     resample=resample,
#     fs=fs,
#     c1=c1,
#     c2=c2,
#     kwin=kwin,
#     krec=krec,
#     kderiv=kderiv)

  base_path=os.getenv('WAVELOC_PATH')

  sig_file=open(os.path.join(base_path,'test_data','test_data_signature.dat'),'w')
  allfiles=glob.glob(os.path.join(base_path,'data',datadir,'*mseed'))
  for filename in allfiles :
    basename=os.path.basename(filename)
    wf=Waveform()
    wf.read_from_file(filename,format='MSEED')
    (maximum, datasum) = wf.compute_signature()
    sig_file.write("%s \t\t %.6f \t %.6f\n"%(basename,maximum,datasum))
    


class ProcessingTests(unittest.TestCase):

  def test_processing(self):
    self.assertTrue(True)

if __name__ == '__main__':

#  import logging
#  logging.basicConfig(level=logging.DEBUG, format='%(levelname)s : %(asctime)s : %(message)s')
 
  unittest.TextTestRunner(verbosity=2).run(suite)
 
