import unittest
import os
from SDS_processing import do_SDS_processing_setup_and_run

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

  do_SDS_processing_setup_and_run(
     datadir=datadir,
     net_list=net_list,
     sta_list=sta_list,
     comp_list=comp_list,
     starttime=starttime,
     endtime=endtime,
     resample=resample,
     fs=fs,
     c1=c1,
     c2=c2,
     kwin=kwin,
     krec=krec,
     kderiv=kderiv)



class ProcessingTests(unittest.TestCase):

  def test_processing(self):
    self.assertTrue(True)

if __name__ == '__main__':

#  import logging
#  logging.basicConfig(level=logging.DEBUG, format='%(levelname)s : %(asctime)s : %(message)s')
 
  unittest.TextTestRunner(verbosity=2).run(suite)
 
