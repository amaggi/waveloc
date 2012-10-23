import os, glob, unittest
from options import WavelocOptions
from OP_waveforms import Waveform
from double_diff import *
import numpy as np

def suite():
  suite = unittest.TestSuite()
  suite.addTest(SyntheticsDoubleDiffTests('test_DD'))
  #suite.addTest(DoubleDiffTests(''))
  return suite


class SyntheticsDoubleDiffTests(unittest.TestCase):

  def setUp(self):
    # Let's consider a 100x100x10 km grid with an homogeneous velocity Vp=6 km/s
    self.Vp=6

    # Seismic stations are all at z=0km, and placed every 2 km in x and y directions, from 0 to 10 km
    # positive z axis upwards
    self.sta={}
    self.sta[1]={'x':0,'y':0,'depth':0,'elev':0,'station':1}
    self.sta[2]={'x':40,'y':0,'depth':0,'elev':0,'station':2}
    self.sta[3]={'x':80,'y':0,'depth':0,'elev':0,'station':3}
    self.sta[4]={'x':20,'y':20,'depth':0,'elev':0,'station':4}
    self.sta[5]={'x':60,'y':20,'depth':0,'elev':0,'station':5}
    self.sta[6]={'x':100,'y':20,'depth':0,'elev':0,'station':6}
    self.sta[7]={'x':0,'y':40,'depth':0,'elev':0,'station':7}
    self.sta[8]={'x':40,'y':40,'depth':0,'elev':0,'station':8}
    self.sta[9]={'x':80,'y':40,'depth':0,'elev':0,'station':9}
    self.sta[10]={'x':20,'y':60,'depth':0,'elev':0,'station':10}
    self.sta[11]={'x':60,'y':60,'depth':0,'elev':0,'station':11}
    self.sta[12]={'x':100,'y':60,'depth':0,'elev':0,'station':12}
    self.sta[13]={'x':0,'y':80,'depth':0,'elev':0,'station':13}
    self.sta[14]={'x':40,'y':80,'depth':0,'elev':0,'station':14}
    self.sta[15]={'x':80,'y':80,'depth':0,'elev':0,'station':15}
    self.sta[16]={'x':20,'y':100,'depth':0,'elev':0,'station':16}
    self.sta[17]={'x':60,'y':100,'depth':0,'elev':0,'station':17}
    self.sta[18]={'x':100,'y':100,'depth':0,'elev':0,'station':18}
    self.sta[19]={'x':0,'y':0,'depth':0,'elev':-5,'station':19}
    self.sta[20]={'x':40,'y':0,'depth':0,'elev':-5,'station':20}
    self.sta[21]={'x':80,'y':0,'depth':0,'elev':-5,'station':21}
    self.sta[22]={'x':20,'y':20,'depth':0,'elev':-5,'station':22}
    self.sta[23]={'x':60,'y':20,'depth':0,'elev':-5,'station':23}
    self.sta[24]={'x':100,'y':20,'depth':0,'elev':-5,'station':24}
    self.sta[25]={'x':0,'y':40,'depth':0,'elev':-5,'station':25}
    self.sta[26]={'x':40,'y':40,'depth':0,'elev':-5,'station':26}
    self.sta[27]={'x':80,'y':40,'depth':0,'elev':-5,'station':27}
    self.sta[28]={'x':20,'y':60,'depth':0,'elev':-10,'station':28}
    self.sta[29]={'x':60,'y':60,'depth':0,'elev':-10,'station':29}
    self.sta[30]={'x':100,'y':60,'depth':0,'elev':-10,'station':30}
    self.sta[31]={'x':0,'y':80,'depth':0,'elev':-10,'station':31}
    self.sta[32]={'x':40,'y':80,'depth':0,'elev':-10,'station':32}
    self.sta[33]={'x':80,'y':80,'depth':0,'elev':-10,'station':33}
    self.sta[34]={'x':20,'y':100,'depth':0,'elev':-10,'station':34}
    self.sta[35]={'x':60,'y':100,'depth':0,'elev':-10,'station':35}
    self.sta[36]={'x':100,'y':100,'depth':0,'elev':-10,'station':36}

    self.area=[0,100,0,100,-10,0]

    # Let's assume 5 seismic events occurring at the same place (x=50,y=50,z=-5) but not at the same time
    self.cluster=[1,2,3,4,5]
    self.N=len(self.cluster)

    # Define true hypocentral parameters
    # positive z axis downwards
    self.locs_true=[]
    self.locs_true.append({'x_mean':50,'y_mean':50,'z_mean':5,'o_time':utcdatetime.UTCDateTime('2010-01-01T12:00:00.0000Z')})
    self.locs_true.append({'x_mean':50,'y_mean':50,'z_mean':5,'o_time':utcdatetime.UTCDateTime('2010-01-01T12:00:00.0000Z')})
    self.locs_true.append({'x_mean':50,'y_mean':50,'z_mean':5,'o_time':utcdatetime.UTCDateTime('2010-01-01T12:00:00.0000Z')})
    self.locs_true.append({'x_mean':50,'y_mean':50,'z_mean':5,'o_time':utcdatetime.UTCDateTime('2010-01-01T12:00:00.0000Z')})
    self.locs_true.append({'x_mean':50,'y_mean':50,'z_mean':5,'o_time':utcdatetime.UTCDateTime('2010-01-01T12:00:00.0000Z')})

    # Measured hypocentral parameters
    # positive z-axis downwards
    err_x=[0,0,0,0,0]
    err_y=[0,0,0,0,0]
    err_z=[0,0,0,0,0]
    err_to=[0,0,0,0,0]

    err_x=[0.2,0.3,-0.2,-0.3,0]
    err_y=[-0.3,-0.1,0.1,0.4,-0.1]
    err_z=[-0.5,-0.25,0.25,0.5,0]
    err_to=[2,4,-2,1,-4]

    self.locs_mes=[]
    for i in range(len(self.locs_true)):
      self.locs_mes.append({'x_mean':self.locs_true[i]['x_mean']+err_x[i],'y_mean':self.locs_true[i]['y_mean']+err_y[i],'z_mean':self.locs_true[i]['z_mean']+err_z[i],'o_time':self.locs_true[i]['o_time']+err_to[i]})

    # Input parameters
    self.threshold=0.8
    self.nbmin=3

    # Compute the traveltimes and arrival times
    self.ttimes_true={}
    self.atimes_true={}
    self.ttimes_mes={}
    self.atimes_mes={}
    for staname in self.sta.keys():
      xsta=self.sta[staname]['x']
      ysta=self.sta[staname]['y']
      zsta=-self.sta[staname]['elev'] # positive z-axis downwards
      self.ttimes_true[staname]=[]
      self.atimes_true[staname]=[]
      self.ttimes_mes[staname]=[]
      self.atimes_mes[staname]=[]
      for j in range(self.N):
        d_true=np.sqrt((xsta-self.locs_true[j]['x_mean'])**2+(ysta-self.locs_true[j]['y_mean'])**2+(zsta-self.locs_true[j]['z_mean'])**2)
        self.ttimes_true[staname].append(d_true/self.Vp)
        self.atimes_true[staname].append(self.locs_true[j]['o_time']+self.ttimes_true[staname][j])
        d_mes=np.sqrt((xsta-self.locs_mes[j]['x_mean'])**2+(ysta-self.locs_mes[j]['y_mean'])**2+(zsta-self.locs_mes[j]['z_mean'])**2)
        self.ttimes_mes[staname].append(d_mes/self.Vp)
        self.atimes_mes[staname].append(self.locs_mes[j]['o_time']+self.ttimes_mes[staname][j])

    self.coeff={}
    self.delay={}
    for staname in self.sta.keys():
      self.coeff[staname]=np.zeros((self.N,self.N))
      up_tr=np.triu_indices(self.N)
      self.coeff[staname][up_tr]=1
      self.delay[staname]=np.zeros((self.N,self.N))
      for i in range(self.N):
        for j in range(i+1,self.N):
          self.delay[staname][i][j]=self.ttimes_true[staname][i]-self.ttimes_true[staname][j]+err_to[j]-err_to[i]

    self.locs_expected=[]
    for i in range(len(self.locs_true)):
      self.locs_expected.append({'x_mean':self.locs_true[i]['x_mean']+np.mean(err_x),'y_mean':self.locs_true[i]['y_mean']+np.mean(err_y),'z_mean':self.locs_true[i]['z_mean']+np.mean(err_z),'o_time':self.locs_true[i]['o_time']+np.mean(err_to)})


  def test_DD(self):

    x,y,z,z_ph,to = coord_cluster(self.cluster,self.locs_mes)

    xcal,ycal,zcal,tocal = do_double_diff(x,y,z,to,self.sta,self.coeff,self.delay,self.cluster,self.threshold,self.ttimes_mes,self.atimes_mes)

    locs_cal=[]
    for i in range(len(xcal)):
      self.assertAlmostEqual(xcal[i],self.locs_expected[i]['x_mean'],2)
      self.assertAlmostEqual(ycal[i],self.locs_expected[i]['y_mean'],2)
      self.assertAlmostEqual(zcal[i],self.locs_expected[i]['z_mean'],2)
      self.assertAlmostEqual(tocal[i],self.locs_expected[i]['o_time'],2)
      locs_cal.append({'x_mean':xcal[i],'y_mean':ycal[i],'z_mean':zcal[i],'o_time':tocal[i]})

if __name__ == '__main__':

  import logging
  logging.basicConfig(level=logging.INFO, format='%(levelname)s : %(asctime)s : %(message)s')
 
  unittest.TextTestRunner(verbosity=2).run(suite())

