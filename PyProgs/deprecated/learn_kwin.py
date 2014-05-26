import os
import numpy as np
from copy import deepcopy
from obspy.core import utcdatetime
from OP_waveforms import Waveform
from scipy.signal import correlate
import logging

logging.basicConfig(level=logging.INFO, format='%(levelname)s : %(asctime)s : %(message)s')

base_path=os.getenv('WAVELOC_PATH')
test_file=os.path.join(base_path,'test_data','raw_data','YA.UV15.00.HHZ.MSEED')

starttime=utcdatetime.UTCDateTime("2010-10-14T00:14:00.0Z")
endtime=utcdatetime.UTCDateTime("2010-10-14T00:18:00.0Z")

c1=4.0
c2=10.0
kwin=3.0

logging.info('Computing full kurtosis kwin=%.2fs'%kwin)
wf_raw=Waveform()
wf_raw.read_from_file(test_file,starttime=starttime,endtime=endtime,rmean=True,taper=True)
wf_raw.bp_filter(freqmin=c1, freqmax=c2)

wf_kurt=deepcopy(wf_raw)
wf_kurt.process_kurtosis(kwin)

tr1=np.array(wf_kurt.trace.data)
tr1 = tr1/np.max(tr1)
tr_len=len(tr1)

misfits=[]
k_factors=np.linspace(0.01,0.2,30)
print k_factors
for k_factor in k_factors:
  k_rec=kwin*k_factor

  logging.info('Computing recursive kurtosis kwin=%.2fs'%k_rec)
  wf_krec=Waveform()
  wf_krec.read_from_file(test_file,starttime=starttime,endtime=endtime,rmean=True,taper=True)
  wf_krec.bp_filter(freqmin=c1, freqmax=c2)
  wf_krec.process_kurtosis(k_rec,recursive=True)

  tr2=np.array(wf_krec.trace.data)
  tr2 = tr2/np.max(tr2)

  logging.info('Computing misfit...')
  misfit=np.dot((tr1-tr2),(tr1-tr2))
  logging.info('misfit = %.3f'%misfit)
  misfits.append((k_factor,misfit))

print misfits
