#!/usr/bin/env python
# encoding: utf-8

"""
The :mod:`OP_waveforms` module provides the :class:`Waveform` class through
which most time-series manipulation is carried out.

.. note:: 
   Requires **scikits.samplerate** for full functionality (on linux systems install with
   ``easy_install scikits.samplerate`` after having installed the
   ``libsamplerate0`` libraries using e.g. ``apt-get install libsamplerate0``). 

"""
import os, sys, glob 
import numpy as np
#import numexpr as ne
#import carray as ca
import logging

# this fakes the import of plt as it was originally
import matplotlib.pyplot as plt

from filters import *

from obspy.core import *
from obspy.signal import *
from obspyaux import *

import scipy.stats as ss

class Waveform(object):
  """
  Wrapper class for obspy streams. Adds many convenience functions, and
  functionality specific to WaveLoc processing.

  **Attributes**
  
  .. attribute:: stream

    An :class:`obspy.core.Stream` object.  Can be accessed directly.  
    Initialises to ``None``.
    If you set this attribute directly, you should also set :attr:`.trace`
    to :attr:`.stream`.traces[0]

  .. attribute:: trace

    An obspy.Trace object.  Can be accessed directly.  
    Initialises to ``None``.
    If you set this directly, make sure it points to :attr:`.stream`.traces[0] 

  .. attribute:: proc

    Text string indicating processing type.  Can be accessed directly.
    Initialises to ``None``.

  **Methods**
  """

  def __init__(self):
    """
    Initialises Waveform to have empty stream, trace and proc (they are all set to None).
    """
    self.stream=None
    self.trace=None
    self.proc=None

  def _get_npts_(self):
    """Read-only property.  Returns number of points in :attr:`trace`. """
    return self.trace.stats.npts
  def _get_delta_(self):
    """Read-only property.  Returns time step between data points in :attr:`trace`. """
    return self.trace.stats.delta
  def _get_station_(self):
    """Read-only property.  Returns station name from ``self.trace``. """
    return self.trace.stats.station
  def _get_channel_(self):
    """Read-only property.  Returns channel name from ``self.trace``. """
    return self.trace.stats.channel
  def _get_t_array_(self):
    """Read-only property.  Returns a ``numpy.array`` containing times in seconds from ``self.trace.stats.starttime`` for data points in ``self.trace``. """
    return np.arange(0,self.trace.stats.npts * self.trace.stats.delta, self.trace.stats.delta)
  def _get_values_(self):
    """Read-only property.  Returns values of ``self.trace.data``. """
    return self.trace.data


  npts=property(_get_npts_)
  delta=property(_get_delta_)
  dt=property(_get_delta_)
  station=property(_get_station_)
  channel=property(_get_channel_)
  comp=property(_get_channel_)
  t_array=property(_get_t_array_)
  values=property(_get_values_)

  def read_from_SDS(self,sds_root,net_name,sta_name,comp_name,starttime=None,endtime=None,rmean=False,taper=False,pad_value=None):
    """
    Read waveform data from an SDS structured archive.  Simple overlaps and
    adjacent traces are merged if possile.

    :param sds_root: root of the SDS archive
    :param net_name: network name
    :param sta_name: station name
    :param comp_name: component name
    :param starttime: Start time of data to be read.
    :param endtime: End time of data to be read.
    :param rmean: If ``True`` removes the mean from the data upon reading. If data are segmented, the mean will be removed from all segments individually.
    :param taper: If ``True`` applies a cosine taper to the data upon reading.  If data are segmented, tapers are applied to all segments individually.
    :param pad_value: If this parameter is set, points between ``starttime`` and the first point in the file, and points between the last point in the file and ``endtime``, will be set to ``pad_value``.  You may want to also use the ``rmean`` and ``taper`` parameters, depending on the nature of the data.

    :type starttime: ``obspy.core.utcdatetime.UTCDateTime`` object 
    :type endtime: ``obspy.core.utcdatetime.UTCDateTime`` object 
    :type rmean: boolean
    :type taper: boolean
    :type pad_value: float

    :raises UserWarning: If there are no data between ``starttime`` and ``endtime`` 

    """
    logging.info("Reading from SDS structure %s %s %s ..."%(net_name,sta_name,comp_name))
    filename=os.path.join(sds_root,net_name,sta_name,"%s.D"%comp_name,"*")
    if os.path.isdir(glob.glob(filename)[0]):
      filename=os.path.join(filename,"*")
    logging.debug(filename)
    st=stream.read(filename,starttime=starttime, endtime=endtime)
    st.merge(method=-1)

    if st.count()>1: # There are gaps after sensible cleanup merging
      logging.info("File contains gaps:")
      st.printGaps()

    if rmean:
      logging.info("Removing the mean from single traces.")
      st=stream_rmean(st)

    if taper:
      logging.info("Tapering single traces.")
      st=stream_taper(st)

    if not pad_value is None:
      try:

        first_tr=st.traces[0]
        # save delta (to save typing)
        delta=first_tr.stats.delta
        if (not starttime is None) and ((first_tr.stats.starttime - starttime) > delta):
          logging.debug("Padding with value %f from %s to first point in file at %s."%(pad_value, starttime.isoformat(), first_tr.stats.starttime.isoformat()))
          # find the number of points from starttime to end of the first trace
          npts_full_trace=int(np.floor((first_tr.stats.endtime-starttime)/delta))+1
          # find the number of points of the padding section
          n_pad=npts_full_trace-first_tr.stats.npts
          # fill the full time range with padd value
          tr_pad=np.zeros(npts_full_trace)+pad_value
          # substitute in the data
          tr_pad[n_pad:]=first_tr.data[:]
          first_tr.data=tr_pad
          first_tr.stats.starttime=starttime
          first_tr.stats.npts=npts_full_trace
          st.traces[0]=first_tr

        last_tr=st.traces[-1]
        # save delta (to save typing)
        delta=last_tr.stats.delta
        if (not endtime is None) and ((endtime - last_tr.stats.endtime) > delta):
          logging.debug("Padding with value %f from last point in file at %s to %s."%(pad_value, last_tr.stats.endtime.isoformat(), endtime.isoformat()))
          # find the number of points from starttime to end of the first trace
          npts_full_trace=int(np.floor((endtime-last_tr.stats.starttime)/delta))+1
          # fill the full time range with padd value
          tr_pad=np.zeros(npts_full_trace)+pad_value
          # substitute in the data
          tr_pad[0:last_tr.stats.npts]=last_tr.data[:]
          last_tr.data=tr_pad
          last_tr.stats.npts=npts_full_trace
          st.traces[-1]=last_tr


      except IndexError:
        logging.warning('No data within time limits requested')
        raise UserWarning('No data within time limits requested.')


    try:
      self.stream=st
      self.trace=st.traces[0]
      self.proc="None"
    except IndexError:
      raise UserWarning('No data within time limits requested.')



  def read_from_file(self,filename, format=None, starttime=None, endtime=None, rmean=False, taper=False, pad_value=None): 
    """
    Read waveform data from file.  Multiple traces are merged if they overlap exactly or are adjacent.

    :param filename: Waveform filename
    :param format: ``obspy`` format type (e.g. 'SAC', 'mseed'...)
    :param starttime: Start time of data to be read.
    :param endtime: End time of data to be read.
    :param rmean: If ``True`` removes the mean from the data upon reading. If data are segmented, the mean will be removed from all segments individually.
    :param taper: If ``True`` applies a cosine taper to the data upon reading.  If data are segmented, tapers are applied to all segments individually.
    :param pad_value: If this parameter is set, points between ``starttime`` and the first point in the file, and points between the last point in the file and ``endtime``, will be set to ``pad_value``.  You may want to also use the ``rmean`` and ``taper`` parameters, depending on the nature of the data.

    :type format: string 
    :type starttime: ``obspy.core.utcdatetime.UTCDateTime`` object 
    :type endtime: ``obspy.core.utcdatetime.UTCDateTime`` object 
    :type rmean: boolean
    :type taper: boolean
    :type pad_value: float

    :raises UserWarning: If there are no data between ``starttime`` and ``endtime`` 

    """
    logging.debug("Reading from %s..."%filename)
    if format is not None:
      st=stream.read(filename,format,starttime=starttime, endtime=endtime)
    else:
      st=stream.read(filename,starttime=starttime, endtime=endtime)

    st.merge(method=-1) 

    if st.count()>1: # There are gaps after intelligent merge
      logging.info("File contains gaps:")
      st.printGaps()

    if rmean:
      st=stream_rmean(st)

    if taper:
      st=stream_taper(st)

    if not pad_value is None:
      try:

        first_tr=st.traces[0]
        # save delta (to save typing)
        delta=first_tr.stats.delta
        if (not starttime is None) and ((first_tr.stats.starttime - starttime) > delta):
          logging.debug("Padding with value %f from %s to first point in file at %s."%(pad_value, starttime.isoformat(), first_tr.stats.starttime.isoformat()))
          # find the number of points from starttime to end of the first trace
          npts_full_trace=int(np.floor((first_tr.stats.endtime-starttime)/delta))+1
          # find the number of points of the padding section
          n_pad=npts_full_trace-first_tr.stats.npts
          # fill the full time range with padd value
          tr_pad=np.zeros(npts_full_trace)+pad_value
          # substitute in the data
          tr_pad[n_pad:]=first_tr.data[:]
          first_tr.data=tr_pad
          first_tr.stats.starttime=starttime
          first_tr.stats.npts=npts_full_trace
          st.traces[0]=first_tr

        last_tr=st.traces[-1]
        # save delta (to save typing)
        delta=last_tr.stats.delta
        if (not endtime is None) and ((endtime - last_tr.stats.endtime) > delta):
          logging.debug("Padding with value %f from last point in file at %s to %s."%(pad_value, last_tr.stats.endtime.isoformat(), endtime.isoformat()))
          # find the number of points from starttime to end of the first trace
          npts_full_trace=int(np.floor((endtime-last_tr.stats.starttime)/delta))+1
          # fill the full time range with padd value
          tr_pad=np.zeros(npts_full_trace)+pad_value
          # substitute in the data
          tr_pad[0:last_tr.stats.npts]=last_tr.data[:]
          last_tr.data=tr_pad
          last_tr.stats.npts=npts_full_trace
          st.traces[-1]=last_tr


      except IndexError:
        logging.warning('No data within time limits requested')
        raise UserWarning('No data within time limits requested.')


    try:
      self.stream=st
      self.trace=st.traces[0]
      self.proc="None"
    except IndexError:
      raise UserWarning('No data within time limits requested.')


  def read_from_file_notaper(self,filename, format=None, starttime=None, endtime=None, rmean=False, pad_value=None): 
    """
    Read waveform data from file without any tapering.  Internally calls :func:`read_from_file` with ``taper=False``. 

    :param filename: Waveform filename
    :param format: ``obspy`` format type (e.g. 'SAC', 'mseed'...)
    :param starttime: Start time of data to be read.
    :param endtime: End time of data to be read.
    :param rmean: If ``True`` removes the mean from the data upon reading. If data are segmented, the mean will be removed from all segments individually.
    :param pad_value: If this parameter is set, points between ``starttime`` and the first point in the file, and points between the last point in the file and ``endtime``, will be set to ``pad_value``.  You may want to also use the ``rmean`` parameter, depending on the nature of the data.

    :type format: string 
    :type starttime: ``obspy.core.utcdatetime.UTCDateTime`` object 
    :type endtime: ``obspy.core.utcdatetime.UTCDateTime`` object 
    :type rmean: boolean
    :type pad_value: float

    .. note::
      This method is deprecated, and will put a warning message on the log to remind you not to use it.  Use :func:`read_from_file` instead.  

    """

    logging.warning("DEPRECATED!  Use read_from_file(...taper=False).")
    self.read_from_file(filename,format=format, starttime=starttime, endtime=endtime,pad_value=pad_value,rmean=rmean,taper=False)


  def cleanup_traces(self,short_data_length=0,rmean=False,taper=False):
    """
    Clean up the stream by removing isolated segments less than ``short_data_length`` seconds in length.
    Perfectly overlapping segments are merged before removal of short isolated segments begins.

    :param short_data_length: Data segments shorter than this number of seconds will be discarded.
    :param rmean: If ``True`` removes the mean from the data after discarding
                  short segments. If data are segmented, the mean will be removed from all
                  segments individually.
    :param taper: If ``True`` removes the mean from the data after discarding
                  short segments. If data are segmented, a taper will be applied to all
                  segments individually.
    :type rmean: boolean
    :type taper: boolean

    """
    # short_data_length = length of data in seconds to discard

    # merge perfectly overlapping or adjacent data
    stream.merge(method=-1)
  
    clean_stream=self.stream.copy()
    clean_stream.clear()
    gaps=self.stream.getGaps()
    self.stream.printGaps()

    for tr in self.stream.traces:
       length=tr.stats.npts * tr.stats.delta
       if length> short_data_length:
         tr2=tr.copy()
         clean_stream.append(tr2)

    clean_stream.printGaps()
    
    if rmean:
      clean_stream.rmean()

    if taper:
      clean_stream.taper()
    
    self.stream=clean_stream
    try:
      self.trace=self.stream[0]
    except IndexError:
      raise UserWarning ('No data left after cleaning')


  def write_to_file_filled(self,filename,format=None,fill_value=0,rmean=False,taper=False):
    """
    Write waveform to file, after merging and filling blanks with ``fill_value``.  

    :param filename: Output filename.
    :param format: ``obspy`` format type (e.g. 'SAC', 'mseed'...)
    :param fill_value: Value used to fill in gaps
    :param rmean: Remove the mean before merging.
    :param taper: Apply taper before merging.

    """
    logging.info("Merging traces before writing file %s\n"%filename)
    # make a copy and write that, so as not to merge the file in memory
    st=self.stream.copy()
    if rmean:
      logging.info("Removing mean before merging.")
      st.rmean()
    if taper:
      logging.info("Applying taper before merging.")
      st.taper()
    st.merge(method=1,fill_value=fill_value)
    for tr in st:
#      tr.data=tr.data.astype("int32")
      tr.data=tr.data.astype("float32")
      
    st.write(filename,format)

  def write_to_file(self,filename,format=None):
    """
    Write waveform to file.  

    :param filename: Output filename.
    :param format: ``obspy`` format type (e.g. 'SAC', 'mseed'...)

    """
 
    for tr in st:
      tr.data=tr.data.astype("int32")
    self.stream.write(filename,format)

  def rmean(self):
    """
    Removes the mean of the stream (iterates over all traces).
    """
 
    self.stream=stream_rmean(self.stream)
    try:
      self.trace=self.stream[0]
    except IndexError:
      raise UserWarning('No data in stream for rmean')

  def taper(self):
    """
    Applies a cosine taper the stream (iterates over all traces).
    """
    self.stream=stream_taper(self.stream)
    try:
      self.trace=self.stream[0]
    except IndexError:
      raise UserWarning('No data in stream for taper')

  def display(self,title="",filename=""):
    """
    Makes a quick and dirty plot of the waveform.
    If filename is given (and is different from "") then the plot of the
    waveform is written to file, otherwise it is shown on the screen.  

    :param title: Title for the plot.
    :param filename: Filename for the plot (format defined by the file extension).
    
    """
    
    plt.clf()
    
    plt.title(title)

    # set the axis labels
    plt.xlabel("Time / s")
    if self.proc=="None":
      plt.ylabel("Raw data")
    elif self.proc=="StaLta":
      plt.ylabel("STA/LTA")
    elif self.proc=="Envelope":
      plt.ylabel("Envelope")
    elif self.proc=="Skewness":
      plt.ylabel("Absolute value Skewness")
    elif self.proc=="Kurtosis":
      plt.ylabel("Kurtosis")
    else:
      plt.ylabel("")
      
    # plot the waveform
    plt.plot(self.t_array,self.values)
    
    # save to file or display to screen
    if not filename=="":
      plt.savefig(filename)
    else:
      plt.show()

  def bp_filter(self,freqmin,freqmax,zerophase=False,rmean=False,taper=False):
    """
    Apply a band-pass filter to the data.  If data are segmented into multiple
    traces, apply the same filter to all traces.
    Calls :func:`obspy.signal.filter` to do the filtering.

    :param freqmin: Low frequency corner
    :param freqmax: High frequency corner
    :param zerophase: If ``True`` applies a non-causal bandpass filter.  If
                      ``False`` applies a causal bandpass filter.
    :param rmean: If ``True`` remove mean before filtering.
    :param taper: If ``True`` apply taper before filtering.

    :type zerophase: boolean
    :type rmean: boolean
    :type taper: boolean

    
    """
    if zerophase:
      logging.info("Non-causal band-pass filtering single traces : %.2fHz - %.2fHz\n"%(freqmin,freqmax))
    else:
      logging.info("Causal band-pass filtering single traces : %.2fHz - %.2fHz\n"%(freqmin,freqmax))

    if rmean:
      self.rmean()
    if taper:
      self.taper()

    for itr in range(self.stream.count()) :
      tr=self.stream.traces[itr]
      xs=filter.bandpass(tr.data,freqmin,freqmax,tr.stats.sampling_rate, zerophase=zerophase)
      tr.data=xs  
      self.stream.traces[itr]=tr

    try:
      self.trace=self.stream.traces[0]
    except IndexError:
      raise UserWarning ('No data in stream at bp_filter.')

  def resample(self,new_samplerate,resample_type='sinc_best'):
    """
    Applies audio-quality resampling in place.  Requires installation of ``scikits.samplerate``.

    :param new_samplerate: New sample rate.
    :param resample_type: Can be ``'sinc_best'``, ...
    """
    try:
      from scikits.samplerate import resample as sci_resample
      old_samplerate=1/np.float(self.delta)
      ratio=new_samplerate/old_samplerate
      for itr in range(self.stream.count()) :
        tr=self.stream.traces[itr]
        xs=sci_resample(tr.data,ratio,resample_type,verbose=True)
        tr.data=xs  
        tr.stats.sampling_rate=new_samplerate
        self.stream.traces[itr]=tr

    except ImportError:
     logging.warn('Cannot import scikits.samlerate.resample')
     


  def decimate(self,factor=1):
    """
    Applies ``obspy`` decimation, after applying an anti-aliasing, non-causal
    filter of 0.4*(new sample rate).

    :param factor: Decimation factor.
    :type factor: integer
    """
    self.trace.filter('lowpass', freq=0.4*self.trace.stats.sampling_rate/float(factor), zerophase=True) 
    self.trace.downsample(decimation_factor=factor, strict_length=False, no_filter=True)

  def get_snr(self,o_time,left_time,right_time):
    """
    Returns signal-to-noise ratio, where signal = max(abs(signal between left_time and right_time)) and noise = median(abs(signal between left_time- and otime)).
    """
    tr_signal=self.trace.slice(left_time,right_time)
    signal_value=np.max(np.abs(tr_signal.data))

#    signal_var=np.std(tr_signal.data)
    tr_noise=self.trace.slice(left_time,o_time)
    noise_value=np.median(np.abs(tr_noise.data))

    if noise_value==0.0:
      snr=0
    else:
      snr=signal_value / noise_value

    return snr

  def compute_signature(self):
    maximum=np.max(self.trace.data)
    datasum=np.sum(self.trace.data)
    return (maximum, datasum)

  def process_envelope(self):
    """
    Runs envelope processing on a waveform.
    """
    xs=filter.envelope(self.values)
    self.trace.data=xs
    self.proc="Envelope"

  def process_none(self):
    """
    No processing on a waveform. 
    """
    self.proc="None"
  
    
  def process_sta_lta(self,stawin,ltawin):
    """
    Runs classic short to long term average ratio processing on waveform.
    
    :param stawin: length of the short term average window in seconds
    :param ltawin: length of the long term average window in seconds
    
    """
    nsta=int(stawin/self.delta)
    nlta=int(ltawin/self.delta)
    xs=trigger.classicStaLta(self.values,nsta,nsta)
    self.trace.data=xs
    self.proc='StaLta'


  def process_skewness(self,Tpad,win):
    """
    Processing waveform using skewness (from statlib package).
    
    Calls filters.skewness(), and overwrites the waveform.
    
    :param Tpad: padding to add on to the start of the waveform if necessary
    :param win:  length of the window on which to calculate the skewness in
                 seconds
    
    """
    # Remove the initial mean, then taper, padd with noise 
    dt=self.dt
    x=self.values

    # factor necessary to replicate as much as possible the initial noise
    # this is determined from the variance of a  box-car distribution of random
    # values
    f_std=3.48

    # chose the initial window length on which detyermine the
    # standard deviation 
    ini_win=3.0
    ini_std=initial_std(x,dt,ini_win)

    # shift the time series accordingly
    x=rmean_ini(self.values,dt,ini_win)

    # when padding with noise the tapering is not necessary 
    x = noise_padding(x,dt,Tpad,ini_std * f_std)
    xs=skewness(x,dt,win)
    
    # At this point, xs starts at point b-Tpad+win
    # so fix b, values and npts accordingly
    self.b=self.b-Tpad+win
    self.npts=len(xs)

    # first and last points are problematic for kurtosis
    # (often contain spikes)
    # so get rid of them to be consistent with kurtosis
    xs.pop()    # remove last point
    xs.pop()    # remove last point
    xs.pop(0) # remove first point
    xs.pop(0) # remove first point
    self.b=self.b+2*self.dt  # shunt b up two time steps
    self.npts=len(xs)      # recalculate number of points
    
    # Save xs values as waveform and set processing type.
    self.values=numpy.array(xs)
    self.proc='Skewness'

  def take_positive_derivative(self, pre_rmean=False, pre_taper=False, post_taper=True):
    """
    Takes positive derivative of a stream
    """

    if pre_rmean:
      self.rmean()

    if pre_taper:
      self.taper()

    self.stream=stream_positive_derivative(self.stream)
    try:
      self.trace=self.stream[0]
    except IndexError:
      raise UserWarning('No data in stream for positive_derivative')

    if post_taper:
      self.taper()

    
  def process_kurtosis(self,win, recursive=False, pre_rmean=False, pre_taper=False, post_taper=True):
    """
    Processing waveform using kurtosis (from statlib package).
    
    Calls filters.kurto(), and overwrites the waveform.
    
    :param win: length of the window (in seconds) on which to calculate the
                kurtosis 
    :param pre_rmean: If ``True`` removes mean of signal before processing.
    :param pre_taper: If ``True`` applies taper to signal before processing.
    :param post_taper: If ``True`` applies taper to signal after processing.

    :type pre_rmean: boolean
    :type pre_taper: boolean
    :type post_taper: boolean
    
    """
   
    logging.info("Applying kurtosis to single traces, window = %.2f s\n"%win)

    dt=self.dt


    if pre_rmean:
      self.rmean()

    if pre_taper:
      self.taper()

    # process each trace independently
    for itr in range(self.stream.count()):
      tr=self.stream.traces[itr]
      x=tr.data
      varx=np.std(x)

      npts=len(tr.data)

      xs=np.zeros(npts)    

      if recursive:
        mean_value=0
        var_value=0
        kurt_value=0
        C=1-dt/win
        for i in range(npts):
          mean_value = C*mean_value + (1-C)*x[i]
          var_value=C*var_value+(1-C)*(x[i]-mean_value)**2
          if var_value>varx: 
            kurt_value=C*kurt_value+(1-C)*(x[i]-mean_value)**4/var_value**2
          else : 
            kurt_value=C*kurt_value+(1-C)*(x[i]-mean_value)**4/varx**2
          xs[i]=kurt_value-3

      else:
        nwin=int(win/dt)+1
        for i in range(nwin,npts - nwin ):
          xs[i+nwin]=ss.kurtosis(x[i:(i+nwin)])

        
      #xs_filt=lowpass(xs,10*tr.stats.delta,1/tr.stats.delta,zerophase=True)
      xs_filt=smooth(xs)
       
      # Save xs values as waveform 
      tr.data=xs_filt

      # put trace back into stream
      self.stream.traces[itr]=tr

     
    # apply taper after kurtosis calculation if required
    if post_taper:
      self.taper()

    # set the process flag
    self.proc='Kurtosis'

#########################
# Functions
#########################
def stream_rmean(st):
  """
  Removes mean from a stream (iterates over all available traces).
  """
  for tr in st:
    t_tr=(tr.data-np.mean(tr.data))
    tr.data=t_tr
  return st

def stream_positive_derivative(st):
  """
  Takes first time derivative of a stream (iterates over all available traces) and keep only positive values (set negative values to zero).
  """
  for tr in st:
    xs=tr.data
    try:
      xtemp=np.gradient(xs)
      for i in range(len(xtemp)):
        if xtemp[i]<0 :
          xtemp[i]=0 
      xs=xtemp
    except IndexError:
      logging.warn('Zero length data segment')
    tr.data=xs
  return st

def stream_taper(st):
  """
  Applies cosine taper to a stream (iterates over all available traces).
  """
  for tr in st:
    try:
      mytaper=cosTaper(tr.stats.npts)
      t_tr=mytaper*(tr.data)
    except ValueError:
      logging.warn('Trace is too short for tapering - multiplying by 0 instead')
      t_tr=0.0*(tr.data)
    tr.data=t_tr
  return st

def read_data_compatible_with_time_dict(filenames, time_dict, starttime, endtime):

  data={}
  deltas=[]

  for filename in filenames:

    st=read(filename,headonly=True)

    delta=st[0].stats.delta
    deltas.append(delta)

    wf_id=st[0].stats.station
    
    if time_dict.has_key(wf_id):
      try:
        wf=Waveform()
        wf.read_from_file(filename,starttime=starttime,endtime=endtime,pad_value=0)
        # just put the bare data into the data dictionary
        data[wf_id]=wf.values
      except UserWarning,msg:
        logging.warn("No data data found between limits for file %s. Ignoring station %s."%(filename,wf_id))
    else:
      logging.warn('Station %s not present in time_grid.  Ignoring station %s.'%(wf_id,wf_id))

  # cheque uniqueness of delta
  u=np.unique(np.array(deltas))
  if len(u) > 1 :
    logging.error('Sampling frequency differs between stations.  Fix this before migrating.')
    for i in xrange(len(deltas)):
      logging.error('Delta %.2f for file %s'%(deltas[i],filenames[i]))
    raise UserWarning
  
  return data, u[0]

def process_all_data_kurtosis(files, start_time, end_time, filter_c1, filter_c2, kurt_window):
  """
  Applies band pass filter and kurtosis to a list of files.  

  .. warning : 
    Deprecated! Will disappear soon 
  """

  logging.warn('process_all_data_kurtosis is depracated.  Will disappear soon.')
  
  for file in files:
    wf=Waveform()
    file_filtered="%s.filt.sac"%file
    file_kurtosis="%s.filt_kurt.sac"%file
    if options.verbose:
      print file
    wf.read_from_file(file,'MSEED',starttime=start_time,endtime=end_time)
    wf.bp_filter(filter_c1,filter_c2,rmean=True,taper=True)
    wf.write_to_file(file_filtered, format='SAC')
    wf.process_kurtosis(kurt_window,post_taper=True)
    wf.write_to_file_filled(file_kurtosis, format='SAC')
 




