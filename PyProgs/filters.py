#!/usr/bin/env python
# encoding: utf-8

"""
Functions for timeseries filtering.

Created by Alessia Maggi and Alberto Michelini.
"""

import numpy as np
from math import *
from numpy import ones, zeros, linalg, matrix, array, dot, append, abs, arange, mean
from itertools import islice
from scipy.signal import lfilter, hilbert
from scipy import sqrt, power

#from statlib import stats
# changed routine for the kurtosis (from scipy.stats)
import scipy.stats as ss

 
def window(seq, n=2):
    "Returns a sliding window (of width n) over data from the iterable"
    "   s -> (s0,s1,...s[n-1]), (s1,s2,...,sn), ...                   "
    it = iter(seq)
    result = tuple(islice(it, n))
    if len(result) == n:
        yield result    
    for elem in it:
        result = result[1:] + (elem,)
        yield result

def sw_kurtosis1(x,n):
    """
    Returns kurtosis calculated over data set x using sliding windows of length
    n points.  Length of returned array is len(x)-n+1.
    """
    npts=len(x)
    xs=np.empty(npts-n+1,dtype=float)
    xs[:]=0.
    for i in xrange(n,npts - n ):
        xs[i]=ss.kurtosis(x[i:(i+n)])
    return xs

def sw_kurtosis2(x,n):
    """
    Returns kurtosis calculated over data set x using sliding windows of length
    n points.  Length of returned array is len(x)-n+1.
    This version uses filter.window() to slice the data array and is faster
    than sw_kurtosis1.
    """
    npts=len(x)
    windows=window(x,n)
    xs=np.empty(npts-n+1,dtype=float)
    xs[:]=0.
    k_array=np.empty((npts-n+1,n),dtype=float)
    i=0
    for w in windows:
        k_array[i,0:n]=w
        i=i+1
    xs=ss.kurtosis(k_array,axis=1)
    return xs

def rec_kurtosis_old(x,C):
    npts=len(x)
    varx=np.std(x)
    mean_value=0
    var_value=0
    kurt_value=0
    xs=np.empty(npts)
    for i in xrange(npts):
      mean_value = C*mean_value + (1-C)*x[i]
      var_value=C*var_value+(1-C)*(x[i]-mean_value)**2
      if var_value>varx: 
        kurt_value=C*kurt_value+(1-C)*(x[i]-mean_value)**4/var_value**2
      else : 
        kurt_value=C*kurt_value+(1-C)*(x[i]-mean_value)**4/varx**2
      xs[i]=kurt_value-3

    return xs

def rec_kurtosis(x,C1):
    """
    Recursive Kurtosis calculated using Chassande-Mottin (2002)
    """
    npts = len(x)
    kappa4 = np.empty(npts,dtype=float)

    a1 = 1-C1
    C2 = (1-a1*a1)/2.0
    bias = -3*C1 - 3.0

    mu1_last=0.0
    mu2_last=1.0
    k4_bar_last=0.0
    

    for i in xrange(npts):
        mu1 = a1*mu1_last + C1*x[i]
        dx2 = (x[i]-mu1_last)*(x[i]-mu1_last)
        mu2 = a1*mu2_last + C2*dx2
        dx2 = dx2 / mu2_last
        k4_bar = (1+C1 - 2*C1*dx2)*k4_bar_last + C1 * dx2 * dx2
        kappa4[i] = k4_bar + bias
        mu1_last=mu1
        mu2_last=mu2
        k4_bar_last=k4_bar

    return kappa4

def rec_dx2(x,C1):
    """
    Recursive dx2
    """
    npts = len(x)
    dx2_out = np.empty(npts,dtype=float)

    a1 = 1-C1
    C2 = (1-a1*a1)/2.0

    mu1_last=0
    mu2_last=1
    

    for i in xrange(npts):
        mu1 = a1*mu1_last + C1*x[i]
        dx2 = (x[i]-mu1_last)*(x[i]-mu1_last)
        mu2 = a1*mu2_last + C2*dx2
        dx2 = dx2 / mu2_last
        dx2_out[i] = dx2
        mu1_last=mu1
        mu2_last=mu2

    return dx2_out


def lfilter_zi(b,a):
    """
    Computes the zi state from the filter parameters. See [Gust96].
    
    Based on:
    [Gust96] Fredrik Gustafsson, Determining the initial states in forward-backward 
    filtering, IEEE Transactions on Signal Processing, pp. 988--992, April 1996, 
    Volume 44, Issue 4
    """
    
    n=max(len(a),len(b))
    
    zin = (  np.eye(n-1) - np.hstack( (-a[1:n,np.newaxis], 
                                 np.vstack((np.eye(n-2), zeros(n-2)))))) 
    
    zid=  b[1:n] - a[1:n]*b[0] 
    
    zi_matrix=linalg.inv(zin)*(matrix(zid).transpose())
    zi_return=[]

    #convert the result into a regular array (not a matrix)
    for i in range(len(zi_matrix)):
      zi_return.append(float(zi_matrix[i][0]))
    
    return array(zi_return)

    
def filtfilt(b,a,x):
    """
    What does this function do?  Alberto??
    """
    
    #For now only accepting 1d arrays
    ntaps=max(len(a),len(b))
    edge=ntaps*3
        
    if x.ndim != 1:
        raise ValueError, "Filiflit is only accepting 1 dimension arrays."

    #x must be bigger than edge
    if x.size < edge:
        raise ValueError, "Input vector needs to be bigger than 3 * max(len(a),len(b)."
        
    if len(a) < ntaps:
	a=np.r_[a,zeros(len(b)-len(a))]

    if len(b) < ntaps:
	b=np.r_[b,zeros(len(a)-len(b))]         
    
    zi=lfilter_zi(b,a)
    
    #Grow the signal to have edges for stabilizing 
    #the filter with inverted replicas of the signal
    s=np.r_[2*x[0]-x[edge:1:-1],x,2*x[-1]-x[-1:-edge:-1]]
    #in the case of one go we only need one of the extrems 
    # both are needed for filtfilt
    
    (y,zf)=lfilter(b,a,s,-1,zi*s[0])

    (y,zf)=lfilter(b,a,np.flipud(y),-1,zi*y[-1])
    
    return np.flipud(y[edge-1:-edge+1])
  
  
    
def envelope(x,N):
	"""
	Determines the envelope of a signal using the Hilbert transform.
	
	Uses scipy.signal.hilbert() to do the Hilbert transform.
	"""

	Hx=hilbert(x,N)
	
	# determine the analytic function as f(t) -iFhi(t)
	xana = x - 1j * Hx

    # determine the envelope
	Henv=sqrt((xana * xana.conjugate()).real)
	return(Henv)



def stalta(x,dt,STAwin,LTAwin):
    """
    Determines the short to long time average ratio of a timeseries. 
    
    dt = sampling interval in seconds. 
    STAwin = short time average window in seconds.
    LTAwin = long time average in seconds.
    """

    # find number of samples in averaging windows
    nSTA=int(STAwin/dt)
    nLTA=int(LTAwin/dt)
    
    xabs=abs(x)

    ratio=[]
    for i in range(len(x) - nLTA - nSTA + 1):
        LTA=mean(xabs[i:(i+nLTA)])
        STA=mean(xabs[(i+nLTA):(i+nLTA+nSTA)])
        ratio.append(STA/LTA)
        
    return(ratio)

def allens_stalta(x,dt,STAwin,LTAwin):
    """
    Determines the short to long time average ratio of a timeseries. 
    using Allen (1978) caracteristic function (see report)
    dt = sampling interval in seconds. 
    STAwin = short time average window in seconds.
    LTAwin = long time average in seconds.
    """
    dx=np.gradient(x,dt)
    E=x*x+dx*dx

    # find number of samples in averaging windows
    nSTA=int(STAwin/dt)
    nLTA=int(LTAwin/dt)
    
    Eabs=abs(E)

    ratio=[]
    for i in range(len(E) - nLTA - nSTA + 1):
        LTA=mean(Eabs[i:(i+nLTA)])
        STA=mean(Eabs[(i+nLTA):(i+nLTA+nSTA)])
        ratio.append(STA/LTA)
        
    return(ratio)


def skewness(x,dt,LENwin):
    """
    Determines the absolute value of the skewness of a timeseries. 
    
    dt = sampling interval in seconds. 
    LENwin =  time window (in secs) over which the skewness is determined 
    """

    # find number of samples in averaging windows
    nLEN=int(LENwin/dt)+1
    #xabs=abs(x)

    absskew=[]
    for i in range(len(x) - nLEN +1 ):
        a=ss.skew(x[i:(i+nLEN)])
        absskew.append(abs(a))
    return(absskew)

def kurto(x,dt,LENwin):
    """
    Determines the kurtosis of a timeseries. 
    
    dt = sampling interval in seconds. 
    LENwin =  time window (in secs) over which the kurtosis is determined 
    """

    # find number of samples in averaging windows
    nLEN=int(LENwin/dt)+1
    #xabs=abs(x)

   
    kurtos=[]
    for i in range(len(x) - nLEN +1):
        #a=stats.lkurtosis(x[i:(i+nLEN)])
        # changed routine for kurtosis
        a=ss.kurtosis(x[i:(i+nLEN)],fisher=False)
        kurtos.append(a)
    return(kurtos)


def kurto_improved(x,dt,LENwin):
    """
    Determines the kurtosis of a timeseries. 
    
    dt = sampling interval in seconds. 
    LENwin =  time window (in secs) over which the skewness is determined 

    described by kuperkoch et al. 2010: calculate the kurtosis recursively
    Results are not satisfying, but one may want to improve this... it might
    save some time!
   
    """

    # find number of samples in averaging windows
    nLEN=int(LENwin/dt)+1
    #xabs=abs(x)
    kurtos=[]
    first_window=ss.kurtosis(x[0:(0+nLEN)],fisher=False)
    kurtos.append(first_window)
    i=1
    while i<(len(x) - nLEN +1):
      new_value=kurtos[i-1]-(x[i-1])**4+(x[i-1+nLEN])**4
      kurtos.append(new_value)
      i+=1
    return(kurtos)

def variance(x,dt,LENwin):
    """
    Determines the variance of a timeseries. 
    
    dt = sampling interval in seconds. 
    LENwin =  time window (in secs) over which the variance is determined 
    """

    # find number of samples in averaging windows
    nLEN=int(LENwin/dt)+1
    #xabs=abs(x)

    vari=[]
    for i in range(len(x) - nLEN +1):
        a=np.var(x[i:(i+nLEN)])
        vari.append(a)
    return(vari)

def gradient(x,dt,LENwin):
    """
    Determines the gradient of a timeseries. 
    
    dt = sampling interval in seconds. 
    LENwin =  time window (in secs) over which the gradient is determined 
    """

    # find number of samples in averaging windows
    nLEN=int(LENwin/dt)+1
    #xabs=abs(x)
    grad=[]
    for i in range(len(x) - nLEN +1):
        a=np.gradient(x[i:(i+nLEN)],dt)
        b=np.average(a)
        grad.append(b)
    return(grad)


def variance_under_rotation(x,y,dt,LENwin):
    """
    Determines the variance under rotation of a timeseries. 
    
    dt = sampling interval in seconds. 
    LENwin =  time window (in secs) over which the variance under rotation is determined 
    """

    # find number of samples in averaging windows
    nLEN=int(LENwin/dt)+1
    #xabs=abs(x)

    var=[]
    for i in range(len(x) - nLEN +1):
        a=varrot(x[i:(i+nLEN)],y[i:(i+nLEN)])
        var.append(a)
    return(var)

def varrot(x,y):
  """
  for each angle, calculates the mean of the projection over the whole array
  then calculates the mean over each angle.

  Then, for each angle, calculates the variance, and then the mean over each angle
  """
  av=[]
  for theta in range(0,180,10):
    av.append(np.average(projection(x,y,theta)))
  varrot_mean=np.average(av)

  val_av=[]
  for theta in range(0,180,10):
    val_av.append(np.average(projection(x,y,theta)))
  val=np.average(val_av)
  return (val)
    
def projection(x,y,theta):
  """
   return an array
  """
  t=radians(theta)
  a=np.array(x)
  b=np.array(y)
  proj=a*cos(t)+b*sin(t)
  return (proj)

  
  



def smooth(x,window_len=11,window='hanning'):
    """smooth the data using a window with requested size.
    
    This method is based on the convolution of a scaled window with the signal.
    The signal is prepared by introducing reflected copies of the signal 
    (with the window size) in both ends so that transient parts are minimized
    in the begining and end part of the output signal.
    
    input:
        x: the input signal 
        window_len: the dimension of the smoothing window; should be an odd integer
        window: the type of window from 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'
            flat window will produce a moving average smoothing.

    output:
        the smoothed signal
        
    example:

    t=linspace(-2,2,0.1)
    x=sin(t)+randn(len(t))*0.1
    y=smooth(x)
    
    see also: 
    
    numpy.hanning, numpy.hamming, numpy.bartlett, numpy.blackman, numpy.convolve
    scipy.signal.lfilter
 
    TODO: the window parameter could be the window itself if an array instead of a string
    NOTE: length(output) != length(input), to correct this: return y[(window_len/2-1):-(window_len/2)] instead of just y.
    """

    if x.ndim != 1:
        raise ValueError, "smooth only accepts 1 dimension arrays."

    if x.size < window_len:
        raise ValueError, "Input vector needs to be bigger than window size."


    if window_len<3:
        return x


    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError, "Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'"


    s=np.r_[x[window_len-1:0:-1],x,x[-1:-window_len:-1]]
    #print(len(s))
    if window == 'flat': #moving average
        w=np.ones(window_len,'d')
    else:
        w=eval('np.'+window+'(window_len)')

    y=np.convolve(w/w.sum(),s,mode='valid')
    #return y
    #return y[(window_len/2-1):-(window_len/2)]
    return y[(window_len/2):-(window_len/2)]


