#!/usr/bin/env python
# encoding: utf-8

import os, sys, optparse, glob
import numpy as np
import scipy.signal as si
import logging
import matplotlib.pyplot as plt
from obspy.core import read, utcdatetime, trace
from OP_waveforms import *
from obspy.signal import *
import logging
import cPickle
from locations_trigger import read_locs_from_file
from correlation import BinaryFile
logging.basicConfig(level=logging.INFO, format='%(levelname)s : %(asctime)s : %(message)s')

# ***************************************************************************************
# FUNCTIONS
eps = np.finfo(float).eps

def nextpow2(n):
    m_f = np.log2(n)
    m_i = np.ceil(m_f)
    return 2**m_i


def get_h_parameters(NFIR, fcut):
    """NFIR : length of FIR filter
    fcut: fraction of Nyquist for filter"""
    h = si.firwin(NFIR+1,fcut) * np.exp(2*1j*np.pi*np.arange(NFIR+1)*0.125)
    n = np.arange(2,NFIR+2)
    g = h[(1-n)%NFIR]*(-1)**(1-n)
    NFIR = np.fix((3./2.*NFIR))
    h1 = si.firwin(NFIR+1,2./3*fcut)*np.exp(2j*np.pi*np.arange(NFIR+1)*0.25/3.)
    h2 = h1*np.exp(2j*np.pi*np.arange(NFIR+1)/6.)
    h3 = h1*np.exp(2j*np.pi*np.arange(NFIR+1)/3.)  
    return (h, g, h1, h2, h3)


def plot_kurtogram(Kwav, freq_w, nlevel, Level_w, Fs, fi, index, opt1=None, opt2=None):
  I=index[0]

  fig = plt.figure()
  fig.set_facecolor('white')
  plt.imshow(Kwav,aspect='auto',extent=(0,freq_w[-1],range(2*nlevel)[-1],range(2*nlevel)[0]),interpolation='none')
  xx=np.arange(0,int(freq_w[len(freq_w)-1]),step=5)
  plt.xticks(xx)
  plt.yticks(range(2*nlevel),np.round(Level_w*10)/10)
  #plt.plot(Fs*fi,I,'yo')
  plt.xlabel("Frequency (Hz)")
  plt.ylabel("Level k")
  if opt2==1:
    plt.title("Level %.1f, Bw=%.2f Hz, fc=%.2f Hz"%(np.round(10*Level_w[I])/10,Fs*2**(-(Level_w[I]+1)),Fs*fi))
  else:
    plt.title("Level %.1f, Bw=%.2f Hz, fc=%.2f Hz"%(np.round(10*Level_w[I])/10,Fs*2**(-(Level_w[I]+1)),Fs*fi))
  plt.colorbar()
  #plt.show()


def getBandwidthAndFrequency(nlevel, Fs, level_w, freq_w, level_index, freq_index):

  #f1 = freq_w[freq_index]
  l1 = level_w[level_index]
  fi = (freq_index)/3./2**(nlevel+1)
  fi += 2.**(-2-l1)
  bw = Fs * 2 **-(l1) /2
  fc = Fs * fi

  return bw, fc, fi, l1


def get_GridMax(grid):
    index = np.argmax(grid)
    M = np.amax(grid)
    index = np.unravel_index(index,grid.shape)
    return M, index


def Fast_Kurtogram(x, nlevel,verbose=False, Fs=1, NFIR=16, fcut=0.4, opt1=None, opt2=None):
    # Fast_Kurtogram(x,nlevel,Fs)
    # Computes the fast kurtogram of signal x up to level 'nlevel'
    # Maximum number of decomposition levels is log2(length(x)), but it is 
    # recommended to stay by a factor 1/8 below this.
    # Fs = sampling frequency of signal x (default is Fs = 1)
    # opt1 = 1: classical kurtosis based on 4th order statistics
    # opt1 = 2: robust kurtosis based on 2nd order statistics of the envelope
    # (if there is any difference in the kurtogram between the two measures, this is
    # due to the presence of impulsive additive noise)
    # opt2 = 1: the kurtogram is computed via a fast decimated filterbank tree
    # opt2 = 2: the kurtogram is computed via the short-time Fourier transform
    # (option 1 is faster and has more flexibility than option 2 in the design of the
    # analysis filter: a short filter in option 1 gives virtually the same results as option 2)
    # 
    # -------------------
    # J. Antoni : 02/2005
    # Translation to Python: T. Lecocq 02/2012
    # -------------------

    N = len(x)
    N2 = np.log2(N) - 7
    if nlevel > N2:
       logging.error('Please enter a smaller number of decomposition levels')
       #sys.exit()

    if opt2 is None:
        #~ opt2 = int(raw_input('Choose the kurtosis measure (classic = 1  robust = 2): '))
        opt2 = 1
    if opt1 is None:
        #~ opt1  = int(raw_input('Choose the algorithm (filterbank = 1  stft-based = 2): '))
        opt1  = 1
    # Fast computation of the kurtogram
    ####################################
    
    if opt1 == 1:
        # 1) Filterbank-based kurtogram
        ############################
        # Analytic generating filters

        h, g, h1, h2, h3 = get_h_parameters(NFIR, fcut)

        if opt2 == 1:
            Kwav = K_wpQ(x,h,g,h1,h2,h3,nlevel,verbose,'kurt2')	# kurtosis of the complex envelope
        else:
            Kwav = K_wpQ(x,h,g,h1,h2,h3,nlevel,verbose,'kurt1') # variance of the envelope magnitude
       
        # keep positive values only!
        Kwav[Kwav <= 0] = 0
        Level_w = np.arange(1,nlevel+1)
        Level_w = np.array([Level_w, Level_w + np.log2(3.)-1])
        Level_w = sorted(Level_w.ravel())
        Level_w = np.append(0,Level_w[0:2*nlevel-1])
        freq_w = Fs*(np.arange(0,3*2.0**nlevel)/(3*2**(nlevel+1)) + 1.0/(3.*2.**(2+nlevel)))

        M, index=get_GridMax(Kwav)
        level_index=index[0]
        freq_index=index[1]
        bw, fc, fi, l1 = getBandwidthAndFrequency(nlevel, Fs, Level_w, freq_w, level_index, freq_index)

        if verbose:
          plot_kurtogram(Kwav, freq_w, nlevel, Level_w, Fs, fi, index, opt1, opt2)

    else:
        logging.error('stft-based is not implemented')

    # Signal filtering !
    c = [];
    test = 1
    lev = l1
    while test == 1:
        test=0
        if opt2 == 1:
          c,s,threshold,Bw,fc = Find_wav_kurt(x,h,g,h1,h2,h3,nlevel,lev,fi,Fs=Fs,verbose=verbose)
          kx = kurt(c,'kurt2')
        else:
          c,s,threshold,Bw,fc = Find_wav_kurt(x,h,g,h1,h2,h3,nlevel,lev,fi,Fs=Fs,verbose=verbose)

    # Apply appropriate filter and compare with previous results (i.e. 4-10 Hz filter) 
    f_lower=Fs*np.round((fc-Bw/2.)*10**3)/10**3
    f_upper=Fs*np.round((fc+Bw/2.)*10**3)/10**3

    return Kwav, Level_w, freq_w, c, f_lower, f_upper

   

def K_wpQ(x,h,g,h1,h2,h3,nlevel,verbose,opt,level=0):
    # K = K_wpQ(x,h,g,h1,h2,h3,nlevel)
    # Calculates the kurtosis K of the complete quinte wavelet packet transform w of signal x, 
    # up to nlevel, using the lowpass and highpass filters h and g, respectively.
    # The WP coefficients are sorted according to the frequency decomposition.
    # This version handles both real and analytical filters, but does not yield WP coefficients
    # suitable for signal synthesis.
    #
    # -----------------------
    # J. Antoni : 12/2004 
    # Translation: T. Lecocq 02/2012
    # -----------------------
    L = np.floor(np.log2(len(x)))
    if level == 0:
        if nlevel >= L:
            logging.error('nlevel must be smaller')
        level=nlevel
    x = x.ravel()
    KD, KQ = K_wpQ_local(x,h,g,h1,h2,h3,nlevel,verbose,opt,level)
    K = np.zeros((2*nlevel,3*2**nlevel))
    #print "******************************************************"
    #print KD.shape, KQ.shape, K.shape
    #~ K = KD
    K[0,:]=KD[0,:]
    for i in range(1,nlevel):
        K[2*i-1,:] = KD[i,:]
        K[2*i,:] = KQ[i-1,:]
   
    K[2*nlevel-1,:] = KD[nlevel,:]
    #print "K Final Shape", K.shape
    return K


def K_wpQ_local(x,h,g,h1,h2,h3,nlevel,verbose,opt,level):
    #if verbose:
      #print "LEVEL", level
    a,d = DBFB(x,h,g)

    N = len(a)
    d = d*np.power(-1.,np.arange(1,N+1)) # indices pairs (commencent à 0) multipliés par -1
    K1 = kurt(a[len(h)-1:],opt)
    K2 = kurt(d[len(g)-1:],opt)

    if level > 2:
        a1,a2,a3 = TBFB(a,h1,h2,h3)
        d1,d2,d3 = TBFB(d,h1,h2,h3)
        Ka1 = kurt(a1[len(h)-1:],opt)
        Ka2 = kurt(a2[len(h)-1:],opt)
        Ka3 = kurt(a3[len(h)-1:],opt)
        Kd1 = kurt(d1[len(h)-1:],opt)
        Kd2 = kurt(d2[len(h)-1:],opt)
        Kd3 = kurt(d3[len(h)-1:],opt)
    else:
        Ka1 = 0
        Ka2 = 0
        Ka3 = 0
        Kd1 = 0
        Kd2 = 0
        Kd3 = 0
    
    if level==1:
        K =np.array([K1*np.ones(3),K2*np.ones(3)]).flatten()
        KQ = np.array([Ka1,Ka2,Ka3,Kd1,Kd2,Kd3])
        #if verbose:
          #print "level = 1"
          #print 'K.shape',K.shape
          #print 'KQ.shape',KQ.shape
    if level > 1:
        #if verbose:
          #print "entering rec with level %i"%(level-1)
          #print "doing A"
        Ka,KaQ = K_wpQ_local(a,h,g,h1,h2,h3,nlevel,verbose,opt,level-1)
        #if verbose:
          #print "doing D"
        Kd,KdQ = K_wpQ_local(d,h,g,h1,h2,h3,nlevel,verbose,opt,level-1)
        #if verbose:
          #print "out of rec level %i" % (level -1)
        #print Ka.shape, Kd.shape
        K1 = K1*np.ones(np.max(Ka.shape))
        K2 = K2*np.ones(np.max(Kd.shape))
        K12 = np.append(K1,K2)
        Kad = np.hstack((Ka, Kd))
        #print ">", K12.shape, Kad.shape
        K = np.vstack((K12,Kad))

        Long = 2./6*np.max(KaQ.shape)
        Ka1 = Ka1*np.ones(Long)
        Ka2 = Ka2*np.ones(Long)
        Ka3 = Ka3*np.ones(Long)
        Kd1 = Kd1*np.ones(Long)
        Kd2 = Kd2*np.ones(Long)
        Kd3 = Kd3*np.ones(Long)
        tmp = np.hstack((KaQ,KdQ))
        #print "HEEEERE"
        #print tmp.shape
        KQ = np.concatenate((Ka1,Ka2,Ka3,Kd1,Kd2,Kd3))
        KQ = np.vstack((KQ, tmp))

        #~ if tmp.shape[0] != KQ.shape[0]:
            #~ tmp = tmp.T
        #~ for i in range(tmp.shape[0]):
            #~ KQ = np.vstack((KQ,tmp[i]))
        
        #print "4", K.shape, KQ.shape       
 
    if level == nlevel:
        K1 = kurt(x,opt)
        K = np.vstack((K1*np.ones(np.max(K.shape)), K))
        #print "K shape", K.shape

        a1,a2,a3 = TBFB(x,h1,h2,h3)
        Ka1 = kurt(a1[len(h)-1:],opt)
        Ka2 = kurt(a2[len(h)-1:],opt)
        Ka3 = kurt(a3[len(h)-1:],opt)
        Long = 1./3*np.max(KQ.shape)
        Ka1 = Ka1*np.ones(Long)
        Ka2 = Ka2*np.ones(Long)
        Ka3 = Ka3*np.ones(Long)
        #print KQ.shape
        tmp = np.array(KQ[0:-2])
        #print "level==nlevel"
        
        KQ = np.concatenate((Ka1,Ka2,Ka3))
        KQ = np.vstack((KQ,tmp))
 
    #print "i'm leaving level=%i and K.shape="%level,K.shape, "and KQ.shape=",KQ.shape
    return K, KQ


def kurt(x, opt):
    if opt=='kurt2':
        if np.all(x==0):
            K=0
            E=0
            return K
        x = x - np.mean(x)
        E = np.mean(np.abs(x)**2)
        if E < eps:
            K=0
            return K

        K = np.mean(np.abs(x)**4)/E**2

        if np.all(np.isreal(x)):
            K = K - 3
        else:
            K = K - 2

    if opt=='kurt1':
        if np.all(x==0):
            K=0
            E=0
            return K
        x = x - np.mean(x)
        E = np.mean(np.abs(x))
        if E < eps:
            K=0
            return K

        K = np.mean(np.abs(x)**2)/E**2

        if np.all(np.isreal(x)):
            K=K-1.57
        else:
            K=K-1.27

    return K


def DBFB(x,h,g):
    # Double-band filter-bank.
    #   [a,d] = DBFB(x,h,g) computes the approximation
    #   coefficients vector a and detail coefficients vector d,
    #   obtained by passing signal x though a two-band analysis filter-bank.
    #   h is the decomposition low-pass filter and
    #   g is the decomposition high-pass filter.
    
    N = len(x)
    La = len(h)
    Ld = len(g)

    # lowpass filter
    a = si.lfilter(h,1,x)
    a = a[1::2]
    a = a.ravel()

    # highpass filter
    d = si.lfilter(g,1,x)
    d = d[1::2]
    d = d.ravel()
    return (a,d)


def TBFB(x,h1,h2,h3):
    # Trible-band filter-bank.
    #   [a1,a2,a3] = TBFB(x,h1,h2,h3) 
    
    N = len(x)
    La1 = len(h1)
    La2 = len(h2)
    La3 = len(h3)

    # lowpass filter
    a1 = si.lfilter(h1,1,x)
    a1 = a1[2::3]
    a1 = a1.ravel()

    # passband filter
    a2 = si.lfilter(h2,1,x)
    a2 = a2[2::3]
    a2 = a2.ravel()

    # highpass filter
    a3 = si.lfilter(h3,1,x)
    a3 = a3[2::3]
    a3 = a3.ravel()
    return (a1,a2,a3)


def Find_wav_kurt(x,h,g,h1,h2,h3,nlevel,Sc,Fr,Fs=1,verbose=False):
    # [c,Bw,fc,i] = Find_wav_kurt(x,h,g,h1,h2,h3,nlevel,Sc,Fr,opt2)
    # Sc = -log2(Bw)-1 with Bw the bandwidth of the filter
    # Fr is in [0 .5]
    #
    # -------------------
    # J. Antoni : 12/2004
    # -------------------
    level = np.fix((Sc))+ ((Sc%1) >= 0.5) * (np.log2(3)-1)
    Bw = 2**(-level-1)
    freq_w = np.arange(0,2**level) / 2**(level+1) + Bw/2.
    J = np.argmin(np.abs(freq_w-Fr))
    fc = freq_w[J]
    i = int(np.round(fc/Bw-1./2))
    if level % 1 == 0:
        acoeff = binary(i, level)
        bcoeff = []
        temp_level = level
    else:
        i2 = np.fix((i/3.))
        temp_level = np.fix((level))-1
        acoeff = binary(i2,temp_level)
        bcoeff = i-i2*3
    acoeff = acoeff[::-1]
    c = K_wpQ_filt(x,h,g,h1,h2,h3,acoeff,bcoeff,temp_level)

    t = np.arange(len(x))/float(Fs)
    tc = np.linspace(t[0],t[-1],len(c))
    s=np.real(c*np.exp(2j*np.pi*fc*Fs*tc))

    sig = np.median(np.abs(c))/np.sqrt(np.pi/2.)
    threshold = sig*raylinv(np.array([.999,]),np.array([1,]))

    if verbose:
      plot_envelope(x, Fs, c, fc, level)

    return c,s,threshold,Bw,fc 


def getFTSquaredEnvelope(c):
  nfft = int(nextpow2(len(c)))
  env = np.abs(c)**2
  S = np.abs(np.fft.fft( (env.ravel()-np.mean(env))*np.hanning(len(env))/len(env),nfft))
  return S


def plot_envelope(x, Fs, c, fc, level):
  #spec = int(raw_input('	Do you want to see the envelope spectrum (yes = 1 ; no = 0): '))
  spec=0

  fig = plt.figure()
  fig.set_facecolor('white')
  t = np.arange(len(x))/Fs
  tc = np.linspace(t[0],t[-1],len(c))
  plt.subplot(2+spec,1,1)
  plt.plot(t,x/np.max(x),'k')
  plt.plot(tc,np.abs(c)/np.max(np.abs(c)),'r')
  plt.title('Original Signal (4-10 Hz)')

  plt.subplot(2+spec,1,2)
  plt.plot(tc,np.abs(c),'k')
  #plt.plot(tc,threshold*np.ones(len(c)),'--r')
  plt.title("Envelope of the filtered signal, Bw=Fs/2^%.1f, fc=%.2f Hz"%(np.round(level*10)/10,Fs*fc))
  #plt.title("Envelope of the filtered signal, Bw=Fs/2^%.1f, fc=%.2f Hz, Kurt=%.1f, alpha=.1%%"%(np.round((level+1)*10)/10,Fs*fc,np.round(np.abs(10*kx))/10),size='small')
  plt.xlabel('time [s]')
  if spec == 1:
    #print nextpow2(len(c))
    nfft = int(nextpow2(len(c)))
    env = np.abs(c)**2
    S = np.abs(np.fft.fft(env.ravel()-np.mean(env)*np.hanning(len(env))/len(env),nfft))
    f = np.linspace(0, 0.5*Fs/2**level,nfft/2)
    plt.subplot(313)
    plt.plot(f,S[:nfft/2],'k')
    plt.title('Fourier transform magnitude of the squared envelope')
    plt.xlabel('frequency [Hz]')
    #plt.show()


def binary(i,k):
    # return the coefficients of the binary expansion of i:
    # i = a(1)*2^(k-1) + a(2)*2^(k-2) + ... + a(k)

    if i>=2**k:
        logging.error('i must be such that i < 2^k !!')
    
    a = np.zeros(k)
    temp = i
    for l in np.arange(k-1,-1,-1):
        a[k-l-1] = np.fix(temp/2**l)
        temp = temp - a[k-l-1]*2**l
    return a


def K_wpQ_filt(x,h,g,h1,h2,h3,acoeff,bcoeff,level=0):
    # c = K_wpQ_filt(x,h,g,h1,h2,h3,acoeff,bcoeff,level)
    # Calculates the kurtosis K of the complete quinte wavelet packet transform w of signal x,
    # up to nlevel, using the lowpass and highpass filters h and g, respectively.
    # The WP coefficients are sorted according to the frequency decomposition.
    # This version handles both real and analytical filters, but does not yield WP coefficients
    # suitable for signal synthesis.
    #
    # -----------------------
    # J. Antoni : 12/2004
    # -----------------------
    nlevel = len(acoeff)
    L = np.floor(np.log2(len(x)))
    if level==0:
        if nlevel >= L:
            logging.error('nlevel must be smaller !!')
        level = nlevel
    x = x.ravel()
    if nlevel == 0:
        if bcoeff==[]:
            c = x
        else:
            c1, c2, c3 = TBFB(x,h1,h2,h3)
            if bcoeff == 0:
                c = c1[len(h1)-1:]
            elif bcoeff == 1:
                c = c2[len(h2)-1:]
            elif bcoeff == 2:
                c = c3[len(h3)-1:]
    else:
        c = K_wpQ_filt_local(x,h,g,h1,h2,h3,acoeff,bcoeff,level)
    return c


def K_wpQ_filt_local(x,h,g,h1,h2,h3,acoeff,bcoeff,level):
    #print level, x[:9]
    a,d = DBFB(x,h,g)         # perform one analysis level into the analysis tree
    N = len(a)
    d = d*np.power(-1.,np.arange(1,N+1))
    level = int(level)
    if level == 1:
        if bcoeff==[]:
          if acoeff[level-1] == 0:
             c = a[len(h)-1:]
          else:
             c = d[len(g)-1:]
        else:
            if acoeff[level-1] == 0:
                c1,c2,c3 = TBFB(a,h1,h2,h3)
            else:
                c1,c2,c3 = TBFB(d,h1,h2,h3)
            if bcoeff == 0:
                c = c1[len(h1)-1:]
            elif bcoeff == 1:
                c = c2[len(h2)-1:]
            elif bcoeff == 2:
                c = c3[len(h3)-1:]
    if level > 1:
        if acoeff[level-1] == 0:
            c = K_wpQ_filt_local(a,h,g,h1,h2,h3,acoeff,bcoeff,level-1)
        else:
            c = K_wpQ_filt_local(d,h,g,h1,h2,h3,acoeff,bcoeff,level-1)
    #print 'kurt', kurt(c,'kurt2')
    return c


def raylinv(p,b):
    #RAYLINV  Inverse of the Rayleigh cumulative distribution function (cdf).
    #   X = RAYLINV(P,B) returns the Rayleigh cumulative distribution
    #   function with parameter B at the probabilities in P.

    #~ if nargin <  1:
        #~ logging.error('Requires at least one input argument.')

    # Initialize x to zero.
    x = np.zeros(len(p))
    # Return NaN if the arguments are outside their respective limits.
    k = np.where(((b <= 0)| (p < 0)| (p > 1)))[0]
    
    if len(k) != 0:
        tmp  = np.NaN
        x[k1] = tmp(len(k))

    # Put in the correct values when P is 1.
    k = np.where(p == 1)[0]
    #print k
    if len(k)!=0:
        tmp  = Inf
        x[k] = tmp(len(k))

    k = np.where(((b > 0) & (p > 0) & (p < 1)))[0]
    #print k
    
    if len(k)!=0:
        pk = p[k]
        bk = b[k]
        #print pk, bk
        x[k] = np.sqrt((-2*bk ** 2) * np.log(1 - pk))
    return x
# --------------------------------------------------------------------------
def plot_trace(x,xfilt,kurtx,tr,c,info,f_lower,f_upper,snr,snr_ref,snr_kurt,kmax,kmax_ref,tstack):
  dt=info['dt']
  sta=info['station']  

  fig=plt.figure()
  fig.set_facecolor('white')

  ax1=fig.add_subplot(311, title="Signal filtered between 4 - 10 Hz")
  ax1.plot(x/np.max(np.abs(x)),'k')
  ax1.plot(kurtx/np.max(np.abs(kurtx)),'y')
  ax1.text(10*dt,-0.8,"SNR: %.1f \nKmax: %.2f"%(snr_ref,kmax_ref))

  ax2=fig.add_subplot(312, title="Signal filtered between %.1f - %.1f Hz"%(f_lower, f_upper))
  ax2.plot(xfilt/np.max(np.abs(xfilt)),'k')
  ax2.plot(tr/np.max(np.abs(tr)),'r')
  ax2.text(10*dt,-0.8,"SNR: %.1f \nKmax: %.2f"%(snr,kmax))

  ax3=fig.add_subplot(313, title="Kurtosis comparison")
  ax3.plot(tr,'r',label="new")
  ax3.plot(kurtx,'y',label="initial")
  handles, labels = ax3.get_legend_handles_labels()
  ax3.legend(handles, labels)

  plt.setp(ax1.get_xticklabels(), visible=False)
  plt.setp(ax2.get_xticklabels(), visible=False)
  fig.suptitle("%s - %s"%(sta,tstack))

# --------------------------------------------------------------------------
def write_file(info, tstart, tend, tr):
  print "Writing a new kurtosis file..."
  t_orig=info['tdeb']
  dt=info['dt']
  ind1=int((tstart-t_orig)/dt)
  ind2=int((tend-t_orig)/dt)
  if ind1 < 0:
    tr=tr[-ind1:]
    ind1=0
  if ind2 > len(info['new_kurt']):
    ind2=len(info['new_kurt'])-1

  if len(info['new_kurt'][ind1:ind2+1])==len(tr):
    info['new_kurt'][ind1:ind2+1]=tr
  return info
# --------------------------------------------------------------------------
def waveval(xall,tstart,tend,dt,tdeb):
  tstart=tstart-tdeb
  tend=tend-tdeb
  istart=int(round(tstart*1./dt))
  iend=int(round(tend*1./dt))
  val=xall[istart-1:iend]
  return val
# --------------------------------------------------------------------------
def kurto(origin_time, info, opdict):
  verbose=opdict['verbose']
  kwin=opdict['kwin']

  start_time=origin_time-5.0
  end_time=origin_time+20.0
  dt=info['dt']
  
  # Trace
  x=waveval(info['data_ini'],start_time,end_time,dt,info['tdeb_data'])
  if not x.any() and x.all():
    return info
 
  # Initial kurtosis (trace filtered between 4-10Hz)
  kurtx=waveval(info['kurt_ini'],start_time,end_time,dt,info['tdeb_kurt'])
  kurtx=smooth(kurtx)

  N=len(x)
  N2=np.log2(N)-7
  nlevel=int(np.fix(N2))
  
  snr_ref=np.max(np.abs(x))/np.mean(np.abs(x))
  snr_kurt_ref=np.max(np.abs(kurtx))/np.mean(np.abs(kurtx))
  kmax_ref=np.max(kurtx) # maximum of the kurtosis

  # Compute the kurtogram and keep best frequencies
  Kwav, Level_w, freq_w, c, f_lower, f_upper = Fast_Kurtogram(x, nlevel,verbose,Fs=1/dt,opt2=1)

  # Comparison of the kurtosis computed in the new frequency band and the old one (criterion : snr, kmax)
  # 1. Read the initial data
  wf=Waveform()
  wf.read_from_file(info['data_file'], starttime=start_time-kwin, endtime=end_time+kwin)

  nbpts=int(kwin*1./dt)

  # 2. Filter the trace with kurtogram frequencies
  wf.bp_filter(f_lower,f_upper)
  x_filt=wf.values
  x_filt=x_filt[nbpts:-nbpts]

  # 3. Compute the kurtosis
  wf.process_kurtosis(kwin,recursive=opdict['krec'])
  new_kurtx=wf.values
  new_kurtx=new_kurtx[nbpts+1:-nbpts-1]

  snr=np.max(np.abs(x_filt))/np.mean(np.abs(x_filt))
  snr_kurt=np.max(np.abs(new_kurtx))/np.mean(np.abs(new_kurtx))
  kmax=np.max(new_kurtx)

  if snr > snr_ref and kmax >= kmax_ref:
    info['filter'].append((round(f_lower*100)/100,round(f_upper*100)/100))
    if info.has_key('new_kurt_file'):
      info=write_file(info, start_time, end_time, new_kurtx)
  else:
    info['filter'].append((0,50))

  if verbose and snr > 3:
    print "snr:", snr, " ; snr_ref:", snr_ref
    print "snr new kurtosis:", snr_kurt, " ; snr kurtosis reference:", snr_kurt_ref
    print "kurtosis max, kurt_ref :", kmax, kmax_ref
    plot_trace(x,x_filt,kurtx,new_kurtx,c,info,f_lower,f_upper,snr, snr_ref, snr_kurt, kmax, kmax_ref,origin_time)
    plt.show()

  return info
# *************************************************************************
def read_kurtogram_frequencies(filename):
  a=BinaryFile(filename)
  freqs=a.read_binary_file()

  for staname in sorted(freqs):
    print "%s %.1f %.1f"%(staname,np.mean(freqs[staname][:,0]),np.mean(freqs[staname][:,1]))
    fig=plt.figure()
    fig.set_facecolor('white')
    plt.hist([freqs[staname][:,0],freqs[staname][:,1]],35,histtype='stepfilled',alpha=.2,color=('b','g'),label=['f_low','f_up'])
    plt.title(staname)
    plt.xlabel('Frequency (Hz)')
    plt.figtext(0.15,0.85,"Lower f = %.1f Hz"%np.mean(freqs[staname][:,0]))
    plt.figtext(0.15,0.8,"Upper f = %.1f Hz"%np.mean(freqs[staname][:,1]))
    plt.show()
# *************************************************************************
def do_kurtogram_setup_and_run(opdict):

  base_path=opdict['base_path']

  # data
  data_dir=os.path.join(base_path,'data',opdict['datadir'])
  data_glob=opdict['dataglob']
  data_files=glob.glob(os.path.join(data_dir,data_glob))
  data_files.sort()

  kurt_glob=opdict['kurtglob']
  kurt_files=glob.glob(os.path.join(data_dir,kurt_glob))
  kurt_files.sort()

  # output directory
  out_dir=os.path.join(base_path,'out',opdict['outdir'])

  # location file
  locdir=os.path.join(out_dir,'loc')
  locfile=os.path.join(locdir,'locations.dat')
  # Read locations
  locs=read_locs_from_file(locfile)
  nb_tot_event=len(locs)

  # create a file containing the best filtering parameters for each event and each station
  kurto_file=os.path.join(out_dir,'kurto')

  tdeb=utcdatetime.UTCDateTime(opdict['starttime'])
  tfin=utcdatetime.UTCDateTime(opdict['endtime'])

  # write filenames in a dictionary
  kurtdata={}
  for filename in kurt_files:
    try:
      wf=Waveform()
      wf.read_from_file(filename)
      sta=wf.station
      kurtdata[sta]=filename
    except UserWarning:
      logging.info('No data around %s for file %s.'%(o_time.isoformat(),filename))

  data={}
  for filename in data_files:
    try:
      wf=Waveform()
      wf.read_from_file(filename)
      sta=wf.station
      data[sta]=filename
    except UserWarning:
      logging.info('No data around %s for file %s.'%(o_time.isoformat(),filename))

  # --------------------------------------------------------------------------
  # Create an empty dictionnary that will contain the filtering parameters
  param={}

  for station in sorted(data):
    wf1=Waveform()
    wf1.read_from_file(data[station],starttime=tdeb,endtime=tfin)

    wf2=Waveform()
    wf2.read_from_file(kurtdata[station],starttime=tdeb,endtime=tfin)
  
    info={}
    info['data_file']=data[station]
    info['station']=station
    info['tdeb_data']=wf1.starttime
    info['tdeb_kurt']=wf2.starttime
    info['kurt_file']=kurtdata[station]
    info['data_ini']=wf1.values
    info['kurt_ini']=wf2.values
    info['dt']=wf1.dt
    info['filter']=[]
 
    logging.info('Processing station %s'%info['station']) 
    
    if opdict['new_kurtfile']:
      new_filename='filt_kurtogram.sac'
      info['new_kurt_file']=os.path.join("%s%s"%(data[station].split(data_glob[1:])[0],new_filename))
      trace_kurt_fin=Waveform()
      trace_kurt_fin.read_from_file(kurt_file)
      info['new_kurt']=trace_kurt_fin.values

    for loc in locs:
      origin_time=loc['o_time']
      if opdict['verbose']:
        print "***************************************************************"
        print logging.info(origin_time)

      if origin_time > tdeb and origin_time < tfin:
        info=kurto(origin_time, info, opdict)
      else:
        break

    info['filter']=np.matrix(info['filter'])
    sta=info['station']
    param[sta]=info['filter']

    if info.has_key('new_kurt_file'):
      trace_kurt_fin.values[:]=info['new_kurt']
      trace_kurt_fin.write_to_file_filled(info['new_kurt_file'],format='MSEED',fill_value=0)

  # Write the dictionnary 'param' in a binary file
  a=BinaryFile(kurto_file)
  a.write_binary_file(param)
  read_kurtogram_frequencies(kurto_file)

if __name__ == '__main__':
  from options import WavelocOptions

  logging.basicConfig(level=logging.INFO, format="%(levelname)s : %(asctime)s : %(message)s")

  wo = WavelocOptions()
  args=wo.p.parse_args()

  wo.set_all_arguments(args)
  wo.verify_kurtogram_options()

  do_kurtogram_setup_and_run(wo.opdict)
