#!/usr/bin/env python
# encoding: utf-8

import numpy as np
import scipy.signal as si
import logging
import matplotlib.pyplot as plt
import sys, time
import scipy.stats as ss

# ***************************************************************************************
# FUNCTIONS
eps = np.finfo(float).eps

def nextpow2(n):
    m_f = np.log2(n)
    m_i = np.ceil(m_f)
    return 2**m_i
    
def Fast_Kurtogram(x, nlevel,options_verbose, Fs=1, opt1=None, opt2=None):
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
        N = 16			
        fc = .4	# a short filter is just good enough!
        h = si.firwin(N+1,fc) * np.exp(2*1j*np.pi*np.arange(N+1)*0.125)
        n = np.arange(2,N+2)
        g = h[(1-n)%N]*(-1)**(1-n)
        N = np.fix((3./2.*N))
        h1 = si.firwin(N+1,2./3*fc)*np.exp(2j*np.pi*np.arange(N+1)*0.25/3.)
        #~ plt.plot(h1)
        #~ plt.show()
        h2 = h1*np.exp(2j*np.pi*np.arange(N+1)/6.)
        h3 = h1*np.exp(2j*np.pi*np.arange(N+1)/3.)
       
        if opt2 == 1:
            Kwav = K_wpQ(x,h,g,h1,h2,h3,nlevel,options_verbose,'kurt2')	# kurtosis of the complex envelope
        else:
            Kwav = K_wpQ(x,h,g,h1,h2,h3,nlevel,options_verbose,'kurt1') # variance of the envelope magnitude
       
        # keep positive values only!
        Kwav[Kwav <= 0] = 0
        Level_w = np.arange(1,nlevel+1)
        Level_w = np.array([Level_w, Level_w + np.log2(3.)-1])
        Level_w = sorted(Level_w.ravel())
        Level_w = np.append(0,Level_w[0:2*nlevel-1])
        freq_w = Fs*(np.arange(0,3*2.0**nlevel)/(3*2**(nlevel+1)) + 1.0/(3.*2.**(2+nlevel)))

        #I,J,M=max_IJ(Kwav)
        index = np.argmax(Kwav)
        index = np.unravel_index(index,Kwav.shape)
        I,J,M=index[0],index[1],np.max(Kwav)
        f1 = freq_w[index[1]]
        #print "index : ",index
        l1 = Level_w[index[0]]
        fi = (index[1])/3./2**(nlevel+1)
        fi += 2.**(-2-l1)

        if options_verbose:
          fig = plt.figure()
          fig.set_facecolor('white')
          #~ plt.subplot(ratio='auto')
          plt.imshow(Kwav,aspect='auto',extent=(0,freq_w[-1],range(2*nlevel)[-1],range(2*nlevel)[0]),interpolation='bilinear')
          xx=np.arange(0,int(freq_w[len(freq_w)-1]),step=5)
          plt.xticks(xx)
          plt.yticks(range(2*nlevel),np.round(Level_w*10)/10)
          #plt.plot(Fs*fi,I,'yo')
          plt.xlabel("Frequency (Hz)")
          plt.ylabel("Level k")
          print freq_w[-1]
          if opt2==1:
            plt.title("Level %.1f, Bw=%.2f Hz, fc=%.2f Hz"%(np.round(10*Level_w[I])/10,Fs*2**(-(Level_w[I]+1)),Fs*fi))
            #plt.title("fb-kurt.2 - Kmax=%.1f, level %.1f, Bw=%.2f Hz, fc=%.2f Hz"%(np.round(10*M)/10,np.round(10*Level_w[I])/10,Fs*2**(-(Level_w[I]+1)),Fs*fi))
          else:
            plt.title("Level %.1f, Bw=%.2f Hz, fc=%.2f Hz"%(np.round(10*Level_w[I])/10,Fs*2**(-(Level_w[I]+1)),Fs*fi))
            #plt.title("fb-kurt.1 - Kmax=%.1f, level %.1f, Bw=%.2f Hz, fc=%.2f Hz"%(np.round(10*M)/10,np.round(10*Level_w[I])/10,Fs*2**(-(Level_w[I]+1)),Fs*fi))
          plt.colorbar()
          fig.savefig('kurtogram.png')
    else:
        logging.error('stft-based is not implemented')
    
    # Signal filtering !
    c = [];
    #~ test = int(raw_input('Do you want to filter out transient signals from the kurtogram (yes = 1 ; no = 0): '))
    test = 1
    #~ fi = fi * Fs
    lev = l1
    while test == 1:
        test=0
        #~ fi = float(raw_input(['	Enter the optimal carrier frequency (btw 0 and ',num2str(Fs/2),') where to filter the signal: ']));
        #~ fi = fi/Fs;
        #~ if opt1 == 1:
            #~ lev = raw_input(['	Enter the optimal level (btw 0 and ',num2str(nlevel),') where to filter the signal: ']);
        if opt2 == 1:
            c,Bw,fc = Find_wav_kurt(x,h,g,h1,h2,h3,nlevel,options_verbose,lev,fi,'kurt2',Fs)
            kx = kurt(c,'kurt2')
        #~ else
            #~ [c,Bw,fc] = Find_wav_kurt(x,h,g,h1,h2,h3,nlevel,lev,fi,'kurt1',Fs);
        #test = int(raw_input('Do you want to keep on filtering out transients (yes = 1 ; no = 0): '))

    # Apply appropriate filter and compare with previous results (i.e. 4-10 Hz filter) 
    f_lower=Fs*np.round((fc-Bw/2.)*10**3)/10**3
    f_upper=Fs*np.round((fc+Bw/2.)*10**3)/10**3
    
    return c, f_lower, f_upper

#def max_IJ(X):
    # [I,J] = max_IJ(X) returns the row and column indices of the maximum in the matrix X.
    # J. Antoni : 07/2004
#    M=np.max(X)
#    I,J=np.where(X==M)
#    I,J=I[0],J[0]
#    return I,J,M

def K_wpQ(x,h,g,h1,h2,h3,nlevel,options_verbose,opt,level=0):
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
    KD, KQ = K_wpQ_local(x,h,g,h1,h2,h3,nlevel,options_verbose,opt,level)
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

def K_wpQ_local(x,h,g,h1,h2,h3,nlevel,options_verbose,opt,level):
    #if options_verbose:
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
        #if options_verbose:
          #print "level = 1"
          #print 'K.shape',K.shape
          #print 'KQ.shape',KQ.shape
    if level > 1:
        #if options_verbose:
          #print "entering rec with level %i"%(level-1)
          #print "doing A"
        Ka,KaQ = K_wpQ_local(a,h,g,h1,h2,h3,nlevel,options_verbose,opt,level-1)
        #if options_verbose:
          #print "doing D"
        Kd,KdQ = K_wpQ_local(d,h,g,h1,h2,h3,nlevel,options_verbose,opt,level-1)
        #if options_verbose:
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

def Find_wav_kurt(x,h,g,h1,h2,h3,nlevel,options_verbose,Sc,Fr,opt,Fs=1):
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

    kx = kurt(c,opt)
    #print "kx", kx

    sig = np.median(np.abs(c))/np.sqrt(np.pi/2.)
    #print sig
    threshold = sig*raylinv(np.array([.999,]),np.array([1,]))
    #print "threshold", threshold
    #spec = int(raw_input('	Do you want to see the envelope spectrum (yes = 1 ; no = 0): '))
    spec=0

    if options_verbose:
      fig = plt.figure()
      fig.set_facecolor('white')
      t = np.arange(len(x))/Fs
      tc = np.linspace(t[0],t[-1],len(c))
      plt.subplot(2+spec,1,1)
      plt.plot(t,x/np.max(x),'k')
      plt.plot(tc,np.abs(c)/np.max(np.abs(c)),'r')
      plt.title('Original Signal')

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
    return [c,Bw,fc]

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

if __name__=='__main__':
  # *************************************************************************
  # MAIN PROGRAM
  import os, sys, optparse, glob
  from obspy.core import read, utcdatetime, trace
  from OP_waveforms import *
  from obspy.signal import *
  import logging
  logging.basicConfig(level=logging.INFO, format='%(levelname)s : %(asctime)s : %(message)s')

  # get path
  base_path=os.getenv('WAVELOC_PATH')

  # Read command line
  p=optparse.OptionParser()
  p.add_option('--outdir', '-o', action='store', help='output subdirectory')
  p.add_option('--datadir', '-d', action='store', help='data subdirectory')
  p.add_option('--data_glob',action='store',help="data glob")
  p.add_option('--kurt_glob',action='store',help="kurtosis glob")
  p.add_option('--stations','-s',action='store',default='channels_HHZ.dat',help='station list (found in $WAVELOC_PATH/lib)')
  p.add_option('--starttime',action='store',help="start time for data e.g. 2010-10-14T00:00:00.0Z")
  p.add_option('--endtime',action='store',help="end time for data e.g. 2010-10-14T10:00:00.0Z")
  p.add_option('--new_file',action='store',help="name of new kurtosis files")
  p.add_option('--refine',action='store_true',help="apply kurtogram twice")
  p.add_option('--verbose','-v',action='store_true',help='print debugging information to stout')
  (options,arguments)=p.parse_args()

  out_dir="%s/out/%s"%(base_path,options.outdir)
  loc_filename="%s/loc/locations.dat"%out_dir
  data_dir="%s/data/%s"%(base_path,options.datadir)
  sta_filename="%s/lib/%s"%(base_path,options.stations)

  with open(loc_filename,'r') as locfile:
    loc_lines=locfile.readlines()
    locfile.close()
  # --------------------------------------------------------------------------
  # Read station file
  with open(sta_filename,'r') as sta:
    sta_lines=sta.readlines()
    sta.close()
  # --------------------------------------------------------------------------

  tdeb=utcdatetime.UTCDateTime(options.starttime)
  tfin=utcdatetime.UTCDateTime(options.endtime)
  inc=20
  kurt_window=3.0
  tref=time.time()

  for sta_line in sta_lines:
    sta=sta_line.split()[1]
    print "##############", sta, "##############"

    filepath="%s/*%s%s"%(data_dir,sta,options.data_glob)
    kurtpath="%s/*%s%s"%(data_dir,sta,options.kurt_glob)
    new_file="%s/kurto_%s.sac"%(data_dir,sta)
    t1=tdeb
    a=np.zeros((tfin-tdeb)*100)

    while t1 < tfin:
      if t1 < tfin-inc:
        t2=t1+inc
      else:
        t1=tfin-inc
        t2=tfin

      wf=Waveform()
      wf.read_from_file(filepath,starttime=t1,endtime=t2)
      dt=wf.delta
      x=wf.values

      # find indexes
      i1=int((t1-tdeb)*1./dt)
      i2=int((t2-tdeb)*1./dt)

      N=len(x)
      N2=np.log2(N)-7
      nlevel=int(np.fix(N2))
      c,flower,fupper = Fast_Kurtogram(x, nlevel,options.verbose,Fs=1/dt,opt2=1)

      wf.bp_filter(flower,fupper)
      wf.process_kurtosis(kurt_window)
      if i1 == 0:
        a[i1:i2]=wf.values[1:]
      elif i2 == (tfin-tdeb)*100:
        a[i1:i2]=wf.values
      else:
        a[i1:i2+1]=wf.values

      t1=t1+inc

    wf=Waveform()
    wf.read_from_file(kurtpath,starttime=tdeb,endtime=tfin)
    wf.values[:]=a
    wf.write_to_file_filled(new_file,format='SAC')

  print "temps de calcul:",time.time()-tref
