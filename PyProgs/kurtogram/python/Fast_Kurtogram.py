#!/usr/bin/env python
# encoding: utf-8

import numpy as np
import scipy.signal as si
import logging
import matplotlib.pyplot as plt
import sys

eps = np.finfo(float).eps

def nextpow2(n):
    m_f = np.log2(n)
    m_i = np.ceil(m_f)
    return 2**m_i
    
def Fast_Kurtogram(x, nlevel, Fs=1, opt1=None, opt2=None):
    # Fast_Kurtogram(x,nlevel,Fs)
    # Computes the fast kurtogram of signal x up to level 'nlevel'
    # Maximum number of decomposition levels is log2(length(x)), but it is 
    # recommed to stay by a factor 1/8 below this.
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
        fc = .4					# a short filter is just good enough!
        h = si.firwin(N+1,fc) * np.exp(2*1j*np.pi*np.arange(N+1)*0.125)
        n = np.arange(2,N+2)
        print n
        g = h[(1-n)%N]*(-1)**(1-n)
        N = np.fix((3./2.*N))
        print N
        h1 = si.firwin(N+1,2./3*fc)*np.exp(2j*np.pi*np.arange(N+1)*0.25/3.)
        #~ plt.plot(h1)
        #~ plt.show()
        h2 = h1*np.exp(2j*np.pi*np.arange(N+1)/6.)
        h3 = h1*np.exp(2j*np.pi*np.arange(N+1)/3.)  
        
        if opt2 == 1:
            Kwav = K_wpQ(x,h,g,h1,h2,h3,nlevel,'kurt2')				# kurtosis of the complex envelope
        #~ else:
            #~ Kwav = K_wpQ(x,h,g,h1,h2,h3,nlevel,'kurt1')				# variance of the envelope magnitude

        
        # keep positive values only!
        Kwav[Kwav <= 0] = 0
        fig = plt.figure()
        #~ plt.subplot(ratio='auto')
        Level_w = np.arange(1,nlevel+1)
        Level_w = np.array([Level_w, Level_w + np.log2(3.)-1])
        Level_w = sorted(Level_w.ravel())
        Level_w = np.append(0,Level_w[0:2*nlevel-1])
        freq_w = Fs*(np.arange(0,3*2.0**nlevel-1+1))/(3*2**(nlevel+1)) + 1.0/(3.*2.**(2+nlevel))
        plt.imshow(Kwav,aspect='auto',extent=(freq_w[0],freq_w[-1],Level_w[0],Level_w[-1]),interpolation='bilinear')
        
        index = np.argmax(Kwav)
        index = np.unravel_index(index,Kwav.shape)
        f1 = freq_w[index[1]]
        l1 = Level_w[index[0]+1]
        fi = (index[1])/3./2**(nlevel+1)
        fi += 2.**(-2-l1)
        print fi, l1, Fs*fi
        plt.colorbar()
        plt.show()
    else:
        logging.error('stft-based is not implemented')
    
    #Ajouter le signal filtering !
    c = [];
    #~ test = int(raw_input('Do you want to filter out transient signals from the kurtogram (yes = 1 ; no = 0): '))
    test = 1
    #~ fi = fi * Fs
    lev = l1
    while test == 1:
        #~ fi = float(input(['	Enter the optimal carrier frequency (btw 0 and ',num2str(Fs/2),') where to filter the signal: ']));
        #~ fi = fi/Fs;
        #~ if opt1 == 1:
            #~ lev = input(['	Enter the optimal level (btw 0 and ',num2str(nlevel),') where to filter the signal: ']);
        if opt2 == 1:
            c,Bw,fc = Find_wav_kurt(x,h,g,h1,h2,h3,nlevel,lev,fi,'kurt2',Fs)
        #~ else
            #~ [c,Bw,fc] = Find_wav_kurt(x,h,g,h1,h2,h3,nlevel,lev,fi,'kurt1',Fs);
        test = int(raw_input('Do you want to keep on filtering out transients (yes = 1 ; no = 0): '))

    

def K_wpQ(x,h,g,h1,h2,h3,nlevel,opt,level=0):
    # K = K_wpQ(x,h,g,h1,h2,h3,nlevel)
    # Calculates the kurtosis K of the complete quinte wavelet packet transform w of signal x, 
    # up to nlevel, using the lowpass and highpass filters h and g, respectively. 
    # The WP coefficients are sorted according to the frequency decomposition.
    # This version handles both real and analytical filters, but does not yiels WP coefficients
    # suitable for signal synthesis.
    #
    # -----------------------
    # J Antoni : 12/2004 
    # Translation: T. Lecocq 02/2012
    # -----------------------   
    L = np.floor(np.log2(len(x)))
    if level == 0:
        if nlevel >= L:
            logging.error('nlevel must be smaller')
        level=nlevel
    x = x.ravel()
    print "THIS"
    print h, g
    KD, KQ = K_wpQ_local(x,h,g,h1,h2,h3,nlevel,opt,level)
    K = np.zeros((2*nlevel,3*2**nlevel))
    print "******************************************************"
    print KD.shape, KQ.shape, K.shape
    #~ K = KD
    for i in range(nlevel-1):
        print K[2*i,:].shape
        K[2*i,:] = KD[i+1,:]
        print K[2*i+1,:].shape
        K[2*i+1,:] = KQ[i,:]
       

    K[2*nlevel-1,:] = KD[nlevel,:]
    print "K Final Shape", K.shape
    return K

def K_wpQ_local(x,h,g,h1,h2,h3,nlevel,opt,level):
    print "LEVEL", level
    a,d = DBFB(x,h,g)
    
    N = len(a)
    d = d*np.power(-1.,np.arange(1,N+1))
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
    
    if level ==1:
        print "level = 1"
        K =np.array([K1*np.ones(3),K2*np.ones(3)]).flatten()
        print 'K.shape',K.shape
        KQ = np.array([Ka1,Ka2,Ka3,Kd1,Kd2,Kd3])
        print 'KQ.shape',KQ.shape
    if level > 1:
        print "entering rec with level %i"%(level-1)
        print "doing A"
        Ka,KaQ = K_wpQ_local(a,h,g,h1,h2,h3,nlevel,opt,level-1)
        print "doing D"
        Kd,KdQ = K_wpQ_local(d,h,g,h1,h2,h3,nlevel,opt,level-1)
        print "out of rec level %i" % (level -1)
        print Ka.shape, Kd.shape
        K1 = K1*np.ones(np.max(Ka.shape))
        K2 = K2*np.ones(np.max(Kd.shape))
        K12 = np.append(K1,K2)
        Kad = np.hstack((Ka, Kd))
        print ">", K12.shape, Kad.shape
        K = np.vstack((K12,Kad))

        Long = 2./6*np.max(KaQ.shape)
        Ka1 = Ka1*np.ones(Long)
        Ka2 = Ka2*np.ones(Long)
        Ka3 = Ka3*np.ones(Long)
        Kd1 = Kd1*np.ones(Long)
        Kd2 = Kd2*np.ones(Long)
        Kd3 = Kd3*np.ones(Long)
        tmp = np.hstack((KaQ,KdQ))
        print "HEEEERE"
        print tmp.shape
        KQ = np.concatenate((Ka1,Ka2,Ka3,Kd1,Kd2,Kd3))
        KQ = np.vstack((KQ, tmp))
        #~ if tmp.shape[0] != KQ.shape[0]:
            #~ tmp = tmp.T
        #~ for i in range(tmp.shape[0]):
            #~ KQ = np.vstack((KQ,tmp[i]))
        
        print "4", K.shape, KQ.shape
        

    
    if level == nlevel:
        K1 = kurt(x,opt)
        K = np.vstack((K1*np.ones(np.max(K.shape)), K))
        print "K shape", K.shape

        a1,a2,a3 = TBFB(x,h1,h2,h3)
        Ka1 = kurt(a1[len(h)-1:],opt)
        Ka2 = kurt(a2[len(h)-1:],opt)
        Ka3 = kurt(a3[len(h)-1:],opt)
        Long = 1./3*np.max(KQ.shape)
        Ka1 = Ka1*np.ones(Long)
        Ka2 = Ka2*np.ones(Long)
        Ka3 = Ka3*np.ones(Long)
        print KQ.shape
        tmp = np.array(KQ[0:-2])
        print "level==nlevel"
        
        KQ = np.concatenate((Ka1,Ka2,Ka3))
        KQ = np.vstack((KQ,tmp))
    
    print "i'm leaving level=%i and K.shape="%level,K.shape, "and KQ.shape=",KQ.shape
    return K, KQ

def kurt(x, opt):
    if opt=='kurt2':
        if np.all(x==0):
            K=0
            E=0
            return K
        x -= np.mean(x)
        E = np.mean(np.abs(x)**2)
        if E < eps:
            K=0
            return K
        K = np.mean(np.abs(x)**4)/E**2
        if np.all(np.isreal(x)):
            K = K - 3
        else:
            K = K - 2
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

def Find_wav_kurt(x,h,g,h1,h2,h3,nlevel,Sc,Fr,opt,Fs=1):
    # [c,Bw,fc,i] = Find_wav_kurt(x,h,g,h1,h2,h3,nlevel,Sc,Fr,opt2)
    # Sc = -log2(Bw)-1 with Bw the bandwidth of the filter
    # Fr is in [0 .5]
    #
    # -------------------
    # J. Antoni : 12/2004
    # -------------------
    level = np.fix((Sc))+ ((Sc%1) >= 0.5) * (np.log2(3)-1)
    Bw = 2**(-level-1)
    freq_w = np.arange(0,2**(level-1)) / 2**(level+1) + Bw/2.
    J = np.argmin(np.abs(freq_w-Fr))
    fc = freq_w[J]
    i = int(np.round(fc/Bw-1./2))
    if level % 1 == 0:
        acoeff = binary(i, level)
        bcoeff = np.array([])
        temp_level = level
    else:
        i2 = np.fix((i/3.))
        temp_level = np.fix((level))-1
        acoeff = binary(i2,temp_level)
        bcoeff = i-i2*3
    acoeff = acoeff[::-1]
    c = K_wpQ_filt(x,h,g,h1,h2,h3,acoeff,bcoeff,temp_level)
    print c
    kx = kurt(c,opt)
    
    print "kx", kx
    
    sig = np.median(np.abs(c))/np.sqrt(np.pi/2.)
    print sig
    threshold = sig*raylinv(np.array([.999,]),np.array([1,]))
    print "threshold", threshold
    spec = int(raw_input('	Do you want to see the envelope spectrum (yes = 1 ; no = 0): '))
    
    fig = plt.figure()
    t = np.arange(len(x))/Fs
    tc = np.linspace(t[0],t[-1],len(c))
    plt.subplot(2+spec,1,1)
    plt.plot(t,x,'k',label='Original Signal')
    plt.subplot(2+spec,1,2)
    plt.plot(tc,np.abs(c),'k')
    plt.plot(tc,threshold*np.ones(len(c)),'--r')
    #~ plt.title('Envlp of the filtr sgl, Bw=Fs/2^{'+(level+1)+'}, fc='+(Fs*fc)+'Hz, Kurt='+(np.round(np.abs(10*kx))/10)+', \alpha=.1%']
    plt.xlabel('time [s]')
    if spec == 1:
        print nextpow2(len(c))
        nfft = int(nextpow2(len(c)))
        env = np.abs(c)**2
        S = np.abs(np.fft.fft(env.ravel()-np.mean(env)*np.hanning(len(env))/len(env),nfft))
        f = np.linspace(0, 0.5*Fs/2**level,nfft/2)
        plt.subplot(313)
        plt.plot(f,S[:nfft/2],'k')
        plt.title('Fourier transform magnitude of the squared envelope')
        plt.xlabel('frequency [Hz]')
        plt.show()
    return [c,Bw,fc]

def binary(i,k):
    # return the coefficients of the binary expansion of i:
    # i = a(1)*2^(k-1) + a(2)*2^(k-2) + ... + a(k)

    if i>=2**k:
        logging.error('i must be such that i < 2^k !!')
    
    a = np.zeros(k+1)
    print a.shape
    temp = i
    for l in np.arange(k,0,-1):
        a[k-l] = np.fix(temp/2**l)
        temp = temp - a[k-l]*2**l
    return a
    
def K_wpQ_filt(x,h,g,h1,h2,h3,acoeff,bcoeff,level=0):
    # c = K_wpQ_filt(x,h,g,h1,h2,h3,acoeff,bcoeff,level)
    # Calculates the kurtosis K of the complete quinte wavelet packet transform w of signal x, 
    # up to nlevel, using the lowpass and highpass filters h and g, respectively. 
    # The WP coefficients are sorted according to the frequency decomposition.
    # This version handles both real and analytical filters, but does not yiels WP coefficients
    # suitable for signal synthesis.
    #
    # -----------------------
    # J Antoni : 12/2004 
    # -----------------------   
    nlevel = len(acoeff)
    L = np.floor(np.log2(len(x)))
    if level==0:
        if nlevel >= L:
            logging.error('nlevel must be smaller !!')
        level = nlevel
    x = x.ravel()
    if nlevel == 0:
        #if np.empty(bcoeff):
        if len(bcoeff)==0:
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

def  K_wpQ_filt_local(x,h,g,h1,h2,h3,acoeff,bcoeff,level):
    print level, x[:10]
    a,d = DBFB(x,h,g)         # perform one analysis level into the analysis tree
    N = len(a)                       
    d = d*np.power(-1.,np.arange(1,N+1))
    level = int(level)
    if level == 1:
        #if np.empty(bcoeff):
        if len(bcoeff)==0:
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
        print "acoeff", acoeff[level-1]
        if acoeff[level-1] == 0:
            c = K_wpQ_filt_local(a,h,g,h1,h2,h3,acoeff,bcoeff,level-1)
        else:
            c = K_wpQ_filt_local(d,h,g,h1,h2,h3,acoeff,bcoeff,level-1)
    print 'kurt', kurt(c,'kurt2')
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
    print k
    if len(k)!=0:
        tmp  = Inf
        x[k] = tmp(len(k))

    k = np.where(((b > 0) & (p > 0) & (p < 1)))[0]
    print k
    
    if len(k)!=0:
        pk = p[k]
        bk = b[k]
        print pk, bk
        x[k] = np.sqrt((-2*bk ** 2) * np.log(1 - pk))
    return x


if __name__ == "__main__":
    from scipy.io.matlab import loadmat
    v1 = loadmat(r"../Matlab/VOIE1.mat")
    x = v1['v1']
    Fs = 1
    nlevel= 7
    c = Fast_Kurtogram(x, nlevel, Fs)
