"""
Quick and dirty functions for time-series manipulation.
"""
from datetime import datetime, timedelta
from numpy import average,empty,zeros,ones,std
from numpy.random import ranf
from filters import stalta
#from statlib import stats
#import scipy.stats as ss


#from grids_paths import GeoPoint


def first_sample(a):
    """
    Finds  the standard dev. only of the first Tm seconds from the timeseries.
    """
    smp0=a[0]

    return(smp0)

def initial_std(a,dt,Tm):
    """
    Finds  the standard dev. only of the first Tm seconds from the timeseries.
    """
    nsamps=int(Tm/dt)
    #std=stats.lstdev(a[:nsamps])
    s=std(a[:nsamps])

    return(s)

def initial_mean(a,dt,Tm):
    """
    Finds  the mean only of the first Tm seconds from the timeseries.
    """
    nsamps=int(Tm/dt)
    mean=average(a[:nsamps])

    return(mean)
    
def rmean_ini(a,dt,Tm):
    """
    Removes the mean only of the first Tm seconds from the timeseries.
    """
    nsamps=int(Tm/dt)
    mean=average(a[:nsamps])
    a=a-mean

    return(a)

def rmean(a):
    """
    Removes the mean of the time series.
    """
    mean=average(a)
    a=a-mean

    return(a)



def taper(x,tb=5.0):
    """
    Applies a sin taper symmetrically to the start and end of a timeseries.
    """
    from numpy import pi,sin
    te=tb
    pio2= pi/2.
    npts=len(x)
 
    # Apply taper to start of signal
    n=int(npts * tb/100.)
    for i in range(n):
        arg=pio2 * i / float(n)
        x[i]=x[i]*sin(arg)

    # Apply taper to end of signal
    n=int(npts * te/100.)
    for i in range(n):
        arg=pio2 * float(i)/float(n)
        x[(npts-1)-i] = x[(npts-1) - i] * sin(arg)

    return(x)

  
    
def noise_padding(x,dt,Tpad,amp=1e-8):
    """
    Adds random noise at the beginning of a timeseries
    to avoid 0 division in  rmean.  
    """

    # Create the random noise padding
    Npad=int(Tpad/dt)
    padding=zeros(Npad)
    for i in range(Npad):
        padding[i]=(ranf() - 0.5) * amp


    # insert the padding before the timeseries
    Nxx=len(x) + Npad
    xx=empty(Nxx)
    xx[0:Npad]=padding
    xx[Npad:]=x

    return(xx)

def constant_padding(x,dt,Tpad,amp=1e-6):
    """
    Adds a constant value to the (tapered) timeseries
    to avoid jumps in stalta.
    """

    # Create the constant noise padding
    Npad=int(Tpad/dt)
    padding=zeros(Npad)


    # insert the padding before the timeseries
    Nxx=len(x) + Npad
    xx=empty(Nxx)
    xx[0:Npad]=padding
    xx[Npad:]=x
    xx=xx+amp

    return(xx)


def stalta_proc(x,dt,Tpad,stawin,ltawin):
    """
    Processes a timeseries using a short term to long term average ratio.
    
    Adds a noise padding of length Tpad.  Uses filters.stalta to calculate
    the stalta.
    """
    x=rmean_ini(x,dt,3.0)
    x=taper(x)

    x = noise_padding(x,dt,Tpad)

    xs=stalta(x,dt,stawin,ltawin)

    return(xs)


def direction(a,b,w_level_dir):
  """
  Returns the direction from a to b  'UP', 'DOWN', 'FLAT' 
  """
  EPS=1.0e-8
  
  if  a < w_level_dir or b < w_level_dir  : 
    direction = 'FLAT'
  elif  abs(a-b) < EPS :
    direction = 'FLAT'
  elif b > a :
    direction = 'UP'
  else :
    direction = 'DOWN'
  
  return(direction)
  
  
def find_maxima(x,w_level):
  maxima=[]

  prev_dir=direction(x[0],x[1], w_level[0])
  
  for i in range(1,len(x)-1):
    next_dir=direction(x[i],x[i+1],w_level[i])
    if (next_dir == 'DOWN' and prev_dir == 'UP') :
      maxima.append(i)
    prev_dir=next_dir

  return(maxima)

def find_minima(x,w_level):
  minima=[]

  prev_dir=direction(x[0],x[1], w_level[0])
  
  for i in range(1,len(x)-1):
    next_dir=direction(x[i],x[i+1],w_level[i])
    if (next_dir == 'UP' and prev_dir == 'DOWN') :
      minima.append(i)
    prev_dir=next_dir

  return(minima)


# given an array, the index of a local maximum, and a water level, give
# the index of the closest points in the array that dip below the water
# level
def bracket_maximum(x, imax , w_level) :

  i1=0
  i2=len(x)-1

  for i in range(imax,0,-1):
    if(x[i]<w_level[i]) :
      i1=i
      break
    
  for i in range(imax,len(x)):
    if(x[i]<w_level[i]) :
      i2=i
      break
    
  return(i1,i2)
  

def jday_to_month_day(year,jday):
  """ Returns the calendar month (1-12) and day (1-31) for a given four-digit year
  and julian day (1-366).  Makes use of time.py.

  Usage: (month, day) = jday_to_month_day(year,jday)"""

  import time 

  mytime=time.strptime("%d %d"%(jday, year),"%j %Y")
  month=mytime.tm_mon
  day=mytime.tm_mday

  return (month,day)


#######################################################


def month_day_to_jday(year,month,day):
  """ Returns the julian day (1-366) for a given four-digit year, month (1-12) and day (1-31).
  Makes use of time.py.

  Usage: jday = month_day_to_jday(year,month,day)"""

  import time 

  mytime=time.strptime("%d %d %d"%(day,month, year),"%d %m %Y")
  jday=mytime.tm_yday

  return jday


#######################################################
def datetime_diff(dtime1,dtime2):
  epoch_dtime=datetime(1970,1,1)
  dt1=dtime1-epoch_dtime
  dt2=dtime2-epoch_dtime
  dt1_seconds=dt1.days*24*3600 + dt1.seconds + dt1.microseconds/1000000.0
  dt2_seconds=dt2.days*24*3600 + dt2.seconds + dt2.microseconds/1000000.0
  return (dt1_seconds-dt2_seconds)

#######################################################
def string_to_datetime(string):
  """
  Takes a string formatted as follows yyyymmddhhmmss and returns a datetime.datetime object.
  """
  year=int(string[0:4])
  month=int(string[4:6])
  day=int(string[6:8])
  hour=int(string[8:10])
  min=int(string[10:12])
  sec=int(string[12:14])

  return datetime(year,month,day,hour,min,sec)

def string_to_datetime_jdy(string):
  """
  Takes a string formatted as follows yyyyjjjhhmmss and returns a datetime.datetime object.
  """
  year=int(string[0:4])
  jday=int(string[4:7])
  hour=int(string[7:9])
  min=int(string[9:11])
  sec=int(string[11:13])
  
  (month,day)=jday_to_month_day(year,jday)

  return datetime(year,month,day,hour,min,sec)

#######################################################
def datetime_to_string(date_time):
  """ 
  Takes a datetime object and returns a string formatted yyyymmddhhmmss
  """
  return "%04d%02d%02d%02d%02d%02d%03d"%(date_time.year,date_time.month,\
    date_time.day,date_time.hour,date_time.minute,int(date_time.second),date_time.microsecond)

def datetime_to_string_jdy(date_time):
  """ 
  Takes a datetime object and returns a string formatted yyyyjjjddhhmmss
  """
  return "%04d%03d02d%02d%02d"%(date_time.year,\
    month_day_to_jday(date_time.year,date_time.month),date_time.hour,\
    date_time.minute,int(date_time.second))
#######################################################

def datetime_at_time(ref_time,t):
  delta_seconds=timedelta(seconds=t)
  return ref_time + delta_seconds


def cut_times(station, event, Vp, Vs, dtP = 10, dtS = 10):
  from grids_paths import GeoPoint
  """
  station and event are GeoPoints
  Calculates the cutting times for a station - event pair on the basis
  of a velocity model given by Vp and Vs
  Truncates the waveforms 10s before the P phase and 10s after 
  the S phase to account for uncertainties in location and velocity model
  """
  sta = GeoPoint(station.lat, station.lon)
  evt = GeoPoint(event.lat, event.lon)
  sta = station
  evt = event
  D = evt.dist_km_from(sta.lat, sta.lon)
  tP = D/Vp
  tS = D/Vs
  t1 = tP - dtP
  t2 = tS + dtS
  #t1 and t2 are offsets from the event origin time
  return (t1, t2)




