#!/usr/bin/env python
# encoding: utf-8
"""
AM_subs.py

Contains generally useful functions 

Created by Alessia Maggi on 2006-12-31.
Additions by Joanne BUCKENMEYER
Copyright (c) 2007 Alessia Maggi, Joanne Buckenmeyer. All rights reserved.

Functions in this module:
jday_to_month_day()
month_day_to_jday()
string_to_datetime()
string_to_datetime_jdy()
hilbert()
displacement_to_strain()
"""
import copy, datetime
#######################################################

def displacement_to_strain(x1,x0,x2,xs,seis1,seis0,seis2):

  npts=seis0.npts

  out_disp=copy.deepcopy(seis0)
  out_strain=copy.deepcopy(seis0)

  r0=x2-x1
  r1=x0-x1
  r2=x2-x0
  dx=xs-x0

  try:
    c0 = 1.+dx*(r2-r1-dx)/r1/r2
    c1 = dx*(-r2+dx)/r0/r1
    c2 = dx*( r1+dx)/r0/r2

    k0 = (r2-r1-2.*dx)/r1/r2
    k1 = (-r2 + 2.*dx)/r0/r1
    k2 = ( r1 + 2.*dx)/r0/r2


  except ZeroDivisionError:
    print "Bad geometry, x1,x0,x2 = %f %f %f"%(x1,x0,x2)

  for i in range(npts):
    out_disp.u_t[i]   = c1*seis1.u_t[i] + c0*seis0.u_t[i] + c2*seis2.u_t[i]
    out_strain.u_t[i] = k1*seis1.u_t[i] + k0*seis0.u_t[i] + k2*seis2.u_t[i]

  return (out_disp,out_strain)

#######################################################

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
  epoch_dtime=datetime.datetime(1970,1,1)
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

  return datetime.datetime(year,month,day,hour,min,sec)

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

  return datetime.datetime(year,month,day,hour,min,sec)

#######################################################
def datetime_to_string(date_time):
  """ 
  Takes a datetime object and returns a string formatted yyyymmddhhmmss
  """
  return "%04d%02d%02d%02d%02d%02d"%(date_time.year,date_time.month,\
    date_time.day,date_time.hour,date_time.minute,int(date_time.second))

def datetime_to_string_jdy(date_time):
  """ 
  Takes a datetime object and returns a string formatted yyyyjjjddhhmmss
  """
  return "%04d%03d02d%02d%02d"%(date_time.year,\
    month_day_to_jday(date_time.year,date_time.month),date_time.hour,\
    date_time.minute,int(date_time.second))
#######################################################

def hilbert(data,N):
  """
  Calculates the hilbert transform of a time-series, using numpy.fft
  on a size N array.  Returns a real array.

  Usage h = hilbert(data,N)
  """

  from numpy import fft 
  from numpy import imag 

  # calculate ft of data
  ft1=fft.fft(data,N)
  # zero out the negative frequency bins
  ft1[N/2+1:]=0.0
  # double the magnitude of the remaining spectrum
  ft1=2*ft1
  # inverse ft
  h1=fft.ifft(ft1,N)
  # return the imaginary part of the resulting array
  return imag(h1)


  ###########################################################3

def my_arccos(arg):
  import numpy 
  if arg < -1.0:
    arg=-0.9999999
  elif arg > 1.0:
    arg=0.9999999
  return numpy.arccos(arg)
    

def distaz(elat,elon,slat,slon):
  import numpy 
  """
  Sets azimuth (azm), backazimuth (bzm), distance in degrees (ddg),
  and distance in km (dkm) from local elat, elon, slat, slon

  Usage : (azm, bzm, ddg, dkm) = distaz(elat,elon,slat,slon)
  """

  c1=57.29578
  c2=1.570796
  c3=.9932313
  c4=.0174533
  a=c2-numpy.arctan2(c3*numpy.tan(c4*elat),1)
  b=-elon
  c=c2-numpy.arctan2(c3*numpy.tan(c4*slat),1)
  d=-slon
  delo=b-d
  if (delo < -180.) :
    delo=360.+delo
  if (delo >  180.) :
    delo=delo-360.
  delo=delo*c4
  del1=my_arccos(numpy.cos(a)*numpy.cos(c)+numpy.sin(a)*numpy.sin(c)*numpy.cos(delo))
  ddg=c1*del1
  dkm=(6371.227*(1.+.0033785*(1./3.-numpy.cos((a+c)/2.)**2)))*del1
  e=my_arccos((numpy.cos(c)-numpy.cos(a)*numpy.cos(del1))/(numpy.sin(a)*numpy.sin(del1)))
  s=my_arccos((numpy.cos(a)-numpy.cos(c)*numpy.cos(del1))/(numpy.sin(c)*numpy.sin(del1)))
  if (delo < 0.) :
    azm=360.-c1*e
    bzm=c1*s
  else :
    azm=c1*e
    bzm=360.-c1*s
			
  return (azm,bzm,ddg,dkm)



def latlong_at_distaz(latE,longE,azimut,distance):
  from math import pi,sin,cos,asin
  """ 
  Recherche des coordonnées d'un point en rapport avec un second :
  La fonction renvoi la latitude et la longitude du point cherché.
  La fonction connait la latitude et la longitude d'un point ainsi
  que la distance et l'azimut avec le point cherché.

  Usage (latS,lonS)=latlong(latE,longE,azimut,distance_km)
  """

  R = 6371            # Rayon de la Terre (en km)
  d = 1.* distance/R  # distance en radian sur la surface de la Terre

  azimut=azimut*pi/180.0

  latE=latE*pi/180.0
  longE=longE*pi/180.0

  latS =  asin(cos(d)*sin(latE)+sin(d)*cos(latE)*cos(azimut))
  longS = longE + asin(sin(azimut)*sin(d)/cos(latS))  
  
  latS=latS*180.0/pi
  longS=longS*180.0/pi

  return (latS, longS)




if __name__ == '__main__':
  pass
