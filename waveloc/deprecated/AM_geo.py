#!/usr/bin/env python
# encoding: utf-8

"""
Classes pertaining to geographical calculations.

Created by Alessia Maggi 
"""
# original import
import sys, os, numpy 
#
# this fakes the import of pylab as it was originally
#import matplotlib.pyplot as pylab
#
from numpy import sqrt,arange,pi,array,cos,sin,dot,matrix,mat,rint,linspace,average,median,std
#from statlib import stats
#import scipy.stats as ss
#from scipy import interpolate
#
#from pylab import *
from mpl_toolkits.basemap.pyproj import Proj
from mpl_toolkits.basemap.pyproj import Geod
#from mpl_toolkits.basemap import Basemap as Basemap
#from matplotlib.mlab import griddata
#from waveloc_funcs import *
#from flexwin_funcs import *

proj_default=Proj(proj='utm',zone='10',ellps='WGS84')
geod_default=Geod(ellps='WGS84')

class GeoPoint(object):
  """
  A point in 3D space, as general as can be.
  Abstract class - do not instantiate.

  Contains methods common to daughter classes GeoPointXY and GeoPointLL.
  """

  def __init__(self):
    raise NotImplemented

  def dist_km_from(self,alat,alon):
    """
    Returns the distance in km from a latitude and longitude.
    
    Parameters:
    alat          latitude
    alon          longitude
    """
    (az12,az21,dist) = self.geod.inv(self.lon,self.lat,alon,alat)
    return dist/1000.0

  def dist3D_km_from(self,alat,alon,adepth):
    """
    Returns the distance in km from a latitude and longitude and depth.
    Approximation valid only for near-surface depths.
    
    Parameters:
    alat          latitude
    alon          longitude
    adepth        depth (meters)
    """
    (az12,az21,dist) = self.geod.inv(self.lon,self.lat,alon,alat)
    x=dist/1000.0
        
    return sqrt(x*x+(self.z-adepth)*(self.z-adepth))


  def azimuth_from(self,alat,alon):
    """
    Returns the azimuth in degrees from a latitude and longitude.
    
    Parameters:
    alat          latitude
    alon          longitude
      
    """
    (az12,az21,dist) = self.geod.inv(self.lon,self.lat,alon,alat)
    return az21
    

  def azimuth_to(self,alat,alon):
    """
    Returns the azimuth in degrees to a latitude and longitude.
    
    Parameters:
    alat          latitude
    lon          longitude
      
    """
    (az12,az21,dist) = self.geod.inv(self.lon,self.lat,alon,alat)
    return az12


class GeoPointXY(GeoPoint):
  """
  A point, to be initialized using x,y coordinates.
 
  Properties:
  x		x coordinate in meters (Easting)
  y		y coordinate in meters (Northing)

  Optional properties:
  z		z coordinate in meters (z positive into the Earth), defaults to 0.0
  value         a value associated to the point (any type or object), defaults to None
  proj 		a pyproj.Proj projection (defaults to UTM, zone 10)
  geod		a proj.Geod object, defaults to WGS84 ellipsoid

  Derived properties:
  lat		latitude (degrees N)
  lon		longitude (degrees E)
  """
  
  def __init__(self,x,y,z=0.0,value=None,proj=proj_default,geod=geod_default):
    self.x=x
    self.y=y
    self.z=z
    self.value=value
    self.proj=proj
    self.geod=geod

  def _get_lat_(self):
    lon,lat=self.proj(self.x,self.y,inverse='True')
    return lat

  def _get_lon_(self):
    lon,lat=self.proj(self.x,self.y,inverse='True')
    return lon

  lat=property(_get_lat_, doc='latitude')
  lon=property(_get_lon_, doc='longitude')


class GeoPointLL(GeoPoint):
  """
  A point, to be initialized using lat,lon coordinates.
 
  Properties:
  lat		latitude (degrees N)
  lon		longitude (degrees E)

  Optional properties:
  z		z coordinate in meters (z positive into the Earth)
  value         a value associated to the point (any type or object)
  proj 		a pyproj.Proj projection (defaults to UTM, zone 10)
  geod		a proj.Geod object, defaults to WGS84 ellipsoid

  Derived properties:
  x		x coordinate in meters (Easting)
  y		y coordinate in meters (Northing)
  """
  
  def __init__(self,lat,lon,z=0.0,value=None,proj=proj_default,geod=geod_default):
    self.lat=lat
    self.lon=lon
    self.z=z
    self.value=value
    self.proj=proj
    self.geod=geod

  def _get_x_(self):
    x,y=self.proj(self.lat,self.lon)
    return x

  def _get_y_(self):
    x,y=self.proj(self.lat,self.lon)
    return y

  x=property(_get_x_,doc='x coordinate in meters (Easting)')
  y=property(_get_y_,doc='y coordinate in meters (Northing)')
  


if __name__=='__main__':

  #p = Proj(proj='utm',zone=40,ellps='WGS84')
  p = Proj(init='epsg:2975')  # For Reunion Island : RGR92 / UTM zone 40S
  g = Geod(ellps='WGS84') # Use WGS84 ellipsoid.
 
  #x=360212 
  #y=7650290 

  P=GeoPointXY(360212,7650290,2378,proj=p)
  print P.lat, P.lon

  print P.dist_km_from(55.0, -21.0)
