#!/usr/bin/env python
# encoding: utf-8

"""
Classes pertaining to grid and path calculations.

Created by Alessia Maggi and Alberto Michelini.
"""
import sys, os, numpy, glob, pickle, pysacio
from AM_geo import *
#
# this fakes the import of pylab as it was originally
import matplotlib.pyplot as pylab
#
from numpy import sqrt,arange,pi,array,cos,sin,dot,matrix,mat,rint,linspace,average,median,std,meshgrid
import scipy.stats as ss
from scipy import interpolate
#
from mpl_toolkits.basemap import Basemap as Basemap
from matplotlib.mlab import griddata
from waveloc_funcs import *
from flexwin_funcs import *
import logging
import numpy as np
import numexpr as ne


# For debugging reversal problem
from obspy.core import Stream, Trace

class Station(object):
  """
  A station.

  Properties:
    name        station name
    location	a GeoPoint location (either GeoPointLL or GeoPointXY)
    depth       depth of the station under the local earth surface (used for borehole stations)

  Derived properties
  x		x coordinate in meters (Easting)
  y		y coordinate in meters (Easting)
  lat		latitude in degrees
  lon		longitude in degrees
  elev	elevation in meters

  """
  def __init__(self,name,loc_tuple,depth=0.0,loc_type='latlon',proj=proj_default,geod=geod_default):
    """
    Initializes station:
    Inputs : 
    - name	 	station name
    - loc_tuple	(lat,lon,elev) or (x,y,elev) - elev is positive upwards
    - depth		depth of the station under local Earth's surface (used for borehole stations)
    - loc_type	'latlon' or 'xy' (defaults to 'latlon')
    - proj		a pyproj.Proj object
    - geod		a pyproj.Geod object
       
    """
    self.name=name
    self.depth=depth # depth is in meters under the local Earth's surface (used for borehole stations)
    
    # set location according to loc_type
    if loc_type=='latlon':
      self.location=GeoPointLL(loc_tuple[0],loc_tuple[1],loc_tuple[2],proj=proj,geod=geod)
    elif loc_type=='xy':
      self.location=GeoPointXY(loc_tuple[0],loc_tuple[1],loc_tuple[2],proj=proj,geod=geod)
    else:
      raise UserWarning('Location type %s unknown.'%loc_type) 

    # turn the elevation (positive) into z
    self.location.z = self.location.z * -1.0

  def _get_elev_(self):
    return -1*self.location.z
  def _get_x_(self):
    return self.location.x
  def _get_y_(self):
    return self.location.y
  def _get_lat_(self):
    return self.location.lat
  def _get_lon_(self):
    return self.location.lon

  elev=property(_get_elev_)
  x=property(_get_x_)
  y=property(_get_y_)
  lat=property(_get_lat_)
  lon=property(_get_lon_)
    

class Channel(Station):
  """
  A channel.  Uniquely identifies a data stream.

  Properties inherited from Station:
  name        station name
  location	a GeoPoint location (either GeoPointLL or GeoPointXY)
  depth       depth of the station under the local earth surface (used for borehole stations)

  Properties:
  comp	component name
  locid	location id (defaults to none)

  Derived properties (inherited from Station)
  x		x coordinate in meters (Easting)
  y		y coordinate in meters (Easting)
  lat		latitude in degrees
  lon		longitude in degrees
  elev	elevation in meters

 """

  def __init__(self,station,comp, locid=None):
    """
    Initializes channel from a station and a component.
    """
    self.location=station.location
    self.name=station.name
    self.depth=station.depth
    self.comp=comp
    self.locid=locid
  
  
class StationList(object):
  """
  Ensemble of Stations.
  
  Properties:
    stations        a dictionary of stations
    
  The keys of the stations dictionary are integers.
  """

  stations={} # stations are in a dictionary


  def __init__(self,proj=proj_default,geod=geod_default):
    """
    Initialize to an empty dictionary.
    """
    self.stations={}
    self.proj=proj
    self.geod=geod

  def _get_nsta_(self):
    return len(self.stations)
  nsta=property(_get_nsta_)
    

  def read_from_file(self,filename,sta_id_start=0):
    """
    Read the ensemble of stations from a file.
    
    Parameters:
      filename      name of file to be read
    
    Optional parameters:
      sta_id_start  starting key for the stations to be read in
                    default value is 0
    
    Keys will be sequential integers starting at sta_id_start.
    
    Expects the following file format (treated as whitespece
    separated columns):
    LOCSRCE ACER LATLON 40.786700 15.942700 0 0.690000
    or
    GTSRCE FJS  XYZ 367.405 7650.849 0.0 2.123
                        "lat"  "lon"  "depth" "elevation in km"
    
    """
    
    # clear the dictionary
    self.stations.clear() 
    
    # initialize the key counter
    sta_id=sta_id_start-1

    # read the file and populate the dictionary   
    file=open(filename,'r')
    lines=file.readlines()
    for line in lines:
      vals=line.split()
      name=vals[1]
      depth=float(vals[5])
      # column 2 contains indication of type of location
      if vals[2]=='LATLON': 
        loc_type='latlon'
        sta=Station(name,(float(vals[3]),float(vals[4]),float(vals[6])*1000),depth=depth,loc_type=loc_type,proj=self.proj,geod=self.geod)
      elif vals[2]=='XYZ':  
        loc_type='xy'
        sta=Station(name,(float(vals[3])*1000,float(vals[4])*1000,float(vals[6])*1000),depth=depth,loc_type=loc_type,proj=self.proj,geod=self.geod)
      else : 
        raise UserWarning('Location type %s not recognised in file %s.'%(vals[2],filename))
        
      sta_id=sta_id+1
      self.stations[sta_id]=sta
    
    file.close()


  
  def id_by_station_name(self,name):
    """
    Return the integer key corresponding to a station name.
    
    Name must be stripped of all whitespace to match station names in dictionary.
    """
    A=[ (sta.name==name,sta_id) for sta_id,sta in self.stations.iteritems() ] 
    A.sort()
    A.reverse()
    # Now the first element contains the station id for the name as the second element
    # of the tuple, unless the first element of the first element of the tuple is false
    ok=A[0][0]
    if ok:
      return A[0][1]
    else: 
      raise UserWarning("Station name %s not present in station list" % name)


  def filter_by_station_names(self,name_list):
    
    new_list=StationList(self.proj,self.geod)
    
    for name in name_list:
      try:
        name_key=self.id_by_station_name(name)
        new_list.stations[name_key]=self.stations[name_key]
      except:
        print "Station name %s not present in station list.  Ignoring." % name
      
    return new_list    
      

  def write_station_file(self,filename,loc_type='latlon'):
    """
    write stations properties in a NLLoc file
    """
    
    file=open(filename,'w')
    if loc_type=='latlon':
      for sta_id,sta in self.stations.iteritems():
        file.write("GTSRCE %05s LATLON %10.4f %10.4f %4.1f %6.4f\n"%\
        (sta.name,sta.lat,sta.lon,sta.depth,sta.elev/1000.0))
    elif loc_type=='xy':
      for sta_id,sta in self.stations.iteritems():
        file.write("GTSRCE %05s XYZ %10.4f %10.4f %4.1f %6.4f\n"%\
        (sta.name,sta.x/1000.0,sta.y/1000.0,sta.depth,sta.elev/1000.0))
    else:
      raise UserWarning('Unknown loc_type %s in write_station_file.'%loc_type)

    file.close()

  def display(self):
    """
    write stations properties to screen
    """
    
    for sta_id,sta in self.stations.iteritems():
      print("%05s %10.4fN %10.4fE %10.1f %10.1f %4.1f %6.1f"%\
      (sta.name,sta.lat,sta.lon,sta.x,sta.y,sta.depth,sta.elev))


  
class ChannelList(object):
  """
  Ensemble of Channels.
  
  Properties:
    channels        a dictionary of channels
    
  The keys of the channels dictionary are integers.
  """

  channels={} # channels are in a dictionary


  def __init__(self):
    """
    Initialize to an empty dictionary.
    """
    self.channels={}


  def _get_ncha_(self):
    return len(self.channels)
  ncha=property(_get_ncha_)

  def populate_from_station_list_and_data_files(self,sta_list,file_list,cha_id_start=0):

    from obspy.core import read

    self.channels.clear()
   
    cha_id=cha_id_start-1

    for filename in file_list : 
      st=read(filename)
      sta=st.traces[0].stats.station
      comp=st.traces[0].stats.channel
      locid=st.traces[0].stats.location
      try:
        sta_id=sta_list.id_by_station_name(sta)
        cha_id=cha_id+1
        cha=Channel(sta_list.stations[sta_id],comp,locid)
        self.channels[cha_id]=cha
      except UserWarning:
        logging.warn('Station name %s not present in station_list.'%sta)
        pass

    
  def populate_from_station_list(self,station_list,comp_string=["HHZ", "HHN", "HHE"],locid=None, cha_id_start=0):
    """
    Populate from a StationList object
    
    Parameters:
      station_list      the StationList to populate from
    
    Optional parameters:
      comp_string   space separated string of component names to add to 
                    the stations to generate channels
                    default value is "HHZ HHN HHE"
      locid		location id to give to the channel - defaults to None
    
    Keys will be sequential integers starting at cha_id_start.
    
    """
    
    # clear the dictionary
    self.channels.clear() 
    
    # initialize the key counter
    cha_id=cha_id_start-1
    

    # iterate through the stations list to populate channel list 
    for sta_id,sta in station_list.stations.iteritems():
      for comp in comp_string:
        cha_id=cha_id+1
        cha=Channel(sta,comp,locid)
        self.channels[cha_id]=cha
    
  def ids_by_station_and_component(self,name,comp):
    """
    Return the integer keys corresponding to a station name and component combination.
    
    Station and component names must be stripped of all whitespace to match station names in dictionary.
    """
    A=[ (cha.name==name and cha.comp==comp,cha_id) for cha_id,cha in self.channels.iteritems() ] 
    A.sort()
    A.reverse()
    # The first elements contain True or False statements, on which to choose the ids in the second elements.
    
    id_list=[]
    for a in A:
      # a[0] = True -> accept the cha_id
      if a[0]: 
        id_list.append(a[1])
    
    return id_list

  def id_by_station_and_component(self,name,comp):
    """
    Return the integer key corresponding to a station name and component combination if only one exists, otherwise raises an exception.
    
    Station and component names must be stripped of all whitespace to match station names in dictionary.
    """
    id_list=self.ids_by_station_and_component(name,comp)
    n_ids=len(id_list)
    if n_ids==0:
      raise UserWarning ('Station / component combination %s %s not found.'%(name,comp))
    elif n_ids==1:
      return id_list[0]
    else :
      raise UserWarning ('Found %d station / component combinations %s %s not found.'%(n_ids,name,comp))

 
  def ids_by_station(self,name):
    """
    Return the integer keys corresponding to a station name.
    
    Station name must be stripped of all whitespace to match station names in dictionary.
    """
    A=[ (cha.name==name,cha_id) for cha_id,cha in self.channels.iteritems() ] 
    A.sort()
    A.reverse()
    # The first elements contain True or False statements, on which to choose the ids in the second elements.
    
    id_list=[]
    for a in A:
      # a[0] = True -> accept the cha_id
      if a[0]: 
        id_list.append(a[1])
    
    return id_list


  def filter_by_station_names(self,name_list):
    return NotImplemented
    

    
class Grid(object):
  """
  A grid of geographical points.
    
  Properties:
  points          a dictionary of GeoPoints
  proj	    a geographical projection
  geod	    a geoid
  grid_type	    a grid type

  Derived properties:
  npts	    number of points in the grid
  min_lat         minimum latitude of the grid
  max_lat         maximum latitude of the grid
  min_lon         minimum longitude of the grid
  max_lon         maximum longitude of the grid
  min_x
  max_x
  min_y
  max_y
  min_z
  max_z
  min_value
  max_value
  
  The dictionary keys are integers.    
  """


  points={} # dictionary of GeoPoints
  
  def __init__(self,proj=proj_default,geod=geod_default):
    """
    Initialize to an empty grid.
    """
    self.points={}
    self.proj=proj
    self.geod=geod
    self.grid_type=None

  def _get_npts_(self):
    return len(self.points)
  def _get_min_lat_(self):
    lats = [p.lat for p in self.points.values()]
    return min(lats)
  def _get_max_lat_(self):
    lats = [p.lat for p in self.points.values()]
    return max(lats)
  def _get_min_lon_(self):
    lons = [p.lon for p in self.points.values()]
    return min(lons)
  def _get_max_lon_(self):
    lons = [p.lon for p in self.points.values()]
    return max(lons)

  def _get_min_x_(self):
    xs = [p.x for p in self.points.values()]
    return min(xs)
  def _get_max_x_(self):
    xs = [p.x for p in self.points.values()]
    return max(xs)
  def _get_min_y_(self):
    ys = [p.y for p in self.points.values()]
    return min(ys)
  def _get_max_y_(self):
    ys = [p.y for p in self.points.values()]
    return max(ys)
  def _get_min_z_(self):
    zs = [p.z for p in self.points.values()]
    return min(zs)
  def _get_max_z_(self):
    zs = [p.z for p in self.points.values()]
    return max(zs)

  def _get_min_value_(self):
    vals = [p.value for p in self.points.values()]
    return min(vals)
  def _get_max_value_(self):
    vals = [p.value for p in self.points.values()]
    return max(vals)

    
  npts=property(_get_npts_)
  min_lat=property(_get_min_lat_)
  max_lat=property(_get_max_lat_)
  min_lon=property(_get_min_lon_)
  max_lon=property(_get_max_lon_)
  min_x=property(_get_min_x_)
  max_x=property(_get_max_x_)
  min_y=property(_get_min_y_)
  max_y=property(_get_max_y_)
  min_z=property(_get_min_z_)
  max_z=property(_get_max_z_)
  min_value=property(_get_min_value_)
  max_value=property(_get_max_value_)


  def read_from_file(self,filename,loc_type='latlon'):
    """
    Read grid from an ascii file.
    
    Parameters:
      filename      name of file to read
      loc_type	    'latlon' or 'xy'
    
    Expects the following file format:
    loc_type = latlon
    '%d  %10.5f %10.5f %10.5f' = integer_key, longitude, latitude, depth (km)
    loc_type = xy
    '%d  %10.5f %10.5f %10.5f' = integer_key, x (km), y(km), depth (km)
    
    Calls append_from_file().
    """
    self.points={} # re_initialize list if it exists
    self.append_from_file(filename,loc_type)
    

  def append_from_file(self,filename,loc_type):
    """
    Appends points from an ascii file to the current grid.
    
    Parameters:
      filename      name of file to read
      loc_type	    'latlon' or 'xy'
    
    Expects the same file format as read_from_file().
    
    TODO: No checking is performed for double keys.  Must implement this sometime.
    """
    
    # read file and populate dictionary
    file=open(filename,'r')
    lines=file.readlines()

    if loc_type=='latlon':
      for line in lines:
        vals=line.split()
        point_id=int(vals[0])
        try : 
          point=GeoPointLL(float(vals[2]),float(vals[1]),float(vals[3])*1000,proj=self.proj,geod=self.geod) # beware: file is written in lon/lat format !!
        except IndexError: # grid only has two dimensions, so fix depth
          point=GeoPointLL(float(vals[2]),float(vals[1]),0.0,proj=self.proj,geod=self.geod) # beware: file is written in lon/lat format !!
        self.points[point_id]=point    

    elif loc_type=='xy':
      for line in lines:
        vals=line.split()
        point_id=int(vals[0])
        try : 
          point=GeoPointXY(float(vals[1])*1000,float(vals[2])*1000,float(vals[3])*1000,proj=self.proj,geod=self.geod) 
        except IndexError: # grid only has two dimensions, so fix depth
          point=GeoPointXY(float(vals[1])*1000,float(vals[2]),0.0,proj=self.proj,geod=self.geod) 
        self.points[point_id]=point    
 
    else:
      raise UserWarning('Unknown loc_type %s in Grid.append_from_file %f.'%(loc_type,filename))


    file.close()
    
    
  
  def write_to_file(self,filename,loc_type='latlon'):
    """
    Writes points from the grid to a file.
    
    Parameters:
      filename      name of file to write
      loc_type	    'latlon' or 'xy'
    
    Writes the same file format as read_from_file().
    """
    
    # write the file
    file=open(filename,'w')
    for point_id, point in self.points.iteritems():
      if loc_type=='latlon':
        file.write("%d  %10.5f %10.5f %10.5f\n"%(point_id,point.lon,point.lat,point.z))   
      elif loc_type=='xy':
        file.write("%d  %10.5f %10.5f %10.5f\n"%(point_id,point.x,point.y,point.z))   
      else:
        raise UserWarning('Unknown loc_type %s in Grid.write_to_file %f.'%(loc_type,filename))
    file.close()
    
    

  def append_to_file(self,filename,loc_type):
    """
    Appends points from the grid to a file.
    
    Parameters:
      filename      name of file to append to
      loc_type	    'latlon' or 'xy'
    
    Writes the same file format as read_from_file().
    """
    
    # append to file
    file=open(filename,'a')
    for point_id, point in self.points.iteritems():
      if loc_type=='latlon':
        file.write("%d  %10.5f %10.5f %10.5f\n"%(point_id,point.lon,point.lat,point.z))   
      elif loc_type=='xy':
        file.write("%d  %10.5f %10.5f %10.5f\n"%(point_id,point.x,point.y,point.z))   
      else:
        raise UserWarning('Unknown loc_type %s in Grid.write_to_file %f.'%(loc_type,filename))
    file.close()
    
    
  def read_from_NLL_files(self,filename,loc_type='xy'):
    """
    Constructs a regular grid from a NLL hdr file with the following format:
    201 151 101 357.185 7644.549 -3.000 0.100 0.100 0.100 SLOW_LEN
    and a NLL .buf file
   
    Parameters : 
      filename		name of .hdr file (without extension)
      loc_type		location type ('xy' or 'latlon') - default value 'xy'
    """
    from array import array 

    # open file and read header line
    hdr_filename="%s.hdr"%filename
    buf_filename="%s.buf"%filename
    f=open(hdr_filename)
    line=f.read()
    vals=line.split()
    f.close()


    # interpret line
    print vals
    nx=int(vals[0])
    ny=int(vals[1])
    nz=int(vals[2])
    x_orig=float(vals[3])
    y_orig=float(vals[4])
    z_orig=float(vals[5])
    dx=float(vals[6])
    dy=float(vals[7])
    dz=float(vals[8])
    self.grid_type=vals[9]

    # open and read buffer file
    f=open(buf_filename,'rb')
    buf=array('f')
    buf.fromfile(f,nx*ny*nz)
    f.close()

    # first point id = 0
    point_id=0
    # construct regular grid
    buf_id=0
    for i in range(nx):
      x=x_orig+i*dx
      for j in range(ny):
        y=y_orig+j*dy
        for k in range(nz):
          z=z_orig+k*dz
          value=buf[buf_id]
          if loc_type=='xy':
            # NLL XYZ format is in km for x y and z, while internal GeoPointXY format is in meters
            self.points[point_id]=GeoPointXY(x*1000,y*1000,z*1000,value=value,proj=self.proj,geod=self.geod)
            point_id=point_id+1
            buf_id=buf_id+1
          elif loc_type=='latlon':
            self.points[point_id]=GeoPointLL(y,x,z,value=value,proj=self.proj,geod=self.geod) # Beware : lat=y and lon=x ...
            point_id=point_id+1
            buf_id=buf_id+1
    
          
  def display_NLL_XYZ(self,display_value,display_type='z',filename='',title='',geography=False):
    """
    Makes a plot of the grid.  Assumes regular NLL-type grid in X Y Z space.

    Parameters:
    display_value : a constant x y or z value (in km) - the values plotted will be those of the
                    gridpoints closest to this constant value
    display_type  : indicates which of 'x', 'y' or 'z' is kept constant
    filename      : filename for output (if none, plot is sent to screen)
    title         : title of plot (defaults to no title)
    """

    # TODO : rewrite this to use interpolation onto a simple grid for plotting purposes - plots will be faster !!
    # NOTE : see examples in correlation grid plotting

    epsilon=1.0

    if display_type=='z':
      z=display_value*1000 # conversion in meters

      # order the points by distance to the display value
      A=[( abs(p.z-z), p ) for p in self.points.values()]
      A.sort()

      # retain only the points whose distance to the display value is within epsilon of the minimum distance
      dist=A[0][0] 
      B=[ a[1] for a in A if (a[0]-dist)<epsilon ]

      # start the plot
      pylab.clf()
      pylab.title(title)

      if geography : 
        border=0.08
        min_lon=self.min_lon-border
        max_lon=self.max_lon+border
        min_lat=self.min_lat-border
        max_lat=self.max_lat+border
     
        m=Basemap(llcrnrlon=min_lon,llcrnrlat=min_lat,\
                  urcrnrlon=max_lon,urcrnrlat=max_lat,resolution='f',\
                  projection='lcc',lon_0=(max_lon-min_lon)/2, lat_0=(max_lat-min_lat)/2)
        m.drawcoastlines()
        m.drawcountries()
        m.drawmapboundary(fill_color='white')
        m.fillcontinents(color='#cc9966',lake_color='#99ffff')
 
        m.drawparallels(arange(min_lat,max_lat,border),labels=[1,0,0,0])
        m.drawmeridians(arange(min_lon,max_lon,border),labels=[0,0,0,1])

        lons=[ b.lon for b in B ] 
        lats=[ b.lat for b in B ] 
        values=[ b.value for b in B ] 

        x,y=m(lons,lats)
        m.contourf(x,y,values,cmap=pylab.cm.jet)

      else:
        x=[ b.x for b in B ] 
        y=[ b.y for b in B ] 
        values=[ b.value for b in B ] 
       
        pylab.contourf(x,y,values,cmap=pylab.cm.jet)
        

      if not filename=='':
        pylab.savefig(filename)
      else:
        pylab.show()



  def construct_regular_grid(self,lon0,lat0,rotang_deg,Xlength_km,Ylength_km,\
                             xstep_km=5,ystep_km=5,point_id_start=0,proj=proj_default):
    """
    Constructs a regular, rotated grid in latitude and longitude.
    
    Parameters:
    lon0          longitude of rotation point
    lat0          latitude of rotation point
    rotang_deg    rotation angle in degrees
    Xlength_km    length in km along the x axis
    Ylength_km    length in km along the y axis
    xstep_km      separation of points along the x axis
    ystep_km      separation of points along the y axis
    
    Optional parameters:
      point_id_start  key value (integer) for first point; is 0 by default
      
    Uses basemap.pyproj.Proj, with UTM projection onto ellipsoid WGS84.    
    """
    
    # fixup input values for subsequent calculations
    rotang_rad= rotang_deg * pi / 180.0
    Xlength=Xlength_km*1000.0
    Ylength=Ylength_km*1000.0
    Xstep=xstep_km*1000.0
    Ystep=ystep_km*1000.0
    point_id=point_id_start-1   

    # turn center lat, lon to X,Y values using projection
    X0,Y0 = p(lon0, lat0)
    XY0=array([[X0],[Y0]])
    
    # start from midddle of the axis
    Xmid=Xlength/2.0
    Ymid=Ylength/2.0
    
    # create grid and center it
    Xs=arange(X0-Xmid,X0-Xmid+Xlength,Xstep)
    Ys=arange(Y0-Ymid,Y0-Ymid+Ylength,Ystep)
    
    # set up rotation matrix
    R=array( [[ cos(rotang_rad), sin(rotang_rad) ],[-sin(rotang_rad),cos(rotang_rad)]] )
    
    # Create grid and rotate 
    for x in Xs:
      for y in Ys:
        # create the relative grid
        X=array([[x],[y]]) - XY0  
        # rotate it  
        Xrot=dot(R,X) + XY0
        xrot=Xrot[0,0]
        yrot=Xrot[1,0]
        # transform rotated grid to lat/lon coordinates
        (lonrot, latrot) = proj(xrot, yrot, inverse=True) # inverse transform
        # write to dictionary
        point_id=point_id+1
        self.points[point_id]=GeoPointLL(latrot,lonrot,proj=self.proj,geod=self.geod)
    

  def display_italy(self,title="Italy grid",filename="italy_grid.png"):
    """
    Plots grid to file.
    ONLY FOR ITALY GRID - THIS SHOULD DISAPPEAR !!
    """
  
    # begin plotting of the points
    # Lambert Conformal Conic map. 
    # Italy coordinates
    m = Basemap(llcrnrlon=5.,llcrnrlat=35.,urcrnrlon=25.,urcrnrlat=48.,
            projection='lcc',lat_1=38.,lat_2=45.,lon_0=10.,
            resolution ='l',area_thresh=1000.)
    # create figure.
    pylab.clf()
    pylab.title(title)
    m.drawcoastlines()
    m.drawcountries()
    m.drawmapboundary(fill_color='#99ffff')
    m.fillcontinents(color='#cc9966',lake_color='#99ffff')

    m.drawparallels(arange(35,50,2),labels=[1,1,0,0])
    m.drawmeridians(arange(5,24,2),labels=[0,0,0,1])

    # make simple lists of lons and lats
    lons=[]
    lats=[]
    for point_id, point in self.points.iteritems():
        lons.append(point.lon)
        lats.append(point.lat)
    x,y = m(lons,lats)
    m.plot(x,y,'ro',ms=1.0)
    #
    if not filename=="":
      pylab.savefig(filename)
    else:
     show()


class QDGrid(Grid):
  """
  A Quick and Dirty Grid - For regular xyz grids only (no geographical plotting etc.).
    
  Properties inherited from Grid:
    points          a dictionary of GeoPoints (may or may not be filled in correctly)
    proj	    = None
    geod	    = None
    grid_type	    a grid type

  Non-inherited properties:

    buf            sequential list of point values (x-y-z loop order)
    
    nx              number of x-points
    ny              number of y-points
    nz              number of z-points
    dx	    step between x-points
    dy	    step between y-points
    dz	    step betwee. z-points
    x_orig
    y_orig
    z_orig

  Derived properties:
    npts	    number of points in the grid
    min_lat         = None
    max_lat         = None
    min_lon         = None
    max_lon         = None
    min_x
    max_x
    min_y
    max_y
    min_z
    max_z
    min_value
    max_value
    
  The dictionary keys are integers.    
 
  """


  def __init__(self):
    """
    Initialize to an empty grid (overrides Grid.init).
    """
    self.points={}
    self.proj=None
    self.geod=None
    self.grid_type=None
    self.buf=None

  def _get_npts_(self):
    return len(self.buf)

  def _get_min_lat_(self):
    return None
  def _get_max_lat_(self):
    return None
  def _get_min_lon_(self):
    return None
  def _get_max_lon_(self):
    return None

  def _get_min_x_(self):
    return self.x_orig
  def _get_max_x_(self):
    return self.x_orig+self.nx*self.dx

  def _get_min_y_(self):
    return self.y_orig
  def _get_max_y_(self):
    return self.y_orig+self.ny*self.dy

  def _get_min_z_(self):
    return self.z_orig
  def _get_max_z_(self):
    return self.z_orig+self.nz*self.dz

  def _get_min_value_(self):
    return min(self.buf)
  def _get_max_value_(self):
    return max(self.buf)


  npts=property(_get_npts_)
  min_lat=property(_get_min_lat_)
  max_lat=property(_get_max_lat_)
  min_lon=property(_get_min_lon_)
  max_lon=property(_get_max_lon_)
  min_x=property(_get_min_x_)
  max_x=property(_get_max_x_)
  min_y=property(_get_min_y_)
  max_y=property(_get_max_y_)
  min_z=property(_get_min_z_)
  max_z=property(_get_max_z_)
  min_value=property(_get_min_value_)
  max_value=property(_get_max_value_)

  def read_NLL_hdr_file(self,filename):
    """
    Reads a NLL hdr file with the following format:
    201 151 101 357.185 7644.549 -3.000 0.100 0.100 0.100 SLOW_LEN
    and sets appropriate properties
    """
    f=open(filename)
    line=f.read()
    vals=line.split()
    f.close()


    # interpret line
    self.nx=int(vals[0])
    self.ny=int(vals[1])
    self.nz=int(vals[2])
    self.x_orig=float(vals[3])
    self.y_orig=float(vals[4])
    self.z_orig=float(vals[5])
    self.dx=float(vals[6])
    self.dy=float(vals[7])
    self.dz=float(vals[8])
    self.grid_type=vals[9]

 

  def read_from_NLL_files(self,filename):
    """
    Constructs a regular grid from a NLL hdr file with the following format:
    201 151 101 357.185 7644.549 -3.000 0.100 0.100 0.100 SLOW_LEN
    and a NLL .buf file

    Overrides Grid.read_from_NLL_files
   
    Parameters : 
      filename		name of .hdr file (without extension)
    """
    from array import array 

    # open file and read header line
    hdr_filename="%s.hdr"%filename
    buf_filename="%s.buf"%filename

    # read hdr file
    self.read_NLL_hdr_file(hdr_filename)

    # open and read buffer file
    f=open(buf_filename,'rb')
    self.buf=array('f')
    self.buf.fromfile(f,self.nx*self.ny*self.nz)
    f.close()
 
    #print "Read %d points from file %s."%(len(self.buf),buf_filename)


  def construct_empty_grid(self,value=0.0):
    try:
      self.buf=zeros(self.nx*self.ny*self.nz)+value
    except AttributeError:
      raise UserWarning ('Atempting to construct a grid without having the header information.')

  def get_ix_iy_iz(self,ib):
    ix=ib/(self.ny*self.nz)
    yz_part=ib%(self.ny*self.nz)
    iy=yz_part/self.nz
    iz=yz_part%self.nz
    return (ix,iy,iz)
    
  def get_grid_value(self,ix,iy,iz):
    ib=ix*self.ny*self.nz + iy*self.nz + iz
    return self.buf[ib]

  def set_grid_value(self,ix,iy,iz,value):
    ib=ix*self.ny*self.nz + iy*self.nz + iz
    self.buf[ib]=value

  def value_at_point(self,x,y,z,epsilon=0.001):
    """
    Performs n-linear interpolation on the regular grid
    """

    # sanity check for point being within grid
    # use a tollerance value to avoid problems with numerical errors
    if x < self.min_x-epsilon or x > self.max_x+epsilon \
    or y < self.min_y-epsilon or y > self.max_y+epsilon \
    or z < self.min_z-epsilon or z > self.max_z+epsilon :
      logging.debug('(x,y,z) = (%f, %f, %f) : (x,y,z)min = (%f, %f, %f) - (x,y,z)max = (%f, %f, %f)'%(x,y,z,self.min_x,self.min_y,self.min_z,self.max_x,self.max_y,self.max_z))
      raise UserWarning('Point (%f, %f, %f) is outside the grid (tollerance=%f).'%(x,y,z,epsilon))

    # fix up lower and upper bounds if they are still (just) outside the grid
    if x < self.min_x : x=self.min_x
    if y < self.min_y : y=self.min_y
    if z < self.min_z : z=self.min_z

    if x > self.max_x : x=self.max_x
    if y > self.max_y : y=self.max_y
    if z > self.max_z : z=self.max_z
    
    # make arrays of X, Y and Z ranges 
    X=numpy.arange(self.nx)*self.dx+self.x_orig
    Y=numpy.arange(self.ny)*self.dy+self.y_orig
    Z=numpy.arange(self.nz)*self.dz+self.z_orig

    # get the position this point would have in the X,Y,Z arrays if they were extended by 1
    ix=X.searchsorted(x)
    iy=Y.searchsorted(y)
    iz=Z.searchsorted(z)

    # set the interpolation "box" for extreme cases
    if self.nx==1 : # special case of 2D grid
      ix1=0
      ix2=0
    elif ix==0: # lower bound
      ix1=0
      ix2=1
    elif ix==self.nx: # upper bound
      ix1=self.nx-2
      ix2=self.nx-1
    else :	# general case
      ix1=ix-1
      ix2=ix

    if iy==0:	# lower bound
      iy1=0
      iy2=1
    elif iy==self.ny: # upper bound
      iy1=self.ny-2
      iy2=self.ny-1
    else :	# general case
      iy1=iy-1
      iy2=iy

    if iz==0:	# lower bound
      iz1=0
      iz2=1
    elif iz==self.nz: # upper bound
      iz1=self.nz-2
      iz2=self.nz-1
    else :	# general case
      iz1=iz-1
      iz2=iz

    # set up the values
    # bottom four values counterclockwise from x1y1
    v_x1y1z1=self.get_grid_value(ix1,iy1,iz1)
    v_x2y1z1=self.get_grid_value(ix2,iy1,iz1)
    v_x2y2z1=self.get_grid_value(ix2,iy2,iz1)
    v_x1y2z1=self.get_grid_value(ix1,iy2,iz1)
    # top four values counterclockwise from x1y1
    v_x1y1z2=self.get_grid_value(ix1,iy1,iz2)
    v_x2y1z2=self.get_grid_value(ix2,iy1,iz2)
    v_x2y2z2=self.get_grid_value(ix2,iy2,iz2)
    v_x1y2z2=self.get_grid_value(ix1,iy2,iz2)

    # set up interpolators
    # take extra care over the X interpolator in case of 2D grid
    if ix2==ix1:
      tx=0
    else:
      tx=(x-X[ix1])/(X[ix2]-X[ix1])
    ty=(y-Y[iy1])/(Y[iy2]-Y[iy1])
    tz=(z-Z[iz1])/(Z[iz2]-Z[iz1])

    # do bilinear interpolation
    result = (1-tx) * (1-ty) * (1-tz) * v_x1y1z1 + \
                tx  * (1-ty) * (1-tz) * v_x2y1z1 + \
                tx  *    ty  * (1-tz) * v_x2y2z1 + \
             (1-tx) *    ty  * (1-tz) * v_x1y2z1 + \
             (1-tx) * (1-ty) *    tz  * v_x1y1z2 + \
                tx  * (1-ty) *    tz  * v_x2y1z2 + \
                tx  *    ty  *    tz  * v_x2y2z2 + \
             (1-tx) *    ty  *    tz  * v_x1y2z2 

    return result
  
  def display_grid_zcut(self,filename,z=0,title=''):

    xarray=numpy.arange(self.nx)*self.dx+self.x_orig
    yarray=numpy.arange(self.ny)*self.dy+self.y_orig

    print "Making plot grid"
    values=[self.value_at_point(x,y,z) for y in yarray for x in xarray]
    vals=array(values)
    vals.shape=(self.ny,self.nx)

    print "Plotting grid to file %s"%filename
    pylab.clf()
    pylab.title(title)
    pylab.contourf(xarray,yarray,vals,cmap=pylab.cm.jet)
    pylab.colorbar()
    pylab.savefig(filename)

  def display_grid_ycut(self,filename,y=0,title=''):

    xarray=numpy.arange(self.nx)*self.dx+self.x_orig
    zarray=numpy.arange(self.nz)*self.dz+self.z_orig


    print "Making plot grid"
    values=[self.value_at_point(x,y,z) for z in zarray for x in xarray]
    vals=array(values)
    vals.shape=(self.nz,self.nx)

    print "Plotting grid to file %s"%filename
    pylab.clf()
    pylab.title(title)
    pylab.contourf(xarray,zarray,vals,cmap=pylab.cm.jet)
    pylab.colorbar()
    pylab.savefig(filename)

  def display_grid_xcut(self,filename,x=0,title=''):

    yarray=numpy.arange(self.ny)*self.dy+self.y_orig
    zarray=numpy.arange(self.nz)*self.dz+self.z_orig


    print "Making plot grid"
    values=[self.value_at_point(x,y,z) for z in zarray for y in yarray]
    vals=array(values)
    vals.shape=(self.nz,self.ny)

    print "Plotting grid to file %s"%filename
    pylab.clf()
    pylab.title(title)
    pylab.contourf(yarray,zarray,vals,25,cmap=pylab.cm.jet)
    pylab.colorbar()
    pylab.savefig(filename)


  def construct_points_from_buffer(self):
    # first point id = 0
    point_id=0
    # construct regular grid
    buf_id=0
    for i in range(self.nx):
      x=self.x_orig+i*self.dx
      for j in range(self.ny):
        y=self.y_orig+j*self.dy
        for k in range(self.nz):
          z=self.z_orig+k*self.dz
          value=self.buf[buf_id]
          # NLL XYZ format is in km for x y and z, while internal GeoPointXY format is in meters
          self.points[point_id]=GeoPointXY(x*1000,y*1000,z*1000,value=value,proj=None,geod=None)
          point_id=point_id+1
          buf_id=buf_id+1
 
class QDTimeGrid(QDGrid):
  """
  A regular geographical grid of time delays as dictionaries (one point per station)
  Inherits from QDGrid.

  """

  def construct_empty_grid(self,array_length=1):
    try:
      #self.buf=[numpy.zeros(array_length) for i in range(self.nx*self.ny*self.nz)]
      self.buf=[{} for i in range(self.nx*self.ny*self.nz)]
    except AttributeError:
      raise UserWarning ('Atempting to construct a grid without having the header information.')

  def dump_buffer_to_file(self,filename):
    f=open(filename,'w')
    pickle.dump(self.buf, f)
    f.close()

  def load_buffer_from_file(self,filename):
    f=open(filename,'r')
    self.buf=pickle.load(f)
    f.close()

  def populate_from_time_grids(self,grid_filename_base,channel_list,out_path,load_buf=False):

    time_grids={}

    #filename for temporary grid file
    grid_base=os.path.basename(grid_filename_base)
    tmp_buf_filename=os.path.join(out_path,"%s.search.ttimes"%grid_base)

    # read all the full-resolution NLL time files
    logging.debug('Reading full-resolution NLL time files')
    logging.debug('Channels contains %d channels'%len(channel_list.channels))
    for s in channel_list.channels.values():
      try:
        nll_grid_name="%s.%s.time"%(grid_filename_base,s.name)
        grid=QDGrid()
        grid.read_from_NLL_files(nll_grid_name)
        grid_id="%s.%s"%(s.name,s.comp)
        logging.debug("%s %d %s "%(nll_grid_name, grid.npts, grid_id))
        time_grids[grid_id]=grid
      except IOError:
        logging.error('Error reading files %s.* They may not exist.  Will ignore station %s for this run.'%(nll_grid_name,s.name))
        pass

    # set up smaller time grid on search grid only
    logging.info('Setting up local grid. This could take some time, be patient...')
    grid_ids=time_grids.keys()
    logging.debug("Grid keys : ")
    logging.debug(grid_ids)
    
    ixarray=range(self.nx)
    iyarray=range(self.ny)
    izarray=range(self.nz)


    # by default do not calculate grids
    calculate_grids=False

    # if we've been passed an option to force calculation, then go ahead and calculate
    if not load_buf : calculate_grids = True

    # Try to load file - if ok, then fine, else force recalculation
    if load_buf :
      logging.info('Atempting to read travel-time buffer from file %s'%tmp_buf_filename)
      try:
        self.load_buffer_from_file(tmp_buf_filename)
        logging.debug("Loaded grid keys:")
        logging.debug(self.buf[0].keys())
      except IOError:
        logging.info('Cannot load file %s : calculating grid and writing to file.'%tmp_buf_filename)
        calculate_grids=True
        
    # if we need to calculate the grids, just go ahead and do it
    if calculate_grids:
      self.construct_empty_grid()
      for ix in ixarray:
        logging.debug('In outer loop %d of %d...'%(ix,self.nx))
        x=ix*self.dx + self.x_orig
        for iy in iyarray:
          y=iy*self.dy + self.y_orig
          for iz in izarray:
            z=iz*self.dz + self.z_orig

            times={}
            for grid_id in grid_ids:
              times[grid_id]=time_grids[grid_id].value_at_point(x,y,z)
            self.set_grid_value(ix,iy,iz,times)
      logging.info('Writing travel-time buffer to file %s'%tmp_buf_filename)
      self.dump_buffer_to_file(tmp_buf_filename)
    del time_grids

    #logging.debug('Times for first point in geographical grid : %s '%self.buf[0])

  def read_station_coordinates_from_hdr_file(self,nll_grid_name,s_name):
    
    # format of file:
    #1 601 61  0.000000 0.000000 -1.000000  0.500000 0.500000 0.500000 TIME2D FLOAT
    #ANSA 34.064651 1.833333 -0.180000
    #TRANSFORM  SIMPLE LatOrig 50.650500  LongOrig 5.023300  RotCW 0.000000

    f=open("%s.hdr"%nll_grid_name,'r')
    lines=f.readlines()

    line=lines[1]
    sta_x=float(line.split()[1])
    sta_y=float(line.split()[2])
    sta_z=float(line.split()[3])

    return(sta_x,sta_y,sta_z)

  def populate_from_2D_time_grids(self,grid_filename_base,channel_list):

    time_grids={}

    # read all the full-resolution NLL time files
    logging.debug('Reading full-resolution NLL time files')
    logging.debug('Channels contains %d channels'%len(channel_list.channels))
    for s in channel_list.channels.values():
      try:
        nll_grid_name="%s.%s.time"%(grid_filename_base,s.name)
        grid=QDGrid()
        grid.read_from_NLL_files(nll_grid_name)
        # for the 2D grid case, we need the station coordinates in the same coordinates of the search grid
        # they can be found in the .time files
        sta_x,sta_y,sta_z = self.read_station_coordinates_from_hdr_file(nll_grid_name,s.name)

        grid_id="%s.%s"%(s.name,s.comp)
        logging.debug("%s %d %s %.2f %.2f %.2f"%(nll_grid_name, grid.npts, grid_id, sta_x, sta_y, sta_z))
        time_grids[grid_id]=grid
      except IOError:
        logging.error('Error reading file %s.  It may not exist.  Will ignore station %s for this run.'%(grid_filename_base,s.name))
        pass

    # set up smaller time grid on search grid only
    logging.info('Setting up local grid. This could take some time, be patient...')
    grid_ids=time_grids.keys()
    
    ixarray=range(self.nx)
    iyarray=range(self.ny)
    izarray=range(self.nz)


    self.construct_empty_grid()
    for ix in ixarray:
      logging.debug('In outer loop %d of %d...'%(ix,self.nx))
      x=ix*self.dx + self.x_orig
      for iy in iyarray:
        y=iy*self.dy + self.y_orig
        for iz in izarray:
          z=iz*self.dz + self.z_orig

          h_distance_km = sqrt((x-sta_x)**2 + (y-sta_y)**2)
          v_distance_km = z-sta_z
      
          times={}
          for grid_id in grid_ids:
            times[grid_id]=time_grids[grid_id].value_at_point(0,h_distance_km,v_distance_km)
          self.set_grid_value(ix,iy,iz,times)

    del time_grids

    #logging.debug('Times for first point in geographical grid : %s '%self.buf[0])

class QDStackGrid(object):
  def __init__(self,nx,ny,nz,nt):
    self.buf=np.zeros((nx,ny,nz,nt))

class QDCorrGrid(QDGrid):
  """
  A regular geographical grid of shifted and stacked waveforms
  Inherits from QDGrid.

  Note : uses int16 integers to stack correlation values !!
  """

  def construct_empty_grid(self,array_length=1):
    try:
      npts=self.nx*self.ny*self.nz
      self.buf=numpy.zeros(array_length*npts, dtype=numpy.int16)
      self.buf.shape=(npts,array_length)
    except AttributeError:
      raise UserWarning ('Atempting to construct a grid without having the header information.')

  def set_grid_value(self,ix,iy,iz,value=[]):
    ib=ix*self.ny*self.nz + iy*self.nz + iz
    print self.buf[ib].shape
    print numpy.array(value).shape
    self.buf[ib]=numpy.array(value)

  def write_grid_timeslice(self,itime,filename):
    data=self.buf[:,itime]
    data.tofile(filename)

    

class CorrGrid(Grid):
  """
  A geographical grid of averaged correlation functions.
    
  Properties inherited from Grid:
  points          a dictionary of GeoPoints
  min_lat         minimum latitude of the grid
  max_lat         maximum latitude of the grid
  min_lon         minimum longitude of the grid
  max_lon         maximum longitude of the grid
  
  Properties:
  base_grid       Grid object containing the geometry that
                  this grid will follow
  n_timelags      length of the cross_correlation timeseries
                  that form the value of each point in the 
                  CorrGrid
  b               time of first point in the correlation timeseries
  dt              timestep of the correlation timeseries
  ref_time        reference time (datetime object)
  path_length     dictionary of lengths of averaging paths through the
                  correlation matrix, indexed by the keys of the single
                  points in base_grid
  max_corr        time ordered list of GeoPoints containing the maximum correlation 
                  values in the grid and their latitude and longitude
  locations       list of GeoPoints containing possible epicenter locations 
                  the value field contains a (time, corr_value) tuple
                    
    
    
  The dictionary keys are integers.  They are taken from base_grid.  Each
  geographical point in the dictionary has as value a vector containing th
  averaged cross correlation timeseries for that point.  
  """
  
  base_grid=None
  path_length={}
  max_corr=[]
  low_line=None
  high_line=None
  n_timelags=None
  b=None
  dt=None
  ref_time=None
  
  
  def __init__(self,base_grid):
    """
    Initialize the geometry of the grid from that of the base grid.
    """
    self.base_grid=base_grid
    #self.min_lat=base_grid.min_lat
    #self.min_lon=base_grid.min_lon
    #self.max_lat=base_grid.max_lat
    #self.max_lon=base_grid.max_lon

  def t_array(self):
    """
    Returns an array of the times corresponding to the cross-correlation samples.
    """
    # create the time array
    times=arange(self.n_timelags)
    times=times*self.dt + self.b
    return times
    
  def extract_paths_from_matrix(self,sta_grid,corr_matrix):
    """
    Averages the cross-correlation timeseries along the paths corresponding to each
    point.
    
    Parameters:
      sta_grid      a Grid_Sta_GF object
      corr_matrix   a CorrMatrix
    
    The paths for each grid point are taken from sta_grid.corr_keys, which is a 
    dictionary indexed by grid id, whose values are lists of keys of the correlation
    matrix.  
    
    Calls corr_matrix.extract_path_and_sum() to extract the paths given by sta_grid.corr_keys
    for each point.  The path sums (i.e. the summed cross-correlograms along each path)
    are then normalized by the number of steps in each path.  The number of steps per path
    is saved in self.path_length dictionary with the grid id as key.
    
    Each grid point is then assigned the normalized cross-correlogram.  Should the length of 
    the path for a given point be zero, that point is assigned a zero cross-correlogram.   
    """
    
    # do some initialization
    self.b=corr_matrix.b
    self.dt=corr_matrix.dt
    self.ref_time=corr_matrix.ref_time
    self.n_timelags=corr_matrix.n_timelags
    self.path_length={}
    self.points={}
    
    # iterate over points in the base_grid
    for grid_id,point in self.base_grid.points.iteritems():
      
      # extract the number of keys and the summed cross-correlation vector from the
      # cross-correlation matrix using the path indicated by sta_grid.corr_keys for
      # this grid point
      (n_keys,a_vector)=corr_matrix.extract_path_and_sum(sta_grid.corr_keys[grid_id])
      self.path_length[grid_id]=n_keys
      
      # assign the normalized cross-correlogram to this point in the grid, dealing with
      # zero length paths
      try:
        self.points[grid_id]=GeoPointLL(point.lat,point.lon,value=a_vector/n_keys)
      except ZeroDivisionError:
        self.points[grid_id]=GeoPointLL(point.lat,point.lon,value=zeros(self.n_timelags))


  
  def extract_time(self,time):
    """
    Returns the cross-correlation values for the grid for a given time.
    
    Parameters:
    time          time of timestep, defined in seconds wrt to the start
                  of the cross-correlogram (self.b)
                 
    Calls extract_timestep().
    
    The values are returned as a list of tuples:
    [(point_latitude, point_longitude, cross_correlation)]

    """
    t_index=int((time-self.b)/self.dt)
    return self.extract_timestep(t_index)


  def extract_timestep(self,t_index):
    """
    Returns the cross-correlation values for the grid for a given timestep.
    
    Parameters:
    t_index       the ordinal index of the cross_correlogram (the t_index'th
                  timestep will be extracted)
                    
    The values are returned as a list of tuples:
      [(point_latitude, point_longitude, cross_correlation)]
    """
    
    corr_at_timestep=[(p.lat,p.lon,p.value[t_index]) for p in self.points.values()]
    return corr_at_timestep
 
 
    
  def plot_timerange(self,time_step,time_range,nx_steps=100,ny_steps=100,data_list=None,
                     base_filename="timestamp",maxcorr_sweep=False,cbase_filename=""):
    """
    Plots the cross-correlation grids for a given time range.
       
    Parameters:
    time_step     time step within the time range
    time_range    +/- this number of seconds around the maximum
      
    Optional Parameters:
    nx_steps      number of steps in the x directions for grid resampling
                  defaults to 100
    ny_steps      number of steps in the x directions for grid resampling
                  defaults to 100
    data_list     the Data_List object that contains the data used for this 
                  location - will be used to plot station locations; defaults to None
    base_filename the base of the filenames used to store the individual plots for
                  each timestep
    max_corr_timerange  if True, triggers plotting of max_correlation figure with sweeping
                        time indicator - defaults to False
    cbase_filename the base of the filenames used to store the individual plots of max_corr 
                    for each timestep
                    
    As the grid can be irregular, it has to be re-sampled to be plotted, hence the nxy_steps
    parameter. 
    
    The plots are done in geographical projection using the functionalities of matplotlib.basemap.
    
    Calls plot_timestep(). 
    """
    
    # obtain start and end times from times of maxima
    T=[loc.value[0] for loc in self.locations]
    Tmax=max(T)
    Tmin=min(T)
    start_time=max(Tmin-time_range,self.b)
    end_time=min(Tmax+time_range,self.b+self.dt*self.n_timelags)
    
    # turn the time range into a list of times for timestep extraction
    nsteps=int(1+(end_time-start_time)/time_step)
    times=[start_time + i*time_step for i in range(nsteps)]
    
    # use the first time extraction to set the size of the geographical
    # arrays
    start_corr=self.extract_time(start_time)
    
    lats=numpy.array([ t[0] for t in start_corr])
    lons=numpy.array([ t[1] for t in start_corr])
    
    min_lon=min(lons)
    max_lon=max(lons)
    min_lat=min(lats)
    max_lat=max(lats)

    # set up the geographical projection
    m = Basemap(llcrnrlon=min_lon,llcrnrlat=min_lat,urcrnrlon=max_lon,urcrnrlat=max_lat,
            projection='lcc',lat_1=38.,lat_2=45.,lon_0=10.,
            resolution ='i',area_thresh=1000.)

    # project the latitude and logitude arrays
    x, y = m(lons,lats)

    # set up the regridding axes in x,y coordinates
    x_step=(max(x)-min(x))/nx_steps
    y_step=(max(y)-min(y))/ny_steps
    x_i=numpy.arange(min(x),min(x)+nx_steps*x_step,x_step)
    y_i=numpy.arange(min(y),min(y)+ny_steps*y_step,y_step)

    # if the data_list is not None, then set up the latitude and longitude
    # arrays for plotting the stations
    if data_list:
      sta_lat=[]
      sta_lon=[]
      for sta_key in data_list.data.keys():
         sta_lat.append(data_list.cha_list.channels[sta_key].lat)
         sta_lon.append(data_list.cha_list.channels[sta_key].lon)
         #sta_lat.append(data_list.sta_list.stations[sta_key].lat)
         #sta_lon.append(data_list.sta_list.stations[sta_key].lon)
      sta_lat=numpy.array(sta_lat)
      sta_lon=numpy.array(sta_lon)
      # project the station arrays
      xsta,ysta=m(sta_lon,sta_lat)
    else:
      xsta=numpy.empty(1)
      ysta=numpy.empty(1)
    
    # set the contour levels based on the maximum and average of self.max_corr
    A=[p.value for p in self.max_corr]
    cont_max=max(A)
#    cont_min=min(A)
    cont_min=0
    cont_step=(cont_max-cont_min)/20.0
    contours=arange(cont_min,cont_max,cont_step)
    # iterate over times in the time range, set the filename and plot the cross-correlation
    for time in times:
      filename="%s_%07.2f.png"%(base_filename,time)
      self.plot_timestep(m,time,x,y,x_i,y_i,contours,xsta,ysta,filename)
      if maxcorr_sweep:
        filename="%s_max_corr_%07.2f.png"%(cbase_filename,time)
        self.plot_max_corr(title="Maximum correlation",filename=filename,snapshot_time=time)
 
 
    
  def plot_timestep(self,m,time,x,y,x_i,y_i,contours,xsta,ysta,filename=""):
    """
    """
    corr_at_timestep=self.extract_time(time)
    corr=numpy.array([ t[2] for t in corr_at_timestep])
    
 
    corr_i=griddata(x,y,corr,x_i,y_i)
    
    
    pylab.clf()
    date_time=datetime_at_time(self.ref_time,time)
    pylab.title("Timestep at time %04d/%02d/%02d %02d:%02d:%05.2f" % (date_time.year,date_time.month,\
    date_time.day,date_time.hour,date_time.minute,date_time.second))
    m.drawmapboundary(fill_color='#ffffff')
#    m.fillcontinents(color='#cc9966',lake_color='#99ffff')
    cs=m.contour(x_i,y_i,corr_i,contours,linewidths=0.5,colors='k')
    cs=m.contourf(x_i,y_i,corr_i,contours,cmap=pylab.cm.jet)
    m.drawcoastlines()
    m.drawcountries()
    m.drawparallels(arange(35,50,1),labels=[1,1,0,0])
    m.drawmeridians(arange(5,24,1),labels=[0,0,0,1])
    if len(xsta) and len(ysta):
      m.scatter(xsta,ysta, s=50, c='r', marker='^')

    if not filename=="":
      pylab.savefig(filename)
    else:
      pylab.show()
     
  
  def extract_maximum_per_timestep(self):
    self.max_corr=[]
    x=numpy.zeros(self.n_timelags)
    for t_index in range(self.n_timelags): 
      corr_at_timestep=[(p.value[t_index],p.lat,p.lon) for p in self.points.values()] 
      max_corr=max(corr_at_timestep) 
      self.max_corr.append(GeoPointLL(max_corr[1],max_corr[2],value=max_corr[0]))
      x[t_index]=max_corr[0]
    st=Stream()
    st.traces.append(Trace(data=x))
    st.write('test_alberto.sac', format='SAC')

  def write_maxcorr_as_sac_file(self,filename="MaxCorr.sac"):

    import array
    
    # create basic header
    (hf,hi,hs)=pysacio.NewHeader()

    # set some parameters
    pysacio.SetHvalue("b",self.b,hf,hi,hs)
    pysacio.SetHvalue("delta",self.dt,hf,hi,hs)
    pysacio.SetHvalue("npts",len(self.max_corr),hf,hi,hs)
    pysacio.SetHvalue("iftype",pysacio.enum_dict["ITIME"],hf,hi,hs)

    pysacio.SetHvalue("nzyear",self.ref_time.year,hf,hi,hs)
    jday=month_day_to_jday(self.ref_time.year,self.ref_time.month,self.ref_time.day)
    pysacio.SetHvalue("nzjday",jday,hf,hi,hs)
    pysacio.SetHvalue("nzhour",self.ref_time.hour,hf,hi,hs)
    pysacio.SetHvalue("nzmin",self.ref_time.minute,hf,hi,hs)
    pysacio.SetHvalue("nzsec",int(self.ref_time.second),hf,hi,hs)
    msecs=int((self.ref_time.second - int(self.ref_time.second))*1000)
    pysacio.SetHvalue("nzmsec",msecs,hf,hi,hs)

    seis=array.array('f')
    A=[(p.value) for p in self.max_corr]
    seis.extend(A)
    pysacio.WriteSacBinary(filename,hf,hi,hs,seis)


  def write_maxcorr(self,filename="MaxCorr.txt"):

    times=arange(self.n_timelags)
    times=times*self.dt + self.b
    d=self.ref_time

    file=open(filename,'w')
    file.write("# Ref. Time : %04d/%02d/%02d %02d:%02d:%05.2f\n"% \
              (d.year,d.month,d.day,d.hour,d.minute,d.second))
    A=[(p.value,p.lat,p.lon) for p in self.max_corr]
    for i in range(len(A)):
      file.write("%.2f %.2f %.3f %.3f\n"%(times[i],A[i][0],A[i][1],A[i][2]))
    file.close()
      



  def define_maxcorr_thresholds(self,threshold_type,threshold_data):
    #     t_index=int((time-self.b)/self.dt)
    # make list of max_corr values
    av_line=[]
    stdev_line=[]
    
    if threshold_type=='Dynamic' : 
      timescale=threshold_data
    
      A=[(p.value) for p in self.max_corr]
      n_timescale=int(timescale/self.dt)
    
      # if timescale is longer than time series
      if n_timescale >= self.n_timelags:
        av_line=ones(self.n_timelags)*average(A)
        #stdev_line=ones(self.n_timelags)*stats.lstdev(A)
        stdev_line=ones(self.n_timelags)*std(A)      
        self.low_line=array(av_line)
        self.high_line=array(av_line)+array(stdev_line)
      else:
        # the first n_timescale/2 values are constant
        n_by_2=int(n_timescale/2)
        av_line.extend(ones(n_by_2)*average(A[0:n_timescale-1]))
        #stdev_line.extend(ones(n_by_2)*stats.lstdev(A[0:n_timescale-1]))
        stdev_line.extend(ones(n_by_2)*std(A[0:n_timescale-1]))
        # now do the remaining steps
        n_steps=self.n_timelags-n_timescale
        for i in range(n_steps):
          av_line.append(average(A[i:i+n_timescale]))
          #stdev_line.append(stats.lstdev(A[i:i+n_timescale]))
          stdev_line.append(std(A[i:i+n_timescale]))
        # the last n_timescale/2 values are consant
        # how many more points do we need
        points_left=self.n_timelags-len(av_line)
        av_line.extend(ones(points_left)*\
        average(A[self.n_timelags-1-n_timescale:self.n_timelags-1]))
        stdev_line.extend(ones(points_left)*\
        #stats.lstdev(A[self.n_timelags-1-n_timescale:self.n_timelags-1]))
        std(A[self.n_timelags-1-n_timescale:self.n_timelags-1]))
        self.low_line=array(av_line)
        self.high_line=array(av_line)+array(stdev_line)
        
    if threshold_type=='Static' : 
      (high_value,low_value)=threshold_data
      self.high_line=array(ones(self.n_timelags)*high_value)
      self.low_line=array(ones(self.n_timelags)*low_value)
            

  def locate_maxima_flexwin(self,flex_params):
    """
    Locate maxima using flexwin-style reasoning.
    """
    self.locations=[]
    
    (t_kurtosis,w_level,c0,c1,c2,c3a,c3b)=flex_params
    
    # make simple list of values
    A_tmp=[(p.value) for p in self.max_corr]
    
    # max_corr is too rough - smooth it using a time-width of t_kurtosis/4
    t=self.t_array()
    tck = interpolate.splrep(t,A_tmp,s=0.25*t_kurtosis/self.dt)
    A = interpolate.splev(t,tck)
    
    # find local max / min
    z=zeros(len(self.max_corr))
    i_maxima=find_maxima(A,z)
    i_minima=find_minima(A,z)
    
    
    iMLR=setup_MLR(A,i_maxima,i_minima,w_level,c1*t_kurtosis/self.dt,c2*w_level)
    print("After setup : %d maxima\n")%(len(iMLR))
   
    iMLR_w=reject_on_water_level(A,iMLR,i_minima,c0*w_level)
    print("After water level rejection : %d maxima\n")%(len(iMLR_w))
   
    iMLR_c=reject_on_center_maximum(A,iMLR_w,i_maxima)
    print("After center maximum rejection : %d maxima\n")%(len(iMLR_c))
    #print_iMLR(A,iMLR_c,self.b,self.dt)
        
    iMLR_s=reject_on_separation(A,iMLR_c,i_maxima,i_minima,c3a,c3b*t_kurtosis/self.dt)
    print("After separation rejection : %d maxima\n")%(len(iMLR_s))
    #print_iMLR(A,iMLR_s,self.b,self.dt)
    
    iMLR=remove_subwindows(A,iMLR_s,t_kurtosis/self.dt)
    print("After subwindow rejection : %d maxima\n")%(len(iMLR))
    #print_iMLR(A,iMLR,self.b,self.dt)
        
    max_indexes=unique_maxima(A,iMLR,t_kurtosis/self.dt)
    print("After removal of duplicates : %d maxima\n")%(len(max_indexes))
     
    # the indexes of the true maxima are now in max_indexes
    # put the corresponding points (with times) in self.locations 
    for m in max_indexes:
      max_corr=self.max_corr[m]
      t=self.b+m*self.dt
      self.locations.append(GeoPointLL(max_corr.lat, max_corr.lon, value=(t, max_corr.value, datetime_at_time(self.ref_time,t))))
    

  def locate_maxima(self,threshold_type,threshold_data):
    # note: this is tricky
    
    if threshold_type=='Flexwin':
      self.locate_maxima_flexwin(threshold_data)
      return()
      
    self.locations=[]
    
    # initialize thresholds 
    self.define_maxcorr_thresholds(threshold_type,threshold_data)    
    
    # make simple list of values    
    A=[(p.value) for p in self.max_corr]
    
    # set waterlevel
    w_level=self.high_line
        
    # construct list of local maxima > w_level
    indexes_of_maxima=[]
    indexes_of_maxima=find_maxima(A,w_level)
    n_max=len(indexes_of_maxima) 
    
    
    # treat cases n_max=0,1 separately
    if n_max==0:
      return()
    if n_max==1:
      max_index=indexes_of_maxima[0]
      self.locations.append(GeoPointLL(self.max_corr[max_index].lat, self.max_corr[max_index].lon,(self.b+max_index*self.dt,self.max_corr[max_index].value)))
      return()
    
    
    w_level_brackets=self.low_line
    # bracket maxima
    brackets=[]
    for i in indexes_of_maxima:
      brackets.append(bracket_maximum(A,i,w_level_brackets))
    
      
    # for each group of maxima that has the same brackets, keep only the highest one    
    max_indexes=[]
    max_index=indexes_of_maxima[0]
    max_value=A[max_index]
    for ib in range(1,n_max):
      if brackets[ib]==brackets[ib-1]: 
        # we are in the same group
        my_index=indexes_of_maxima[ib]
        my_value=A[my_index]
        if my_value > max_value:
          max_index=my_index
          max_value=my_value
        if ib==n_max-1:
          # we're at the last point, so save this maximum
          max_indexes.append(max_index)
      else: 
        # we move onto the next group, so save this maximum and re-initialize to current index
        max_indexes.append(max_index)
        max_index=indexes_of_maxima[ib]
        max_value=A[max_index]      
            
    # the indexes of the true maxima are now in max_indexes
    # put the corresponding points (with times) in self.locations 
    for m in max_indexes:
      max_corr=self.max_corr[m]
      t=self.b+m*self.dt
      self.locations.append(GeoPointLL(max_corr.lat, max_corr.lon, value=(t, max_corr.value, datetime_at_time(self.ref_time,t))))


  def write_locations_to_file(self,filename):
    """
    Writes locations to a file.
    
    Parameters:
      filename      name of file to write
    
    """
    
    # write the file
    file=open(filename,'w')
    file.write("%d\n" % len(self.locations))
    for loc in self.locations:
      date_time=loc.value[2]
      file.write("%10.5f %10.5f %10.2fs %10.2f %04d/%02d/%02d %02d:%02d:%05.2f\n" % \
      (loc.lat,loc.lon,loc.value[0],loc.value[1], date_time.year,date_time.month,\
      date_time.day,date_time.hour,date_time.minute,date_time.second+date_time.microsecond/1000000.0))   
    file.close()

      

  def plot_max_corr(self,title="Maximum correlation",filename="",snapshot_time=None):
    # clear the figure
    pylab.clf()
    pylab.title(title)
    pylab.xlabel("Time / s")
    pylab.ylabel("Maximum cross-correlation")

    # prepare for the waveform plot
    max_values=[p.value for p in self.max_corr ]
    max_corr=max(max_values)
    min_corr=min(max_values)

    # if a shapshot time is give, then plot it
    if not snapshot_time == None :
      pylab.plot([snapshot_time,snapshot_time],[min_corr,max_corr],'y-',lw=4)
    
    # plot the waveform
    pylab.plot(self.t_array(),max_values)
    
    tmin=self.b
    tmax=self.b+self.dt*self.n_timelags
    # plot low line
    #pylab.plot(self.t_array(),self.low_line,'b-.')
    # plot high line
    #pylab.plot(self.t_array(),self.high_line,'b-')
    
    
    # plot the picked maxima
    for loc in self.locations:
      time=loc.value[0]
      corr=loc.value[1]
      date_time=loc.value[2]
      pylab.plot([time,time],[min_corr,corr],'r-',lw=2)
#      pylab.text(time,corr,"  (%.2fN,%.2fE)\n  %04d/%02d/%02d %02d:%02d:%05.2f" % \
#      (loc.lat,loc.lon,date_time.year,date_time.month,\
#      date_time.day,date_time.hour,date_time.minute,date_time.second), \
#      rotation='horizontal', ha='left', va='top')
      pylab.text(time,corr,"  (%04d/%02d/%02d %02d:%02d:%05.2f)" % \
      (date_time.year,date_time.month,\
      date_time.day,date_time.hour,date_time.minute,date_time.second), \
      rotation='vertical', size=8, ha='center', va='bottom')
    
    # save it to file
    if not filename=="":
      pylab.savefig(filename)
    else:
      pylab.show()
   
  def write_to_file(self,filename_key,filename_value):
    # open the key file in ascii write
    # file format:
    # number of elements, length of arrays, dt, b
    # key_1 lat1 lon1
    # key_2 lat2 lon2
    # key_n
    file_key=open(filename_key,'w')
    file_key.write("%d %d %f %f\n"%(self.npts,self.n_timelags,self.dt,self.b))
    
    file_value=open(filename_value,'wb')
        
    for key,point in self.points.iteritems():
      file_key.write("%d %.4f %.4f\n"%(key,point.lat,point.lon))
      point.value.tofile(file_value)
    
    file_value.close()
    file_key.close()
    
  def read_from_file(self,filename_key,filename_value):
    self.points={}
    # read stuff from the key file
    file_key=open(filename_key,'r')
    lines=file_key.readlines()
    file_key.close()
    
    # extract important numbers
    num_elements=int(lines[0].split()[0])
    self.n_timelags=int(lines[0].split()[1])
    self.dt=float(lines[0].split()[2])
    self.b=float(lines[0].split()[3])
    
    # construct key list
    keys=[]
    lats=[]
    lons=[]
    for line in lines[1:]:
      words=line.split()
      keys.append(int(words[0]))
      lats.append(float(words[1]))
      lons.append(float(words[2]))
    
    # read arrays from file
    file_value=open(filename_value,'rb')
    for key,lat,lon in zip(keys,lats,lons):
      values=fromfile(file_value,dtype='float',count=self.n_timelags)
      self.points[key]=GeoPointLL(lat,lon,value=values)
      
    file_value.close()
         



class Grid_Sta_GF(object):
  sta_id=None
  sta_name=None
  sta_lat=None
  sta_lon=None
  grid_id=None
  grid_lat=None
  grid_lon=None
  dist=None
  gf_id=None
  
  def __init__(self,sta_id,sta_name,sta_lat,sta_lon,grid_id,grid_lat,grid_lon,dist,gf_id):
    self.sta_id=sta_id
    self.sta_name=sta_name
    self.sta_lat=sta_lat
    self.sta_lon=sta_lon
    self.grid_id=grid_id
    self.grid_lat=grid_lat
    self.grid_lon=grid_lon
    self.dist=dist
    self.gf_id=gf_id
   

class Grid_Sta_GF_List(object):
  
  associations=[]
  corr_keys={} # dictionary indexed by grid_id, contains lists of "(gf_id,sta_id)" keys that
                     # will be the keys of the correlation matrix
  
  def __init__(self):
    pass
  
  def create_from_grid_sta_gf(self,grid,cha_list,gf_list,type='Maxdist',args_tuple=None):
    self.associations=[]
    try:
      # iterate through grid points and stations and select the correct gfns
      for grid_id,point in grid.points.iteritems():
        for cha_id,cha  in cha_list.channels.iteritems():
          # call the appropriate selection routine
          ok=False
          if type=='Maxdist':
            (mindist,maxdist,tol)=args_tuple
            (ok,gf_id)=self._select_gf_by_distance_only_(gf_list,point,cha,mindist,maxdist,tol)                              
          else:
            raise UserWarning('Invalid type of selection given to create_from_grid_sta_gf')
            
          # add association of grid point, station, green's function to list
          if ok:
            dist=point.dist_km_from(cha.lat,cha.lon)
            self.associations.append(Grid_Sta_GF(cha_id,cha.name,cha.lat,cha.lon, 
                                          grid_id,point.lat,point.lon,dist,gf_id))
      
    except:
      raise UserWarning('Could not create grid_sta_gf from given grid, cha and gf objects')
      raise
     
   
  def _select_gf_by_distance_only_(self,gf_list,point,sta,mindist,maxdist,tol=0.0):
    dist=point.dist_km_from(sta.lat,sta.lon)
    ok=False
    gf_id=None
    if (dist > mindist-tol) and (dist < maxdist+tol):
      ok=True
      gf_id=gf_list.id_by_closest_distance(dist)
    return (ok,gf_id)
    
  def read_from_file(self,filename):
    self.associations=[]
    file=open(filename,'r')
    lines=file.readlines()
    for line in lines:
      (sta_id,sta_name,grid_id,sta_lon,sta_lat,grid_lon,grid_lat,dist,gf_id)=line.split()
      self.associations.append(Grid_Sta_GF(int(sta_id),sta_name,float(sta_lat),float(sta_lon), 
                           int(grid_id),float(grid_lat),float(grid_lon),float(dist),int(gf_id)))
  
  def write_to_file(self,filename):
    file=open(filename,'w')
    for ass in self.associations:
      file.write("%5d %s %d %10.4f %10.4f %10.4f %10.4f %8.2f %d\n" % 
                (ass.sta_id,ass.sta_name,
                ass.grid_id,
                ass.sta_lon,ass.sta_lat,
                ass.grid_lon,ass.grid_lat,
                ass.dist,ass.gf_id))
    file.close()

  def construct_corr_keys(self):
    self.corr_keys={}
    #extract only relevant information
    A=[(ass.grid_id,ass.sta_id,ass.gf_id) for ass in self.associations]
    A.sort()
    for (grid_id,sta_id,gf_id) in A:
      # construct the key that will be used in the correlation matrix
      # make these keys fixed length, so they can be written to a binary file 
      # and retrieved correctly
      corr_key="(%04d,%04d)"%(gf_id,sta_id)
      if self.corr_keys.has_key(grid_id): # if this grid point is already in the dictionary
        self.corr_keys[grid_id].append(corr_key) # append this key
      else:
        self.corr_keys[grid_id]=[corr_key] # create a list with this key as first element
        

def make_movie(timestamp_directory_name,movie_filename):
    files=glob.glob(timestamp_directory_name + os.sep + "*.png")
    for file in files:
      jpg_name=file + ".jpg"
      os.system("convert -quality 70 %s %s"%(file,jpg_name))
      #os.system("convert -quality 85 %s %s"%(file,jpg_name))
    os.system("convert -delay 15 %s%s%s %s"%(timestamp_directory_name,os.sep,'*.png.jpg',movie_filename))



def migrate_4D_stack(integer_data, delta, search_grid_filename, time_grid):
  # save the list of data keys
  # note : keys of integer data are all included in keys of time_grid, but there may be more times than data
  wf_ids=integer_data.keys()
  time_dict=time_grid.buf[0]
  time_ids=time_dict.keys()
  #logging.debug('Length of time ids %d, %s'%(len(time_ids),time_ids))
  #logging.debug('Length of waveform ids %d, %s'%(len(wf_ids),wf_ids))

  # save the smallest number of points of all the data streams 
  # this will dimension many of the subsequent arrays
  min_npts=min([len(integer_data[key]) for key in wf_ids])
  logging.debug("Stack max time dimension = %d"%min_npts)

  # The stack grid has exactly the same geometry as the time-grid
  #stack_grid=QDStackGrid(time_grid.nx,time_grid.ny,time_grid.nz,min_npts)
  #stack_grid=np.zeros((time_grid.nx,time_grid.ny,time_grid.nz,min_npts),dtype=np.int32)
  #stack_grid=np.zeros((time_grid.nx,time_grid.ny,time_grid.nz,min_npts),dtype=np.float)
  stack_grid=np.zeros((time_grid.nx,time_grid.ny,time_grid.nz,min_npts))
  #stack_grid.read_NLL_hdr_file(search_grid_filename)
  #stack_grid.construct_empty_grid(min_npts)

  # Number of geographical points in the stack
  n_buf=time_grid.nx*time_grid.ny*time_grid.nz
  
  # keep information on the shortest length of stack for later
  shortest_n_len=min_npts

  # set up the stack grid 
  for ib in range(n_buf):

      times=time_grid.buf[ib]
      ix,iy,iz=time_grid.get_ix_iy_iz(ib)


      # find the slice indexes
      i_times=[int(round(times[wf_id]/delta)) for wf_id in wf_ids]
      min_i_time=min(i_times)
      max_i_time=max(i_times)
      start_end_indexes=[(i_time-min_i_time, i_time+min_npts-max_i_time) for i_time in i_times]
      n_lens=[start_end_indexes[i][1]-start_end_indexes[i][0] for i in range(len(wf_ids))]
      n_len=min(n_lens)

      # keep shortest n_len for later
      if n_len < shortest_n_len:
        shortest_n_len=n_len

      # initialize the stack
      #stack=numpy.zeros(min_npts,dtype=np.int32)
      #stack=numpy.zeros(min_npts,dtype=np.float)
      stack=numpy.zeros(min_npts)

      for i in range(len(wf_ids)):
        wf_id=wf_ids[i]
        #stack[0:n_len] += integer_data[wf_id][start_end_indexes[i][0]:start_end_indexes[i][1]]
        stack[0:n_lens[i]] += integer_data[wf_id][start_end_indexes[i][0]:start_end_indexes[i][1]]

      #stack_grid[ix,iy,iz,0:n_len] = stack[0:n_len]
      stack_grid[ix,iy,iz,:] = stack[:]
    
      
  logging.debug('Stacking done.')

######## FIXUP THE CORR GRID START TIMES #########


  # Each stack starts at ref_time - the minimum travel-time and ends at the ref_time + seismogram duration - max travel_time
  # We need to homogenize, and get everything to start and end at the same time


  # deal with the start of the traces
  # start index for slice = min_itime for the single stack - smallest min_itime for all stacks

  logging.debug('Fixing up stack start times')
#  iextreme_min_times=[int(round(min(time_grid.buf[ib])/delta)) for ib in range(time_grid.npts)]
#  iextreme_max_times=[int(round(max(time_grid.buf[ib])/delta)) for ib in range(time_grid.npts)]
  iextreme_min_times=[int(round(min([time_grid.buf[ib][wf_id] for wf_id in wf_ids])/delta))  for ib in range(n_buf) ]
  iextreme_max_times=[int(round(max([time_grid.buf[ib][wf_id] for wf_id in wf_ids])/delta))  for ib in range(n_buf) ]
  iextreme_min_time=min(iextreme_min_times)
  iextreme_max_time=max(iextreme_max_times)

  # fix the length of the stack to the shortest possible length given all the previous travel time information
  norm_stack_len=shortest_n_len-iextreme_max_time

  # iterate over the time-arrays in the time_grid to extract the minimum and fix up the stacks
  for ib in range(n_buf):
    ix,iy,iz=time_grid.get_ix_iy_iz(ib)
    start_index = iextreme_min_times[ib] - iextreme_min_time
#    tmp=stack_grid.buf[ib][:]
    tmp=stack_grid[ix,iy,iz,:]
    try:
      #stack_grid.buf[ib][0:norm_stack_len]=tmp[start_index:start_index+norm_stack_len]
      stack_grid[ix,iy,iz,0:norm_stack_len]=tmp[start_index:start_index+norm_stack_len]
    except ValueError:
#      logging.debug('(norm_stack_len,shortest_n_len,iextreme_max_time) = (%s,%s,%s)'%(norm_stack_len,shortest_n_len,iextreme_max_time))
#      logging.debug("(ib,norm_stack_len,start_index) = (%s,%s,%s)"%(ib,norm_stack_len,start_index))
      logging.error("Length of time slice for migration too short compared with the largest migration time.")
      raise 
   
  logging.debug('Done fixing up stack start times')
  stack_shift_time=delta*iextreme_min_time
  return n_buf, norm_stack_len, stack_shift_time, stack_grid
 
def migrate_3D_stack(integer_data, delta, search_grid_filename, time_grid):
  # save the list of data keys
  # note : keys of integer data are all included in keys of time_grid, but there may be more times than data
  wf_ids=integer_data.keys()
  time_dict=time_grid.buf[0]
  time_ids=time_dict.keys()
  #logging.debug('Length of time ids %d, %s'%(len(time_ids),time_ids))
  #logging.debug('Length of waveform ids %d, %s'%(len(wf_ids),wf_ids))

  # save the smallest number of points of all the data streams 
  # this will dimension many of the subsequent arrays
  min_npts=min([len(integer_data[key]) for key in wf_ids])
  logging.debug("Stack max time dimension = %d"%min_npts)

  # The stack grid has exactly the same geometry as the time-grid
  stack_grid=QDCorrGrid()
  stack_grid.read_NLL_hdr_file(search_grid_filename)
  stack_grid.construct_empty_grid(min_npts)

  # Number of geographical points in the stack
  n_buf=stack_grid.nx*stack_grid.ny*stack_grid.nz
  
  # keep information on the shortest length of stack for later
  shortest_n_len=min_npts

  # set up the stack grid 
  for ib in range(n_buf):

      times=time_grid.buf[ib]

      # find the slice indexes
      i_times=[int(round(times[wf_id]/delta)) for wf_id in wf_ids]
      min_i_time=min(i_times)
      max_i_time=max(i_times)
      start_end_indexes=[(i_time-min_i_time, i_time+min_npts-max_i_time) for i_time in i_times]
      n_lens=[start_end_indexes[i][1]-start_end_indexes[i][0] for i in range(len(wf_ids))]
      n_len=min(n_lens)

      # keep shortest n_len for later
      if n_len < shortest_n_len:
        shortest_n_len=n_len

      # initialize the stack
      stack=numpy.zeros(min_npts,dtype=numpy.int16)

      for i in range(len(wf_ids)):
        wf_id=wf_ids[i]
        stack[0:n_len] += integer_data[wf_id][start_end_indexes[i][0]:start_end_indexes[i][1]]
    
      stack_grid.buf[ib][0:n_len]=stack[0:n_len]
      
  logging.debug('Stacking done..')

######## FIXUP THE CORR GRID START TIMES #########


  # Each stack starts at ref_time - the minimum travel-time and ends at the ref_time + seismogram duration - max travel_time
  # We need to homogenize, and get everything to start and end at the same time


  # deal with the start of the traces
  # start index for slice = min_itime for the single stack - smallest min_itime for all stacks

  logging.debug('Fixing up stack start times')
#  iextreme_min_times=[int(round(min(time_grid.buf[ib])/delta)) for ib in range(time_grid.npts)]
#  iextreme_max_times=[int(round(max(time_grid.buf[ib])/delta)) for ib in range(time_grid.npts)]
  iextreme_min_times=[int(round(min([time_grid.buf[ib][wf_id] for wf_id in wf_ids])/delta))  for ib in range(n_buf) ]
  iextreme_max_times=[int(round(max([time_grid.buf[ib][wf_id] for wf_id in wf_ids])/delta))  for ib in range(n_buf) ]
  iextreme_min_time=min(iextreme_min_times)
  iextreme_max_time=max(iextreme_max_times)

  # fix the length of the stack to the shortest possible length given all the previous travel time information
  norm_stack_len=shortest_n_len-iextreme_max_time

  # iterate over the time-arrays in the time_grid to extract the minimum and fix up the stacks
  for ib in range(n_buf):
    start_index = iextreme_min_times[ib] - iextreme_min_time
    tmp=stack_grid.buf[ib][:]
    try:
      stack_grid.buf[ib][0:norm_stack_len]=tmp[start_index:start_index+norm_stack_len]
    except ValueError:
#      logging.debug('(norm_stack_len,shortest_n_len,iextreme_max_time) = (%s,%s,%s)'%(norm_stack_len,shortest_n_len,iextreme_max_time))
#      logging.debug("(ib,norm_stack_len,start_index) = (%s,%s,%s)"%(ib,norm_stack_len,start_index))
      logging.error("Length of time slice for migration too short compared with the largest migration time.")
      raise 
   
  logging.debug('Done fixing up stack start times')
  stack_shift_time=delta*iextreme_min_time
  return n_buf, norm_stack_len, stack_shift_time, stack_grid
 

if __name__ == '__main__':

  from mpl_toolkits.basemap.pyproj import Proj
  from mpl_toolkits.basemap.pyproj import Geod
  logging.basicConfig(level=logging.DEBUG)

  p=Proj(init='epsg:2975') # For La Réunion
  #piton_stations_filename='../../M2-Piton/coord_stations_piton'

  #sta=StationList(p)
  #sta.read_from_file(piton_stations_filename)
  #sta.display()

  print ('**')

  #sta_filtered=sta.filter_by_station_names(['UV01','UV15','TEST'])
  #sta_filtered.display()
  

  grd=QDGrid()
  #grd.read_from_NLL_files('/Users/alessia/grilles_temps/Slow_len.100m.P.FLR.time')
  grd.read_from_NLL_files('../lib/belgium2D.P.MASA.time')
  print grd.npts 
  print grd.min_x, grd.max_x 
  print grd.min_y, grd.max_y
  print grd.min_z, grd.max_z
  print grd.min_value, grd.max_value
  test_x=(grd.min_x+grd.max_x)/2
  test_x=grd.min_x
  test_y=(grd.min_y+grd.max_y)/2
  test_z=(grd.min_z+grd.max_z)/2
  test_value=grd.value_at_point(test_x,test_y,test_z)
  print test_x,test_y, test_z, test_value 
  #grd.display_grid_zcut('testz.pdf',z=test_z)
  #grd.display_grid_ycut('testy.pdf',y=test_y)
  grd.display_grid_xcut('MASA.pdf',x=test_x)
  #grd.display_NLL_XYZ(0.0,filename='test.png')
