#!/usr/bin/env python

# Standard library imports
import os, glob, numpy
import time

# Enthought library imports
#from mayavi.scripts import mayavi2
from enthought.tvtk.api import tvtk
from enthought.mayavi.api import OffScreenEngine, Engine
from enthought.mayavi.sources.array_source import ArraySource
from enthought.mayavi.sources.vtk_data_source import VTKDataSource
from enthought.mayavi.modules.api import Outline, ScalarCutPlane, Surface, Axes, OrientationAxes, ImagePlaneWidget, ContourGridPlane, IsoSurface, Glyph
from enthought.mayavi.filters.api import TransformData, Delaunay2D
from enthought.mayavi.tools.camera import view

from scipy.io.matlab import loadmat

#from mayavi.modules.image_plane_widget import ImagePlaneWidget

#import my stuff
from grids_paths import QDGrid, StationList


base_path=os.getenv('WAVELOC_PATH')
data_path="%s/out/testing"%base_path
lib_path="%s/lib"%base_path


# stations
stations_file="%s/coord_stations_piton"%lib_path
sta=StationList()
sta.read_from_file(stations_file)

# DEM
dem=loadmat("%s/MNT_PdF.mat"%lib_path)
#print dem
demx=numpy.array(dem['XIsub']).flatten()/1000.0
demy=numpy.array(dem['YIsub']).flatten()/1000.0
demz=numpy.array(dem['ZIsub']).flatten()/1000.0

print "Read DEM"

# data 
data_glob = "*2010-10-14T00:15:40.98*.dat"
data_files=glob.glob(data_path + os.sep + data_glob)

hdr_file="%s/grid.500m.search.hdr"%lib_path

# creat the object to contain the stations
pd = tvtk.PolyData()
pd.points = [[s.x/1000.0, s.y/1000.0, -s.elev/1000.0] for s in sta.stations.values()]

# create the DEM
dem_data=tvtk.PolyData()
dem_data.points = numpy.array([demx, demy, demz]).T

for data_file in data_files : 
  print data_file
  data=QDGrid()
  data.read_NLL_hdr_file(hdr_file)
  data.buf=numpy.fromfile(data_file, dtype=numpy.int16)
  data.buf=numpy.array(data.buf, dtype=numpy.float)
  print data.buf.min(), data.buf.max()
  data.buf.shape = (data.nx, data.ny, data.nz)
  max_ib=numpy.argmax(data.buf)
  ix,iy,iz=data.get_ix_iy_iz(max_ib)


  

  # 'mayavi' is always defined on the interpreter.
  e = OffScreenEngine()
  #e = Engine()
  e.start()
  win = e.new_scene(magnification=1)
  win.scene.isometric_view()
  
  # Make the data and add it to the pipeline.
  src = ArraySource(transpose_input_array=True)
  src.scalar_data = data.buf
  src.spacing=(data.dx, data.dy, -data.dz)
  src.origin=(data.x_orig, data.y_orig, -data.z_orig)
  #print src.all_trait_names()
#  e.add_source(src)

  # Visualize the data.
#  o = Outline()
#  e.add_module(o)

#  lut=e.scenes[0].children[0].children[0].scalar_lut_manager
#  lut.data_range=[-1,600]
#  lut.show_legend = True
#  lut.data_name = 'Stack'


  # Create one ContourGridPlane normal to the 'x' axis.
  cgp = ContourGridPlane()
#  e.add_module(cgp)
  # Set the position to the middle of the data.
  cgp.grid_plane.position = ix
  cgp.contour.filled_contours = True
  cgp.actor.property.opacity = 0.6

  # Another with filled contours normal to 'y' axis.
  cgp = ContourGridPlane()
#  e.add_module(cgp)
  # Set the axis and position to the middle of the data.
  cgp.grid_plane.axis = 'y'
  cgp.grid_plane.position = iy
  cgp.contour.filled_contours = True
  cgp.actor.property.opacity = 0.6

  # Another with filled contours normal to 'z' axis.
  cgp = ContourGridPlane()
#  e.add_module(cgp)
  # Set the axis and position to the middle of the data.
  cgp.grid_plane.axis = 'z'
  cgp.grid_plane.position = iz
  cgp.contour.filled_contours = True
  cgp.actor.property.opacity = 0.6

#  a=Axes()
#  e.add_module(a)

# STATIONS

  #d=VTKDataSource()
  #d.data=pd
  #e.add_source(d)

  #g=Glyph()
  #e.add_module(g)
  #g.glyph.glyph_source.glyph_source=g.glyph.glyph_source.glyph_list[4]
  #g.glyph.glyph_source.glyph_source.radius=0.1

# DEM

  d=VTKDataSource()
  d.data=dem_data
  e.add_source(d)

  df=Delaunay2D()
  e.add_filter(df)

  s=Surface()
  e.add_module(s)
  #s.actor.property.set(representation='p', point_size=100)
  

  #view(azimuth=-60,elevation=60,distance=120)
  output_file="tmp.png"
  win.scene.save(output_file,size=(800,600))



#  s = Surface()
#  e.add_module(s)
#  s.actor.property.set(representation='p', point_size=2)
