#!/usr/bin/env python



# Standard library imports
import os, glob, numpy

# Enthought library imports
#from mayavi.scripts import mayavi2
from enthought.tvtk.api import tvtk
from enthought.mayavi.api import OffScreenEngine, Engine
from enthought.mayavi.sources.array_source import ArraySource
from enthought.mayavi.sources.vtk_data_source import VTKDataSource
from enthought.mayavi.modules.api import Outline, ScalarCutPlane, Surface, Axes, OrientationAxes, ImagePlaneWidget, ContourGridPlane, IsoSurface, Glyph
from enthought.mayavi.filters.api import TransformData
from enthought.mayavi.tools.camera import view

#from mayavi.modules.image_plane_widget import ImagePlaneWidget

#import my stuff
from grids_paths import QDGrid, StationList


def plot_slice_mayavi(dat_filename,output_file,hyp_x,hyp_y,hyp_z,search_grid_file_name,max_stack_value):

  base_path=os.getenv('WAVELOC_PATH')
  lib_path="%s/lib"%base_path

  # grid geometry
  hdr_file=lib_path + os.sep + search_grid_file_name

  # detection
  detection=50

  # stations
  stations_file="%s/coord_stations_piton"%lib_path
  sta=StationList()
  sta.read_from_file(stations_file)

  # create the object to contain the stations
  pd = tvtk.PolyData()
  pd.points = [[s.x/1000.0, s.y/1000.0, -s.elev/1000.0] for s in sta.stations.values()]

  # create the object to contain the stations
  try:
    pd_hyp = tvtk.PolyData()
    pd_hyp.points=[[hyp_x,hyp_y,hyp_z]]
  except TypeError:
    pass

  # read the dat file
  print dat_filename
  data=QDGrid()
  data.read_NLL_hdr_file(hdr_file)
  data.buf=numpy.fromfile(dat_filename, dtype=numpy.int16)
  max_ib=numpy.argmax(data.buf)
  print max_ib
  max_val=data.buf[max_ib]
  ix,iy,iz=data.get_ix_iy_iz(max_ib)
  #data.buf=numpy.array(data.buf, dtype=numpy.float)
  data.buf.shape = (data.nx, data.ny, data.nz)


  # 'mayavi' is always defined on the interpreter.
  e = OffScreenEngine()
  #e = Engine()
  e.start()
  win = e.new_scene(magnification=1)
  e.current_scene.scene.off_screen_rendering = True
  win.scene.isometric_view()
  
  # Make the data and add it to the pipeline.
  src = ArraySource(transpose_input_array=True)
  src.scalar_data = data.buf
  src.spacing=(data.dx, data.dy, -data.dz)
  src.origin=(data.x_orig, data.y_orig, -data.z_orig)
  e.add_source(src)

  # Visualize the data.
  o = Outline()
  e.add_module(o)

  lut=e.scenes[0].children[0].children[0].scalar_lut_manager
  lut.data_range=[-1,max_stack_value]
  lut.show_legend = True
  lut.data_name = 'Stack'


  # Create one ContourGridPlane normal to the 'x' axis.
  cgp = ContourGridPlane()
  e.add_module(cgp)
  # Set the position to the middle of the data.
  if max_val > detection:
    cgp.grid_plane.position = ix
  else:
    cgp.grid_plane.position = data.nx/2
  cgp.contour.filled_contours = True
  cgp.actor.property.opacity = 0.6
  output=cgp.grid_plane.outputs[0]
  x_data=numpy.array(output.point_data.scalars.to_array())

  # Another with filled contours normal to 'y' axis.
  cgp = ContourGridPlane()
  e.add_module(cgp)
  # Set the axis and position to the middle of the data.
  cgp.grid_plane.axis = 'y'
  if max_val > detection:
    cgp.grid_plane.position = iy
  else:
    cgp.grid_plane.position = data.ny/2
  cgp.contour.filled_contours = True
  cgp.actor.property.opacity = 0.6
  output=cgp.grid_plane.outputs[0]
  y_data=numpy.array(output.point_data.scalars.to_array())

  # Another with filled contours normal to 'z' axis.
  cgp = ContourGridPlane()
  e.add_module(cgp)
  # Set the axis and position to the middle of the data.
  cgp.grid_plane.axis = 'z'
  if max_val > detection:
    cgp.grid_plane.position = iz
  else:
    cgp.grid_plane.position = data.nz/2
  cgp.contour.filled_contours = True
  cgp.actor.property.opacity = 0.6
  output=cgp.grid_plane.outputs[0]
  z_data=numpy.array(output.point_data.scalars.to_array())

  a=Axes()
  e.add_module(a)

  d=VTKDataSource()
  d.data=pd
  e.add_source(d)
  
  g=Glyph()
  e.add_module(g)
  g.glyph.glyph_source.glyph_source=g.glyph.glyph_source.glyph_list[4]
  g.glyph.glyph_source.glyph_source.radius=0.1

  d=VTKDataSource()
  d.data=pd_hyp
  e.add_source(d)
  
  g=Glyph()
  e.add_module(g)
  g.glyph.glyph_source.glyph_source=g.glyph.glyph_source.glyph_list[4]
  g.glyph.glyph_source.glyph_source.radius=0.5
  g.actor.property.color=(0.0,0.0,0.0)



  #view(azimuth=-60,elevation=60,distance=120)
  win.scene.save(output_file,size=(800,800))

  e.stop()
  del win
  del e

  return (x_data, y_data, z_data)
