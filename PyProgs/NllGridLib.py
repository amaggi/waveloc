import numpy as np

cPI=np.pi		# PI
cRPD = cPI / 180.0	# radians per degree
c111 = 10000.0/90.0	# km per degree

def read_hdr_file(filename):

  # read header file
  f=open(filename)
  lines=f.readlines()
  f.close()
  
  info={}

  # extract information
  vals=lines[0].split()
  info['nx']=int(vals[0])
  info['ny']=int(vals[1])
  info['nz']=int(vals[2])
  info['x_orig']=float(vals[3])
  info['y_orig']=float(vals[4])
  info['z_orig']=float(vals[5])
  info['dx']=float(vals[6])
  info['dy']=float(vals[7])
  info['dz']=float(vals[8])

  for line in lines:
    if line.split()[0]=='TRANSFORM':
      if line.split()[1]=='NONE': info['proj_name']='TRANS_NONE'
      if line.split()[1]=='SIMPLE': 
        info['proj_name']='TRANS_SIMPLE'
        info['orig_lat']=float(line.split()[3])
        info['orig_lon']=float(line.split()[5])
        info['map_rot']=float(line.split()[7])

  return info

def latlon2rect(proj_name,lat,lon,proj_info={}):

  try:
    if proj_name=='TRANS_GLOBAL' or proj_name=='TRANS_NONE':
      x = lon
      y = lat

    if proj_name=='TRANS_SIMPLE':
      xtemp = lon - proj_info['orig_lon']
      if (xtemp > 180.0)  : xtemp -= 360.0
      if (xtemp < -180.0) : xtemp += 360.0
      xtemp = xtemp * c111 * np.cos(cRPD * lat)
      ytemp = (lat - proj_info['orig_lat']) * c111
      
      angle = -cRPD * proj_info['map_rot']
      x = xtemp * np.cos(angle) - ytemp * np.sin(angle)
      y = ytemp * np.cos(angle) + xtemp * np.sin(angle)


    return x,y

  except NameError:
    raise UserWarning('Unknown projection name %s'%proj_name)


def rect2latlon(proj_name,x,y,proj_info={}):

  try:
    if proj_name=='TRANS_GLOBAL' or proj_name=='TRANS_NONE':
      lon = x
      lat = y


    if proj_name=='TRANS_SIMPLE':

      angle = -cRPD * proj_info['map_rot']
      xtemp = x * np.cos(angle) + y * np.sin(angle)
      ytemp = y * np.cos(angle) - x * np.sin(angle)
      lat = proj_info['orig_lat'] + ytemp / c111
      lon = proj_info['orig_lon'] + xtemp / (c111 * np.cos(cRPD * lat))
     

    return lat,lon

  except NameError:
    raise UserWarning('Unknown projection name %s'%proj_name)
