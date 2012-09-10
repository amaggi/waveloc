import numpy as np

cPI=np.pi		# PI
cRPD = cPI / 180.0	# radians per degree
c111 = 10000.0/90.0	# km per degree

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
