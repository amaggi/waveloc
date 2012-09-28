"""
Classes pertaining to NLL type files

Created by Alessia Maggi and Alberto Michelini.  

"""
import numpy as np
from obspy.core import utcdatetime

cPI=np.pi		# PI
cRPD = cPI / 180.0	# radians per degree
c111 = 10000.0/90.0	# km per degree

def read_stations_file(filename):

  stations={}

  # open and read the file
  f=open(filename,'r')
  lines=f.readlines()
  f.close()

  for line in lines:
    sta={}
    words=line.split()
    sta['station']=words[1]
    sta['loc_type']=words[2]
    if sta['loc_type']=='XYZ':
      sta['x']=np.float(words[3])
      sta['y']=np.float(words[4])
    elif sta['loc_type']=='LATLON':
      sta['lat']=np.float(words[3])
      sta['lon']=np.float(words[4])
    else : raise UserWarning('Unknown loc_type %s in file %s'%(words[2],filename))
    sta['depth']=np.float(words[5])
    sta['elev']=np.float(words[6])
    stations[sta['station']]=sta

  return stations

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

    if len(line.split())==4:
        info['station']=line.split()[0]
        info['sta_x']=float(line.split()[1])
        info['sta_y']=float(line.split()[2])
        info['sta_z']=float(line.split()[3])

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




def qd_read_hyp_file(filename):
  f=open(filename, 'r')
  lines=f.readlines()
  f.close()
  
  for line in lines:
    words=line.split()
    try:
      if words[0]=='HYPOCENTER':
        hypo_x=np.float(words[2])
        hypo_y=np.float(words[4])
        hypo_z=np.float(words[6])
      if words[0]=='GEOGRAPHIC':
        year=np.int(words[2])
        month=np.int(words[3])
        day=np.int(words[4])
        hour=np.int(words[5])
        minute=np.int(words[6])
        seconds=np.float(words[7])
        otime=utcdatetime.UTCDateTime(year,month,day,hour,minute,seconds)
      if words[0]=='STATISTICS':
        sigma_x=np.sqrt(np.float(words[8]))
        sigma_y=np.sqrt(np.float(words[14]))
        sigma_z=np.sqrt(np.float(words[18]))
    except IndexError:
      pass

  return (otime, hypo_x, sigma_x, hypo_y, sigma_y, hypo_z, sigma_z)

def qd_read_picks_from_hyp_file(filename):
  f=open(filename, 'r')
  lines=f.readlines()
  f.close()
  
  
  for iline in range(len(lines)):
    line=lines[iline]
    words=line.split()
    if words[0]=='PHASE':
      iline_phase=iline
      break
  
  phases={}
  for line in lines[iline+1:]:
    words=line.split()
    try:
      if words[4]=='P':
        station=words[0]
        year=np.int(words[6][0:4])
        month=np.int(words[6][4:6])
        day=np.int(words[6][6:8])
        hour=np.int(words[7][0:2])
        minute=np.int(words[7][2:4])
        seconds=np.float(words[8])
        ptime=utcdatetime.UTCDateTime(year,month,day,hour,minute,seconds)
        phases[station]=ptime
    except IndexError:
      pass
     
  return phases

