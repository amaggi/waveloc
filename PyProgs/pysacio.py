"""A SAC file access module for Python
\nVersion 2.0.3, by C.J. Ammon, Penn State
This software is free and is distributed with no guarantees.
For a more complete description, start python and enter,

"import pysacio"
"help(pysacio)"

Suspected limitations: I don't used XY files much - I am
not sure that those access routines are bug free.

Send bug reports (not enhancement/feature requests) to: 
cja12@psu.edu [with PySAC in the subject field]
I don't support this software so don't wait for an answer.
I may not have time...

AM (2007) : Have modified the following routines to read the integer 
header fields correctly on both 32 and 64 bit systems.
- ReadSacHeader() 
- ReadSacFile()
- ReadXYSacFile()
Contact me at alessia@sismo.u-strasbg.fr if this breaks.
"""
#
###############################################################################
# This is a set of python functions for reading and writing
#   sac files and accessing their header values 
#
# Version: 2.0.1, 2001 - 2004
#  Author: Charles J. Ammon, Penn State
#
# You would call these functions from a python script.
# To use them, place the file pysacio.py in your python path
# and include the command 'import pysacio' at the top of your
# script. Then prepend the name of the macro package (pysacio)
# on the routines. For example, to use "ReadSacFile", you 
# would call "pysacio.ReadSacFile".
#
###############################################################################
#
import struct,array,os,string,AM_subs,platform,numpy
#
###############################################################################
#
# These dictionaries are for SAC header access
#   use the function GetHvalue('delta',hf,hi,hs)
#   it will return a float, integer, or string
#
fdict = {'delta':0, 'depmin':1, 'depmax':2, 'scale':3,   \
   'odelta':4, 'b':5, 'e':6, 'o':7, 'a':8, 't0':10,\
   't1':11,'t2':12,'t3':13,'t4':14,'t5':15,'t6':16,\
   't7':17,'t8':18,'t9':19,'f':20,'stla':31,'stlo':32,    \
   'stel':33,'stdp':34,'evla':35,'evlo':36,'evdp':38,'mag':39, \
   'user0':40,'user1':41,'user2':42,'user3':43,\
   'user4':44,'user5':45,'user6':46,'user7':47,\
   'user8':48,'user9':49,'dist':50,'az':51,'baz':52,\
   'gcarc':53,'depmen':56,'cmpaz':57,'cmpinc':58}
#
idict = {'nzyear':0, 'nzjday':1, 'nzhour':2, 'nzmin':3, \
   'nzsec':4, 'nzmsec':5, 'nvhdr':6, 'norid':7, \
   'nevid':8,'npts':9, 'nwfid':11, \
   'iftype':15,'idep':16,'iztype':17,'iinst':19,\
   'istreg':20,'ievreg':21,'ievtype':22,'iqual':23,\
   'isynth':24,'imagtyp':25,'imagsrc':26, \
   'leven':35,'lpspol':36,'lovrok':37,\
   'lcalda':38}
#
sdict = {'kstnm':0,'kevnm':1,'khole':2, 'ko':3,'ka':4,\
    'kt0':5,'kt1':6,'kt2':7,'kt3':8,'kt4':9,\
    'kt5':10,'kt6':11,'kt7':12,'kt8':13,\
    'kt9':14,'kf':15,'kuser0':16,'kuser1':17,\
    'kuser2':18,'kcmpnm':19,'knetwk':20,\
    'kdatrd':21,'kinst':22}

# not finished
enum_dict={'IREAL':0,'ITIME':1,'IRLIM':2,'IAMPH':3,'IXY':4}

TRUE=1
FALSE=0

#
###############################################################################
#
def PysacioVersion():
  """Return the module version as a string."""
  myVersion = "2.0.2"
  return(myVersion)
#
###############################################################################
#
# create an empty header
def NewHeader():

  # create empty headers
  hf = array.array('f') # allocate the array for header floats
  # allocate the array for header ints
  (bits,linkage)= platform.architecture()
  if bits == '64bit':
    hi = array.array('i') # this if 4 bytes on by debian 64bit system
  else :
    hi = array.array('l') # allocate the array for header ints
  hs = array.array('c') # allocate the array for header characters
  
  # fill with null values
  for i in range(70):
    hf.append(-12345.0)
  for i in range(40):
    hi.append(-12345)
  for i in range(len(sdict)+1): # note, need to adjust for second string value being 16 chars
    hs.extend(['-','1','2','3','4','5',' ',' '])
  for i in range(16,23):
    hs[i]=' '

  # set the basic stuff only
  SetHvalue('nvhdr',6,hf,hi,hs)
  SetHvalue('ninf',0,hf,hi,hs)
  SetHvalue('nhst',0,hf,hi,hs)
  SetHvalue('lpspol',0,hf,hi,hs)
  SetHvalue('lcalda',1,hf,hi,hs)
  SetHvalue('unused27',0,hf,hi,hs)

  return [hf,hi,hs]


###############################################################################
#
#   access a sac file header values using a string
#    the type returned depends on the what you ask for
#
#   usage: dt = GetHvalue('delta',hf,hi,hs)
#
#     returns the string 'NULL' if nothing matches the item
#
#
def GetHvalue(item,hf,hi,hs):
  """Get a header value using the header arrays: GetHvalue("npts",hf,hi,hs)
  Return value is 1 if no problems occurred, zero otherwise."""
  #
  # it's trivial to search each dictionary with the key and return
  #   the value that matches the key
  #
  key = string.lower(item) # convert the item to lower case
  #
  if fdict.has_key(key):
    index = fdict[key]
    return(hf[index])
  elif idict.has_key(key):
    index = idict[key]
    return(hi[index])
  elif sdict.has_key(key):
    index = sdict[key]
    length = 8
    #
    if index == 0:
      myarray = hs[0:8]
    elif index == 1:
      myarray = hs[8:24]
    else:
      start = 8 + index*8  # the extra 8 is from item #2
      end   = start + 8
      myarray = hs[start:end]
    #
    return(myarray.tostring())
  else:
    return('NULL')
#
###############################################################################
#
#   alter a sac file header values using a string and a value
#    
#
#   usage: SetHvalue('delta',1.0, hf,hi,hs)
#
#     sets the dt value to 1.0 
#
#
def SetHvalue(item,value,hf,hi,hs):
  """Set a header value using the header arrays: SetHvalue("npts",2048,hf,hi,hs)
  Return value is 1 if no problems occurred, zero otherwise."""
  #
  # it's trivial to search each dictionary with the key and return
  #   the value that matches the key
  #
  key = string.lower(item) # convert the item to lower case
  #
  ok = 0
  if fdict.has_key(key):
    index = fdict[key]
    hf[index] = float(value)
    ok = 1
  elif idict.has_key(key):
    index = idict[key]
    hi[index] = int(value)
    ok = 1
  elif sdict.has_key(key):
    index = sdict[key]
    vlen = len(value)
    if index == 0:
      if vlen > 8:
        vlen = 8
      for i in range(0,8):
        hs[i] = ' '
      for i in range(0,vlen):
        hs[i] = value[i]
    elif index == 1:
      start = 8
      if vlen > 16:
        vlen =16 
      for i in range(0,16):
        hs[i+start] = ' '
      for i in range(0,vlen):
        hs[i+start] = value[i]
    else:
      #
      # if you are here, then the index > 2
      #
      if vlen > 8:
        vlen = 8
      start  = 8 + index*8 
      for i in range(0,8):
        hs[i+start] = ' '
      for i in range(0,vlen):
        hs[i+start] = value[i]
    ok = 1
  
  return(ok)
#
#
###############################################################################
#
# Check for a valid SAC file (returns 1 if SAC File, 0 if not a SAC File)
#
#  Right now the checks are very basic (what else can you do?)
#    I check for a positive dt and npts, i could use the version, but
#    this isn't even listed in the ASCII definition of the header...
#    implements a file-size check where npts*4 + headerbytes (632)
#    is compared with the true size of the file.
#
#  usage:   ok = IsSACfile(name,hf,hi,hs)
#
#           if ok = 1, it is a SAC file, if ok = 0, it's not
#   
def IsSACfile(name,hf,hi,hs):
  """Test for a valid SAC file using arrays: IsSACfile(path,hf,hi,hs)
  Return value is a one if valid, zero if not."""
  #
  ok = 1
  #
  # size check info
  #
  npts = GetHvalue('npts',hf,hi,hs)
  st = os.stat(name) #file's size = st[6] 
  sizecheck = st[6] - (632 + 4 * npts)

  #
  # get the SAC file version number
  #
  version = GetHvalue('nvhdr',hf,hi,hs)
  #
  # if any of these conditions are true,
  #   the file is NOT a SAC file
  #
  if GetHvalue('delta',hf,hi,hs) <= 0:
    ok = 0
  elif sizecheck != 0:
    ok = 0
  elif ((version < 0) or (version > 20) ):
    ok = 0
  #
  return(ok)
#
#
###############################################################################
#
# Check for a valid SAC file (returns 1 if SAC File, 0 if not a SAC File)
#
#  Right now the checks are very basic (what else can you do?)
#    I check for a positive npts, i could use the version, but
#    this isn't even listed in the ASCII definition of the header...
#    implements a file-size check where npts*4 + headerbytes (632)
#    is compared with the true size of the file.
#
#  usage:   ok = IsXYSACfile(name,hf,hi,hs)
#
#           if ok = 1, it is a SAC file, if ok = 0, it's not
#   
def IsXYSACfile(name,hf,hi,hs):
  """Test for a valid SAC file using arrays: IsSACfile(path,hf,hi,hs)
  Return value is a one if valid, zero if not."""
  #
  ok = 1
  #
  # size check info
  #
  npts = GetHvalue('npts',hf,hi,hs)
  st = os.stat(name) #file's size = st[6] 
  sizecheck = st[6] - (632 + 2 * 4 * npts)
  #
  # get the SAC file version number
  #
  version = GetHvalue('nvhdr',hf,hi,hs)
  #
  # if any of these conditions are true,
  #   the file is NOT a SAC file
  #
  if sizecheck != 0:
    ok = 0
  elif ((version < 0) or (version > 20) ):
    ok = 0
  #
  return(ok)
#
###############################################################################
#
#  usage:   [hf, hi, hs, ok] = ReadSacHeader(file_name)
#  returns:
#          The header values (70 floats, 40 integers, 
#          192 characters) are returned in an array of floats,
#          an array of integers, a buffer of bytes, 
#          and a list of strings.
#          if ok = 1, it succeeded, if ok = 0, it failed
#
def ReadSacHeader(fname):
  """\nRead a header value into the header arrays 
\t[hf,hi,hs,ok] = ReadSacHeader(thePath)
The header is split into three arrays - floats, ints, and strings
The "ok" value is one if no problems occurred, zero otherwise.\n
AM (2007): This version uses a call to platform.architecture() to test
for a 64bit system and set the size of the integer arrays to be compatible
with the sac format (which is fixed at 32bit)\n"""
  #
  hf = array.array('f') # allocate the array for header floats
  # allocate the array for header ints
  (bits,linkage)= platform.architecture()
  if bits == '64bit':
    hi = array.array('i') # this is 4 bytes on my debian 64bit system
  else :
    hi = array.array('l') # allocate the array for header ints
  hs = array.array('c') # allocate the array for header characters
  #--------------------------------------------------------------
  # open the file
  #
  try:
    f = open(fname,'r')
    #--------------------------------------------------------------
    # parse the header
    #
    # The sac header has 70 floats, then 40 integers, then 192 bytes
    #    in strings. Store them in array (an convert the char to a
    #    list). That's a total of 632 bytes.
    #--------------------------------------------------------------
    hf.fromfile(f,70)     # read in the float values
    hi.fromfile(f,40)     # read in the int values
    hs.fromfile(f,192)    # read in the char values
    #
    f.close()
    #
    ok = IsSACfile(fname,hf,hi,hs)
    #
    #--------------------------------------------------------------
  except:
    # make sure we close the file
    ok = 0
    #if f.close() == 0: # zero means the file is open
    #  f.close()
  #--------------------------------------------------------------
  #
  return([hf, hi, hs, ok])
#
###############################################################################
#
#  usage:   ok = WriteSacHeader(file_name,hf, hi, hs)
# 
#
def WriteSacHeader(fname,hf, hi, hs):
  """\nWrite a header value to the disk 
\tok = WriteSacHeader(thePath,hf,hi,hs)
The header is split into three arrays - floats, ints, and strings
The "ok" value is one if no problems occurred, zero otherwise.\n"""
  #--------------------------------------------------------------
  # open the file
  #
  ok = 1
  try:
    f = open(fname,'r+') # open file for modification
    f.seek(0,0) # set pointer to the file beginning
    # write the header
    hf.tofile(f)
    hi.tofile(f)
    hs.tofile(f)
    f.close()
    #
    #--------------------------------------------------------------
  except:
    # make sure we close the file
    #if f.closed == 0: # zero means the file is open
    #  f.close()
    ok = 0
  #
  return(ok)
  #--------------------------------------------------------------
#
#
###############################################################################
#
#  usage:   [hf, hi, hs, seis, ok] = ReadSacFile(file_name)
#  returns:
#          The header values (70 floats, 40 integers, 
#          192 characters) are returned in an array of floats,
#          an array of integers, a buffer of bytes, 
#          and a list of strings.
#          seis is an array of floats (the seismogram)
#          if ok = 1, it succeeded, if ok = 0, it failed
#
def ReadSacFile(fname):
  """\nRead read in the header and data in a SAC file 
\t[hf,hi,hs,seis,ok] = ReadSacFile(thePath)
The header is split into three arrays - floats, ints, and strings and the
data points are returned in the array seis
The "ok" value is one if no problems occurred, zero otherwise.\n
AM (2007): This version uses a call to platform.architecture() to test
for a 64bit system and set the size of the integer arrays to be compatible
with the sac format (which is fixed at 32bit)\n"""
  #
  seis = array.array('f') # allocate the array for the points
  hf = array.array('f') # allocate the array for header floats
  # allocate the array for header ints
  (bits,linkage)= platform.architecture()
  if bits == '64bit':
    hi = array.array('i') # this if 4 bytes on by debian 64bit system
  else :
    hi = array.array('l') # allocate the array for header ints
  hs = array.array('c') # allocate the array for header characters
  #--------------------------------------------------------------
  # open the file
  #
  try:
    f = open(fname,'rb')
    #--------------------------------------------------------------
    # parse the header
    #
    # The sac header has 70 floats, then 40 integers, then 192 bytes
    #    in strings. Store them in array (an convert the char to a
    #    list). That's a total of 632 bytes.
    #--------------------------------------------------------------
    hf.fromfile(f,70)     # read in the float values
    hi.fromfile(f,40)     # read in the int values
    hs.fromfile(f,192)    # read in the char values
    #
    # only continue if it is a SAC file
    #
    ok = IsSACfile(fname,hf,hi,hs)
    if ok:
      #--------------------------------------------------------------
      # read in the seismogram points
      #--------------------------------------------------------------
      npts = hi[9]  # you just have to know it's in the 10th place
      #             # actually, it's in the SAC manual
      #
      mBytes = npts * 4
      #
      seis = array.array('f')
      seis.fromfile(f,npts) # the data are now in s
      f.close()
    #
    #--------------------------------------------------------------
  except:
    # make sure we close the file
    #if f.closed == 0: # zero means the file is open
    #  f.close()
    ok = 0
  #--------------------------------------------------------------
  #
  return([hf, hi, hs, seis, ok])
#
#
###############################################################################
#
#  usage:   [hf, hi, hs, seis, ok] = ReadXYSacFile(file_name)
#  returns:
#          The header values (70 floats, 40 integers, 
#          192 characters) are returned in an array of floats,
#          an array of integers, a buffer of bytes, 
#          and a list of strings.
#          seis is an array of floats (the seismogram)
#          if ok = 1, it succeeded, if ok = 0, it failed
#
def ReadXYSacFile(fname):
  """\nRead a SAC XY file (not tested much) 
\t[hf,hi,hs,x,y,ok] = ReadXYSacFile(thePath)
The header is split into three arrays - floats, ints, and strings.
The data are in two floating point arrays x and y.
The "ok" value is one if no problems occurred, zero otherwise.\n
AM (2007): This version uses a call to platform.architecture() to test
for a 64bit system and set the size of the integer arrays to be compatible
with the sac format (which is fixed at 32bit)\n"""
  #
  x = array.array('f')
  y = array.array('f')
  hf = array.array('f') # allocate the array for header floats
  (bits,linkage)= platform.architecture()
  if bits == '64bit':
    hi = array.array('i') # this if 4 bytes on by debian 64bit system
  else :
    hi = array.array('l') # allocate the array for header ints
  hs = array.array('c') # allocate the array for header characters
  #--------------------------------------------------------------
  # open the file
  #
  try:
    f = open(fname,'r')
    #--------------------------------------------------------------
    # parse the header
    #
    # The sac header has 70 floats, then 40 integers, then 192 bytes
    #    in strings. Store them in array (an convert the char to a
    #    list). That's a total of 632 bytes.
    #--------------------------------------------------------------
    hf.fromfile(f,70)     # read in the float values
    hi.fromfile(f,40)     # read in the int values
    hs.fromfile(f,192)
    #
    # only continue if it is a SAC file
    #
    #
    ok = IsXYSACfile(fname,hf,hi,hs)
    if ok:
      #--------------------------------------------------------------
      # read in the seismogram points
      #--------------------------------------------------------------
      npts = hi[9]  # you just have to know it's in the 10th place
      #             # actually, it's in the SAC manual
      #
      mBytes = npts * 4
      #
      y.fromfile(f,npts) # the data are now in s
      x.fromfile(f,npts) # the data are now in s
      #
      f.close()
    #
    #--------------------------------------------------------------
  except:
    # make sure we close the file
    if f.closed == 0: # zero means the file is open
      f.close()
    ok = 0
  #--------------------------------------------------------------
  #
  return([hf, hi, hs, x, y, ok])
#
###############################################################################
#
#   this is set up to accept the header in two arrays and a string
#       and the data as an array.
#
#   usage: ok = WriteSacFile('test',hf, hi, hs, seis)
#
#          if ok = 1, it succeeded, if ok = 0, it failed
#
#
def WriteSacBinary(ofname, hf, hi, hs, seis):
  """\nWrite a SAC file using the head arrays and array seis 
\t[ok] = WriteSacBinary(thePath,hf,hi,hs,seis)
The "ok" value is one if no problems occurred, zero otherwise.\n"""
  try:
    f = open(ofname,'wb+')
    hf.tofile(f)
    hi.tofile(f)
    hs.tofile(f)
    seis.tofile(f)
    f.close()
    ok = 1
  except:
    print 'Error writing file ', ofname
    ok = 0
  return(ok)
#
###############################################################################
#
#   usage: ok = WriteSacBinaryXY('test',hf, hi, hs, seis)
#
#          if ok = 1, it succeeded, if ok = 0, it failed
#
#
def WriteSacBinaryXY(ofname, hf, hi, hs, x,y):
  """\nWrite a SAC file using the head arrays and arrays x,y 
\t[ok] = WriteSacBinary(thePath,hf,hi,hs,x,y)
The "ok" value is one if no problems occurred, zero otherwise.\n"""
  try:
    SetHvalue('iftype',enum_dict['IXY'],hf,hi,hs)
    SetHvalue('leven',FALSE,hf,hi,hs)
    f = open(ofname,'wb+')
    hf.tofile(f)
    hi.tofile(f)
    hs.tofile(f)
    y.tofile(f)
    x.tofile(f)
    f.close()
    ok = 1
  except:
    print 'Error writing file ', ofname
    ok = 0
  return(ok)
#

###############################################################################
#
#  These are some simply utility routines to print out header
#     values if they are not equal to the 'undefined'
#     value of -12345
#  
###############################################################################
def PrintIValue(label='=', value=-12345):
  """Convenience function for printing undefined integer header values"""
  if value != -12345:
    print label, value
###############################################################################
def PrintFValue(label='=', value=-12345.0):
  """Convenience function for printing undefined float header values"""
  if value != -12345.0:
    print '%s %.8g' % (label, value)
###############################################################################
def PrintSValue(label='=', value='-12345'):
  """Convenience function for printing undefined string header values"""
  if value != '-12345':
    print label, value
#
###############################################################################
#
# a convenience function to list common header values (old version)
#  depracated - use the next subroutine that does not need the seis array
#
def ListStdValues(hf,hi,hs,seis): # h is a header list, s is a float list
  """ Convenience function for printing common header values
  ListStdValues(hf,hi,hs,seis)"""
  #
  # Seismogram Info:
  #
  nzyear = GetHvalue('nzyear',hf,hi,hs)
  nzjday = GetHvalue('nzjday',hf,hi,hs)
  [month, date] = AM_subs.jday_to_month_day(nzyear, nzjday)
  print '%s %2.2d/%2.2d/%d (%d) %d:%d:%d.%d' % ('\nReference Time = ',    \
              month, date, \
                    GetHvalue('nzyear',hf,hi,hs), \
                                GetHvalue('nzjday',hf,hi,hs), \
                                GetHvalue('nzhour',hf,hi,hs), \
                                GetHvalue('nzmin',hf,hi,hs),  \
                                GetHvalue('nzsec',hf,hi,hs),  \
                                GetHvalue('nzmsec',hf,hi,hs))
  PrintIValue('Npts  = ',GetHvalue('npts',hf,hi,hs))
  PrintFValue('Delta = ',  GetHvalue('delta',hf,hi,hs)  )
  PrintFValue('Begin = ',  GetHvalue('b',hf,hi,hs)  )
  PrintFValue('End   = ',  GetHvalue('e',hf,hi,hs)  )
  PrintFValue('Min   = ',  GetHvalue('depmin',hf,hi,hs)  )
  PrintFValue('Mean  = ',  GetHvalue('depmen',hf,hi,hs)  )
  PrintFValue('Max   = ',  GetHvalue('depmax',hf,hi,hs)  )
  #
  PrintIValue('Header Version = ',GetHvalue('nvhdr',hf,hi,hs))
  #
  # station Info:
  #
  PrintSValue('Station = ',     GetHvalue('kstnm',hf,hi,hs))
  PrintSValue('Channel = ',     GetHvalue('kcmpnm',hf,hi,hs))
  PrintFValue('Station Lat  = ',GetHvalue('stla',hf,hi,hs))
  PrintFValue('Station Lon  = ',GetHvalue('stlo',hf,hi,hs))
  PrintFValue('Station Elev = ',GetHvalue('stel',hf,hi,hs))
  #
  # Event Info:
  #
  PrintSValue('Event       = ',GetHvalue('kevnm',hf,hi,hs))
  PrintFValue('Event Lat   = ',GetHvalue('evla',hf,hi,hs))
  PrintFValue('Event Lon   = ',GetHvalue('evlo',hf,hi,hs))
  PrintFValue('Event Depth = ',GetHvalue('evdp',hf,hi,hs))
  PrintFValue('Origin Time = ',GetHvalue('o',hf,hi,hs))
  #
  PrintFValue('Azimuth        = ',GetHvalue('az',hf,hi,hs))
  PrintFValue('Back Azimuth   = ',GetHvalue('baz',hf,hi,hs))
  PrintFValue('Distance (km)  = ',GetHvalue('dist',hf,hi,hs))
  PrintFValue('Distance (deg) = ',GetHvalue('gcarc',hf,hi,hs))
#
#
###############################################################################
#
# a convenience function to list common header values (use this version)
#
def ListStdValues(hf,hi,hs): # h is a header list, s is a float list
  """ Convenience function for printing common header values
  ListStdValues(hf,hi,hs,seis)"""
  #
  # Seismogram Info:
  #
  nzyear = GetHvalue('nzyear',hf,hi,hs)
  nzjday = GetHvalue('nzjday',hf,hi,hs)
  [month, date] = AM_subs.jday_to_month_day(nzyear, nzjday)
  print '%s %2.2d/%2.2d/%d (%d) %d:%d:%d.%d' % ('\nReference Time = ',    \
              month, date, \
                    GetHvalue('nzyear',hf,hi,hs), \
                                GetHvalue('nzjday',hf,hi,hs), \
                                GetHvalue('nzhour',hf,hi,hs), \
                                GetHvalue('nzmin',hf,hi,hs),  \
                                GetHvalue('nzsec',hf,hi,hs),  \
                                GetHvalue('nzmsec',hf,hi,hs))
  PrintIValue('Npts  = ',GetHvalue('npts',hf,hi,hs))
  PrintFValue('Delta = ',  GetHvalue('delta',hf,hi,hs)  )
  PrintFValue('Begin = ',  GetHvalue('b',hf,hi,hs)  )
  PrintFValue('End   = ',  GetHvalue('e',hf,hi,hs)  )
  PrintFValue('Min   = ',  GetHvalue('depmin',hf,hi,hs)  )
  PrintFValue('Mean  = ',  GetHvalue('depmen',hf,hi,hs)  )
  PrintFValue('Max   = ',  GetHvalue('depmax',hf,hi,hs)  )
  #
  PrintIValue('Header Version = ',GetHvalue('nvhdr',hf,hi,hs))
  #
  # station Info:
  #
  PrintSValue('Station = ',     GetHvalue('kstnm',hf,hi,hs))
  PrintSValue('Channel = ',     GetHvalue('kcmpnm',hf,hi,hs))
  PrintFValue('Station Lat  = ',GetHvalue('stla',hf,hi,hs))
  PrintFValue('Station Lon  = ',GetHvalue('stlo',hf,hi,hs))
  PrintFValue('Station Elev = ',GetHvalue('stel',hf,hi,hs))
  #
  # Event Info:
  #
  PrintSValue('Event       = ',GetHvalue('kevnm',hf,hi,hs))
  PrintFValue('Event Lat   = ',GetHvalue('evla',hf,hi,hs))
  PrintFValue('Event Lon   = ',GetHvalue('evlo',hf,hi,hs))
  PrintFValue('Event Depth = ',GetHvalue('evdp',hf,hi,hs))
  PrintFValue('Origin Time = ',GetHvalue('o',hf,hi,hs))
  #
  PrintFValue('Azimuth        = ',GetHvalue('az',hf,hi,hs))
  PrintFValue('Back Azimuth   = ',GetHvalue('baz',hf,hi,hs))
  PrintFValue('Distance (km)  = ',GetHvalue('dist',hf,hi,hs))
  PrintFValue('Distance (deg) = ',GetHvalue('gcarc',hf,hi,hs))
#
###############################################################################
#
def GetHvalueFromFile(thePath,theItem):
  """\nQuick access to a specific header item in a specified file.
  GetHvalueFromFile(thePath,theItem)
returns -12345 if a problem occurred.\n"""
  #
  #  Read in the Header
  #
  [hf, hi, hs, ok] = ReadSacHeader(thePath)
  #
  if ok:
    return(GetHvalue(theItem,hf, hi, hs))
  else:
    print "Problem in GetHvalueFromFile."
    return(-12345)
#
###############################################################################
#
def SetHvalueInFile(thePath,theItem,theValue):
  """\nQuick access to change a specific header item in a specified file.
  SetHvalueFromFile(thePath,theItem)
The "ok" value is one if no problems occurred, zero otherwise.\n"""
  #
  ok = 0
  #
  #  Read in the Header
  #
  [hf, hi, hs, r_ok] = ReadSacHeader(thePath)
  #
  if r_ok:
    before = GetHvalue(theItem,hf,hi,hs)
    if before != 'NULL':
      c_ok = SetHvalue(theItem,theValue,hf, hi, hs)
      if c_ok:
        after = GetHvalue(theItem,hf,hi,hs)
        ok = WriteSacHeader(thePath,hf,hi,hs)  
        #
        #print 'Changed ',before,' to ',after,' in ',thePath
  else:
    print "Problem in SetHvalueInFile."
  #
  return(ok)


#
###############################################################################
#
def IsValidSacFile(thePath):
  """\nQuick test for a valid SAC binary file file.
  IsValidSACFile(thePath)
The "ok" value is one if no problems occurred, zero otherwise.\n"""
  #
  #  Read in the Header
  #
  [hf, hi, hs, ok] = ReadSacHeader(thePath)
  #
  if ok:
    ok = IsSACfile(thePath,hf,hi,hs)
  else:
    ok = 0
  #
  return(ok)
  #

#
###############################################################################
#
if __name__ == "__main__":
  pass

