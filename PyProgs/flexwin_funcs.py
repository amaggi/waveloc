"""
Functions for implementing flexwin-style selection of max-corr maxima.
"""

from waveloc_funcs import *
from numpy import abs, exp

def setup_MLR(values,i_maxima,i_minima,w_level_max,width,w_level_prominence):
  """
  Finds all possible windows within trace, with central maximum
  above the water level.  Returns a list of tuples (iM,iL,iR) with
  indexes of central maximum, left-most point and right-most point
  of each window.
  """
  iMLR=[]
  
  imin_first=i_minima[0]
  imin_last=i_minima[-1]
  imax_first=i_maxima[0]
  imax_last=i_maxima[-1]
  
  # create a reversed list of minima and maxima indexes
  i_minima_reversed=i_minima[:]
  i_minima_reversed.reverse()
  i_maxima_reversed=i_maxima[:]
  i_maxima_reversed.reverse()
  
  for imax in i_maxima:
    #only continue if there are available minima on either side
    if imax > imin_first and imax < imin_last :
      # only continue if this maximum is above the water level
      if values[imax] > w_level_max :
           
        imax_index=i_maxima.index(imax)
        imax_index_rev=i_maxima_reversed.index(imax)
        
        # find the first maximum on the right that either goes below the water level
        # or is higher than me
        R_max=imax_last
        for imax2 in i_maxima[imax_index:]:
          if values[imax2] > values[imax]  or imax2==imax_last:
            R_max = imax2
            break
        # find the first maximum on the left that either goes below the water level
        # or is higher than me 
        L_max=imax_first
        for imax2 in i_maxima_reversed[imax_index_rev:]:
          if values[imax2] > values[imax] or imax2==imax_first:
            L_max = imax2
            break

        # find the first minimum right of imax
        R=imin_last
        for imin in i_minima:
          if imin > imax :
            R = imin
            break
        # find the first minimum left of imax
        L=imin_first
        for imin in i_minima_reversed:
          if imin < imax :
            L = imin
            break
            
        L_index=i_minima.index(L)
        R_index=i_minima.index(R)
        
#        # apply prominence criteria : if fail, go to next imax
#        if values[imax]-values[L] < w_level_prominence or values[imax]-values[R] < w_level_prominence: 
#          continue
      
        # find the first minimum left of R_max
        RR=imin_last
        for imin in i_minima_reversed:
          if imin < R_max :
            RR = imin
            break
        # find the first minimum right of L_max
        LL=imin_first
        for imin in i_minima:
          if imin > L_max :
            LL = imin
            break
          
        # iterate over minima, going left and right to get all
        # possible windows around imax         
        
        LL_index=i_minima.index(LL)
        RR_index=i_minima.index(RR)
        
                
        for i_left in i_minima[LL_index:L_index] :
          for i_right in i_minima[R_index:RR_index] :
            if (i_right - i_left) > width:
              iMLR.append((imax,i_left,i_right))

  return iMLR
   
   
def reject_on_water_level(values,iMLR,i_minima,w_level):
  """
  Removes windows with internal minima below the water level
  """
  # for each proto-window, check that no internal minima
  # fall below the water level
  
  iMLR_copy=iMLR[:]
  
  for (iM,iL,iR) in iMLR :
    index_iL=i_minima.index(iL)
    index_iR=i_minima.index(iR)
    for imin in i_minima[index_iL:index_iR] :
      if values[imin] < w_level:
        iMLR_copy.remove((iM,iL,iR))        
        break

  return iMLR_copy      	
  
  
def reject_on_window_width(values,iMLR,width):
  """
  Removes windows that are too short
  """
  
  iMLR_copy=iMLR[:]
  
  for (iM,iL,iR) in iMLR :   
    if (iR - iL) < width:
      iMLR_copy.remove((iM,iL,iR))  
                	
  return iMLR_copy   

 
def reject_on_center_maximum(values,iMLR,i_maxima):
  """
  Removes windows for which the center maximum is not the 
  highest point
  """
  iMLR_copy=iMLR[:]

  
  for (iM,iL,iR) in iMLR :
    center_height=values[iM]

    for imax in i_maxima :
      if imax > iL and imax < iR and values[imax] > center_height:
        iMLR_copy.remove((iM,iL,iR))
        break
      
  return iMLR_copy	
  
  
def reject_on_separation(values,iMLR,i_maxima,i_minima,C3a,separation):
  """
  Removes windows for which there are multiple distinct peaks.
  """
  iMLR_copy=iMLR[:]

  
  for (iM,iL,iR) in iMLR :
    # find the lowest minimum within the window
    min_value=values[iL]
    for imin in i_minima[iL:iR] :
      if values[imin] < min_value:
        min_value=values[imin]
    # find height of central maximum above this minimum
    d_iM=values[iM]-min_value
#    print "Window %d %d %d : max = %.2f, d_max = %.2f"%(iL,iM,iR,values[iM],d_iM)
        
    # now apply separation condition on all other maxima 
    # in the window     
    for imax in i_maxima :
      if imax > iL and imax < iR and not (imax == iM) :
        d_imax = values[imax]-min_value
        d_time = abs(iM-imax)
        # evaluate time decay function
        if d_time >= separation : 
          f_time = exp(-((d_time-separation)/separation)**2)
        else:
          f_time = 1.0
        # check condition
        if d_imax > C3a*d_iM*f_time :
#          print "max=%.2f, imax=%d, d_imax=%.2f, d_time=%.2f, f_time=%.2f : failed"%(values[imax],imax,d_imax,d_time,f_time)
          iMLR_copy.remove((iM,iL,iR))
          break
          
  return iMLR_copy 
     
	
 
def remove_subwindows(values,iMLR,separation):
  """
  If a window is entirely within another, and its central
  maximum is below the central maximum of another window,
  then remove it.
  """  

  iMLR_copy=iMLR[:]
  
  for win in iMLR :
    (iM,iL,iR)=win
    center_height=values[iM]
    # find windows that contain me
    containing_wins=[]
    for cwin in iMLR:
      (iM_cwin, iL_cwin, iR_cwin) = cwin
      if not (cwin==win)  and iL >= iL_cwin and iR<= iR_cwin and abs(iM - iM_cwin) < separation:
        # cwin contains me, so check that its maximum is
        # higher than or equal to  mine
          if values[iM_cwin] > values[iM]:
            # I am no longer useful
            iMLR_copy.remove((iM,iL,iR))
            break
      
  return iMLR_copy  
  

def unique_maxima(values,iMLR,separation):
  
    max_indexes=[iM for (iM,iL,iR) in iMLR]
    max_indexes.sort()
    
    # make this list unique
    indexes_copy=max_indexes[:]
    for imax in indexes_copy:
        num_i=max_indexes.count(imax)
        if num_i > 1:
           max_indexes.remove(imax) 
    
    indexes_copy=max_indexes[:]
    for imax in indexes_copy:
      for imax2 in indexes_copy:
        if not imax == imax2 and abs(imax2-imax)<separation and values[imax2] >= values[imax]:
          max_indexes.remove(imax)
          break
      
           
    return max_indexes
    
def print_iMLR(A,iMLR,b,dt):
  for (iM,iL,iR) in iMLR:
    print "%.2f %.2f %.2f : %.2f"%(b+iL*dt,b+iM*dt,b+iR*dt,A[iM])
  print "---"
  
