def time2nztimesac(stringa):
    """
    this function converts the onbspy strings stats.starttime and stats.endtime 
    (eg, '2009-04-06T01:29:59.840000Z') into
    nzyear, month, day, hour, minute, seconds and milliseconds
    It returns the latter values as integers in a tuple
    """

    starttime_tuple=str(stringa).split("T")
    #print starttime_tuple
    (nzyear,month,day)=starttime_tuple[0].split("-")
    #print int(nzyear),int(month),int(day)
    #
    starttime_homisec=starttime_tuple[1]
    # remove the last character
    if starttime_homisec[-1:]=='Z':
       starttime_homisec=starttime_homisec[:-1]
    #determine the milliseconds
    tmp=starttime_homisec.split(".")
    msec=int(tmp[1])
    msec= int(msec/1.0e3)
    # determine the hour, minute and seconds
    (ho,mi,sec)=tmp[0].split(":")
    #print ho,mi,sec
    return int(nzyear),int(month),int(day),int(ho),int(mi),int(sec),msec
