#!/usr/bin/env python
# encoding: utf-8

import os, sys, optparse, glob
from obspy.core import read


p=optparse.OptionParser()
p.add_option('--mseed_file', '-m', action='store',help="mseed filename")
p.add_option('--sac_file', '-s', action='store',help="sac filename")

(options,arguments)=p.parse_args()

st=read(options.mseed_file)
print(st)
st.write(options.sac_file,format='SAC')
