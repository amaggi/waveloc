#!/usr/bin/env python

import os, glob
from hdf5_grids import nll2hdf5

time_glob=glob.glob('*.time.buf')
for bfile in time_glob:
  nll_name,ext= os.path.splitext(bfile)
  h5_name="%s.hdf5"%nll_name
  print h5_name
  nll2hdf5(nll_name,h5_name)
