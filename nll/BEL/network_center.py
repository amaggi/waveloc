#!/usr/bin/env python
# encoding: utf-8

import os, sys

import numpy as np

station_file='bel_stations.txt'

file=open(station_file,'r')
lines=file.readlines()
file.close()

deglats=[np.float(line.split()[3]) for line in lines ]
deglons=[np.float(line.split()[4]) for line in lines ]
names=[line.split()[1] for line in lines ]

cenlat=np.mean(deglats)
cenlon=np.mean(deglons)

print deglats, deglons, names

print cenlat, cenlon
