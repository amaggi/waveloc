.. WaveLoc Classes and other components

###########
WaveLoc API
###########

This page contains the documentation on each of the modules / classes /
functions that make up the waveloc package. The modules are ordered by subject
type. This is documentation is automatically generated from the doc-strings in
the code itself. It is meant as an aide to developpers or advanced users who
need access to the intrigacies of the code. It is not a tutorial (one will be
coming soon).

TRAVEL-TIMES
------------

Waveloc needs travel-times to migrate the kurtosis waveforms. These can in
principle be obtained in many different ways. For simple Earth models, a 
convenient way to obtain travel-time grids is to take advantage of the 
`NonLinLoc <http://alomax.free.fr/nlloc/>`_ implementation of the eikonal
solver of Podvin & Lecomte (1991).

NllGridLib
==========

.. automodule:: NllGridLib
    :members:

WAVEFORMS
---------

The waveform manipulation routines in waveloc are heavily based on
`obspy <http://www.obspy.org>`_. As waveloc development started before obspy
was fully functionnal, some external functions are used where obspy equivalents
were not available. Updating all the waveform manipulation routines to use the
latest obspy features will be completed at some time.

OP_waveform
===========

.. automodule:: OP_waveforms
    :members:

SDS_processing
==============

.. automodule:: SDS_processing
    :members:

