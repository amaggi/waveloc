.. WaveLoc Classes and other components

###########
WaveLoc API
###########

NllGridLib
----------

Waveloc needs travel-times to migrate the kurtosis waveforms. These can in
principle be obtained in many different ways. For simple Earth models, a 
convenient way to obtain travel-time grids is to take advantage of the 
`NonLinLoc <http://alomax.free.fr/nlloc/>`_ implementation of the eikonal
solver of Podvin & Lecomte (1991).

.. automodule:: NllGridLib
    :members:


Waveforms (OP_waveform)
=======================

The main module for waveform manipulation is :mod:`OP_waveforms`.  This module
is a second-generation waveform processing module, adapted to use the
functionality provided by ``obspy.stream`` objects (see the `obspy website
<http://www.obspy.org>`_ for more information and installation instructions).

.. automodule:: OP_waveforms

.. autoclass:: Waveform



