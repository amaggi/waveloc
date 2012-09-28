.. WaveLoc Classes and other components

###########
WaveLoc API
###########

Waveform Manipulation
=====================

The main module for waveform manipulation is :mod:`OP_waveforms`.  This module
is a second-generation waveform processing module, adapted to use the
functionality provided by ``obspy.stream`` objects (see the `obspy website
<http://www.obspy.org>`_ for more information and installation instructions).


OP_waveforms
------------

.. automodule:: OP_waveforms

.. autoclass:: Waveform
  :members: bp_filter, cleanup_traces, decimate, process_kurtosis,
            read_from_SDS, read_from_file, resample, rmean, taper, write_to_file,
            write_to_file_filled


obspyaux
--------

.. automodule:: obspyaux
  :members:


filters
-------

.. automodule:: filters
  :members:


Grid Manipulation
=================

AM_geo
------
.. automodule:: AM_geo
  :members:


grids_paths
-----------
.. automodule:: grids_paths
  :members: QDGrid, QDTimeGrid


Non-Lin-Loc IO
==============

.. automodule:: NLL_IO
  :members:


Other
=====

AM_subs
-------

.. automodule:: AM_subs
  :members:

waveloc_funcs
-------------
.. automodule:: waveloc_funcs
  :members:
 
