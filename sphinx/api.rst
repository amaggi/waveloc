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

Waveloc Options
---------------

Waveloc is piloted entirely through a set of options / parameters that are set
in a WavelocOptions object. All options will (soon!) be explained on the
tutorial pages.

options
=======

.. automodule:: options
    :members:

Travel-times
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

hdf5_grids
==========
Waveloc manipulates many large arrays, some of which are written / read from
disk. On some machines, these arrays are larger than the available RAM. In
order to avoid the drammatic slow-downs when swapping to disc occurs in these
cases, we use hdf (hierarchical data formats) that provide memory-like access
to files on disc. 

All large arrays that are infrequently accessed (e.g. travel-time arrays) are
written and used in this manner. When there is sufficient RAM available, large
arrays that are frequently accessed (search-point travel-times, 4D-migration
grids) are kept in memory.

.. automodule:: hdf5_grids
    :members:

Waveforms
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

filters
=======

.. automodule:: filters
    :members:

kurtogram
=========

.. automodule:: kurtogram
    :members:

SDS_processing
==============

.. automodule:: make_SDS_data_links
    :members:

.. automodule:: SDS_processing
    :members:

Migration
---------

migration
=========

.. automodule:: migration
    :members:

synth_migration
===============

.. automodule:: synth_migration
    :members:

integrate4D
===========

.. automodule:: integrate4D
    :members:

Location
--------

locations_trigger
=================

.. automodule:: locations_trigger
    :members:

locations_prob
==============

.. automodule:: locations_prob
    :members:

magnitude
=========

.. automodule:: magnitude
    :members:

plot_locations2
==================
.. automodule:: plot_locations2
    :members:

plot_mpl
========
.. automodule:: plot_mpl
    :members:

Clustering
----------
This part of the waveloc code was contributed by Nadège Langet during her PhD
thesis.

correlation
===========

.. automodule:: correlation
    :members:

clustering
===========

.. automodule:: clustering
    :members:

double_diff
===========

.. automodule:: double_diff
    :members:

CZ_color
===========

.. automodule:: CZ_color
    :members:
