.. Tutorial for WaveLoc

========
Tutorial
========

Waveloc options
===============

Preprocess your data
====================

Raw data can be provided in any format that obspy is capable of reading (tested
formats : ``MSEED`` and ``SAC`` for now).  Preferred data structure are SDS archives,
but other structures are acceptable so long as a siutable io interface is coded
into the processing code.  As an alternative to re-programming the processing
code, a helper script, ``make_SDS_data_links.py``, is provided to create SDS
archives from flat files using soft links.

The processing code will look for data in the ``$WAVELOC_PATH/data/SOME_NAME``
directory, where the ``SOME_NAME`` subdirectory contains the SDS archive(s) and will
also contain the processed waveforms.  

Most of the options are self-explanatory.  The choice of the frequency ranges
will influence the ability of the kurtosis to detect a non-stationary signal
(i.e. a first arrival).  We provide a helper-code ``run_kurtogram.py`` that may
help you find the best frequency range for your data by plotting the power in
the kurtosis as a function of central frequency and of width of frequency band.
Run this on several examples of data to get an idea of the best frequency range
for the rest of the analysis.

We suggest you use the recursive
kurtosis calculation given the considerable increase in speed it offers.  The
window for the recursive calculation should be long enough to smooth over the
S-P arrival time, in order to avoid double peaks in the kurtosis.
We also suggest you calculate the gradient of the kurtosis, as using this for 
migration will give more precise results. 

.. warning::
  The option to chose recursive kurtosis and kurtosis gradient calculations may
  soon disappear from the code, as they will likely become the default options.

We provide an example run-script ``run_tutorial.sh`` to show how to run the
processing scripts and the following migration and location scripts.  Feel free
to modify this run-script for your own data.

::

  Add file here

Show figure for output.


Prepare the time and search grids
=================================

WaveLoc works on the basis of pre-determined grids.  For programming
convenience, these grids are regularly spaced, and based on the grid format
used by NonLinLoc.  Details of the format can be found by reading the NLL
documentation at http://alomax.free.fr/nlloc or the relevant routines of the
WaveLoc code :mod:`grids_paths`. 

Setting up one of these grids using the NLL utilities requires first of all
setting up an input file for NLL.  The example file ``nll_tutorial.in`` in the 
``$WAVELOC_PATH/nll/TUTORIAL`` directory contains
all the relevant statements, with lots of comments to explain what they are
for, but here is the short story for a 3-D representation of a 1-D velocity model: 

========= =============================================
Statement Comment
========= =============================================
TRANS     Choice of geographical transformation to use.
          All subsequent grids will be locked to this
          transformation, and counted from its center.
VGOUT     Root name for the output files containing the 
          velocity grid.
VGTYPE    Wave-type for the velocity grids (P or S).
VGGRID    Definition of the extent and grid-step for the
          velocity grid in 3D (relative to the geographical
          projection
LAYER     Multiple statements to describe the 1D layered velocity model
GTFILES   Root name for the output files containing the time grids (one per station).
GTMODE    Defines the mode of the time grids (2D/3D, with or without angles).
GTSRCE    Multiple lines containing the station coordinates.
========= =============================================

Make sure that your velocity/time grids are large enough to contain the stations. If a station is within the grid, you should get this sort of output::

  > Vel2Grid nll_tutorial.in 

  Calculating travel times for source: BLS1  X 33.5463  Y -31.7444  Z -0.1000 ...
  Source:  Velocity: 3.400000 km/sec  GridLoc: ix=103.546295 iy=38.255554 iz=2.900000

  Recursive initialization: level 1
  Homogeneous region: x[13->16] y[15->18] z[0->7]

  Starting F.D. computation...
  Starting F.D. computation...
  time_3d: Computations terminated normally.
  Finished calculation, time grid output files: ./tutorial.P.BLS1.*
  
If a station is outside the grid, you will get this sort of error message::

  > Vel2Grid nll_tutorial.in 
  Calculating travel times for source: OBO2  X 117.5147  Y 12.2111  Z -0.0450 ...
  Grid2Time: ERROR: Source point is not inside model grid.
  Source:  GridLoc: ix=187.514709 iy=82.211113 iz=2.955000
  Grid2Time: ERROR: calculating travel times.
  
Having stations outside of the time grid will not make WaveLoc crash, but you
will not be able to use data from those stations (WaveLoc will tell you and
keep going for the other stations).

Create a separate text file containing only the GTSRCE lines, which will be
read by other scripts that need the station coordinates.  Copy this file to
the ``$WAVELOC_PATH/lib`` directory.

Creating the grid files requires running two NonLinLoc commands (NLL must be installed first, of course)::

  > Vel2Grid nll_tutorial.in
  > Grid2Time nll_tutorial.in

You will also need a search grid, centered around the region (within your
velocity/time grid) in which you expect to find events.  This grid can be finer
than the previous one, in which case travel-times for intermediate points will
be calculated by linear interpolation of the time grids.  However, the
geographical transformation and its center should be the same as for the time
grids. Hint: start with an identical grid, by copying the ``.mod.hdr`` file,
then modify this grid, making it smaller or denser as you prefer.  

..warning:: 
  Computation time for WaveLoc increasces linearly with the total number of
  points in the search grid, so make the grid as small as you can, and no denser
  than you need. 

Run ``grid2hdf5`` to turn the nll .hdr and .buf files into .hdf5 files. 
Copy, move or make symbolic links to the .hdf5 files in the 
``$WAVELOC_PATH/lib`` directory, and you should be set to go.





Run migration
=============

The migration code is run using the script
:mod:`run_waveloc_threading.py`, which uses a multi-threaded architecture
to process multiple data streams simultaneously. 

To see the command-line options for the processing code type::

  > run_waveloc_threading.py --help

::

  Output of help run

The migration code performs three steps : 

# Calculates travel-times for the search grid by interpolating the time-grids for each station (this operation takes some time).
# 

If your data directory contains data from stations which were outside the time grid, you will get

Run migration code.  Give example.

Run location
============

Run location code.

Optionnal run plots
===================

Make pretty pictures
