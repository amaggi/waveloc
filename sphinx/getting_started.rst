.. Tutorial for WaveLoc

===============
Getting Started
===============

This is a quick guide for getting started with WaveLoc.

Download and installation
=========================

You can download the latest distribution here: `waveloc-0.1.0.tar.gz <https://github.com/downloads/amaggi/waveloc/waveloc-0.1.0.tar.gz>`_, or if you're feeling courageous you can download the deveopment version from github: https://github.com/amaggi/waveloc .  

System requirements : 

* python, numpy, h5py (all available through the `Enthought python distribution <http://www.enthought.com/products/getepd.php>`_)
* obspy (available here : http://obspy.org)
* NonLinLoc for time grid calculation (not strictly required)
* other stuff I've surelly fogotten about...

Installation
------------

Untar the distribution, then install in the usual Python manner : ::

  python setup.py install

Set the environment variable ``WAVELOC_PATH`` to point to a directory you want to work in (make sure you have plenty of space on the corresponding disk).

If you want to use ``NonLinLoc`` to calculate the 3D time-grids, download it from Anthony Lomax's website http://alomax.free.fr/nlloc/, and install it where you prefer (you just need to be able to call ``Vel2Grid`` and ``Grid2Time`` correctly).


Running the examples
====================

In order to get you started running waveloc, we have prepared the following tests : 

#. a synthetic test to verify the response of the recording network;
#. a real migration test.

Running the synthetic test
--------------------------

Download the example scripts here : `example_scripts.tgz <https://github.com/downloads/amaggi/waveloc/example_scripts.tgz>`_, and the test data here : `test_data.tgz <https://github.com/downloads/amaggi/waveloc/test_data.tgz>`_ (beware : it is a large file !).  Unpack the test data archive in the  ``§WAVELOC_PATH`` directory, and unpack the example scripts.  Run the ``setup\_tests.py`` script to set up the required directory structure for the examples.  Run the  ``run\_syn\_test.py`` script to run the synthetic test.  The first time you run the script it will take a long time, as the time grids need to be recalculated. After the run, you should the find the following figure in the directory ``§WAVELOC_PATH/out/TEST_Dirac/fig``:
  
.. image:: figures/syn_example.png
  :width: 600px
  :align: center

Running the migration test
--------------------------

