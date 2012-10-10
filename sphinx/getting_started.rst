.. Tutorial for WaveLoc

===============
Getting Started
===============

This is a quick guide for getting started with WaveLoc.

Download and installation
=========================

You can download the latest distribution on the `waveloc downloads page  
<http://github.com/amaggi/waveloc/downloads>`_, or if you're feeling
courageous you can download the deveopment version from github:
http://github.com/amaggi/waveloc.  

System requirements : 

* python, numpy, h5py (all available through the `Enthought python distribution
  <http://www.enthought.com/products/getepd.php>`_)
* obspy (available here : http://obspy.org)
* NonLinLoc for time grid calculation (not strictly required)
* other stuff I've surelly fogotten about...

Installation
------------

Untar the distribution, then install in the usual Python manner : ::

  python setup.py install


If you want to use ``NonLinLoc`` to calculate the 3D time-grids, download it
from Anthony Lomax's website http://alomax.free.fr/nlloc/, and install it where
you prefer (you just need to be able to call ``Vel2Grid`` and ``Grid2Time``
correctly).


Running the examples
====================

In order to get you started running waveloc, we have prepared the following
example scripts : 

#. a synthetic test to verify the response of the recording network;
#. a real migration test.

You should find the relevant scripts in the ``examples`` directory in the
waveloc distribution.

* For the tests, temporarily set the environment variable ``$WAVELOC_PATH`` to
  the directory you would like the tests to run in. 

* Download the test data here : `test_data.tgz
  <https://github.com/downloads/amaggi/waveloc/test_data.tgz>`_ (beware : it is a
  large file !), and unpack the archive in the  ``$WAVELOC_PATH`` directory.

* Run the ``setup_tests.py`` script to set up the required directory structure
  for the examples.  

Running the synthetic test
--------------------------
Run the ``run_syn_test.py`` script to run the synthetic test.  The first time
you run the script it will take a long time, as the time grids need to be
interpolated.  After the run, you should the find the following figure in the
directory ``§WAVELOC_PATH/out/TEST_Dirac/fig``:
  
.. image:: figures/syn_example.png
  :width: 600px
  :align: center

Running the migration test
--------------------------

