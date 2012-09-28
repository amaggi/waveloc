.. Introduction to Waveloc stuff

============
Introduction
============

WaveLoc is a simple detector and preliminary locator for seismic phenomena.  It
is based on the central principle that the coherence of arrivals across a network
contains the necessary information for both detection and preliminary location.
It's main object is to provide an alternative to traditional detection /
location by picking and associating discrete arrival times, a process that
becomes complex during earthquake swarms or aftershock sequences.

As the basic source of information for single-event location are the
(physically coherent) arrival times at seismic stations, we enhance the first
arrival information on the the recorded seismograms by applying a kurtosis
filter.

We then migrate the kurtosis waveforms back to a grid of possible event
locations, and stack.  The maxima of the stack act as both detection and
preliminary location of the events.  

Here is the result of an example run:

.. image:: figures/loc_example.png
  :width: 600px
  :align: center


The method is explained in greater detail below, and in a manuscript currently being written.

The kurtosis filter
===================
The kurtosis is the fourth moment of the data, and indicates how non-gaussian a
distribution is.  

.. math::
  {\rm Kurt}(x_1,\ldots,x_N) = \left\{ \frac{1}{N}\sum_{j=1}^N \left[ \frac{x_j - \bar{x}}{\sigma}\right]^4 \right\}-3, \label{eq:kurt}

Distributions that are more peaked than gaussian
distributions are will have a strong positive value of the kurtosis.
Distributions that are flatter than gaussian distributions will have negative
values.  The arrivals of seismic phases are strongly non-gaussian, and tend to
produce strongly peaked amplitude distributions, and hence positive kurtosis
values. The -3 in the formula above normalises the kurtosis to zero for a
gaussian distribution.  See kurtosis column in the introductory figure for an
example.

As the kurtosis calculation in the above equation can be quite cumbersome to
calculate on large volumes of data, we use instead the following formulation
for a recursive kurtosis calculation:

.. math::
  \bar{x}_i = C\ \bar{x}_{i-1} + (1-C)\ x_i, \\
  \sigma^2_i = C\ \sigma_{i-1} + (1-C)\ (x_i-\bar{x}_i)^2, \\
  {\rm Kurt}_i = C\ {\rm Kurt}_{i-1} + (1-C)\ \left[\frac{(x_i-\bar{x}_i)}{\sigma_i}\right]^4 ,

where the constant C is the ratio of the time-step of the data and a chosen
window length for the kurtosis calculation.
  

Migration
=========


Detection / Location
====================
