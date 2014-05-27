#!/usr/bin/env python
# encoding: utf-8
"""
The module SDS_processing contains functions to apply Waveloc specific
processing on data contained in an SDS file structure.

"""

import os
from obspy.core import utcdatetime
from OP_waveforms import Waveform
import logging


def read_channel_file(fname):
    """
    Reads a channel file containing network-code, station-name, component-name
    triplets. Can be used to provide fine control of data to be processed.

    :param fname: File to be read
    :type fname: string

    :returns: List of tuples. Each tuple contains network-name, station-name,
        channel-name as strings.

    """
    f = open(fname, 'r')
    lines = f.readlines()
    f.close()
    triplet_list = []
    for line in lines:
        net, sta, comp = line.split()
        triplet_list.append((net, sta, comp))
    return triplet_list


def do_SDS_processing_setup_and_run(opdict):
    """
    Does all the processing for an SDS archive. All options are given within a
    *WavlocOptions* object. Steps are:

    * reading data
    * filtering
    * resampling if requested
    * applying kurtosis processing
    * calculating kurtosis gradient if requested
    * convolving with gaussian if requested
    * writing processed files at each stage

    :param opdict: Dictionary of options. Refers the the opdict attribute of
        *WavelocOptions* objects.
    :type opdict: dictionary

    """

    base_path = opdict['base_path']
    data_dir = os.path.join(base_path, 'data', opdict['datadir'])

    filter_c1 = opdict['c1']
    filter_c2 = opdict['c2']
    kurt_window = opdict['kwin']

    dataglob = opdict['dataglob']
    kurtglob = opdict['kurtglob']
    gradglob = opdict['gradglob']
    if opdict['gauss']:
        gaussglob = opdict['gaussglob']

    # start and end time to process
    start_time = utcdatetime.UTCDateTime(opdict['starttime'])
    end_time = utcdatetime.UTCDateTime(opdict['endtime'])

    # if have a channel file then read it
    if 'channel_file' in opdict:
        fname = os.path.join(base_path, 'lib', opdict['channel_file'])
        triplet_list = read_channel_file(fname)
    else:
        # else make triplet list from net, sta, comp lists
        triplet_list = []
        net_list = opdict['net_list'].split(',')
        sta_list = opdict['sta_list'].split(',')
        comp_list = opdict['comp_list'].split(',')
        for net in net_list:
            for sta in sta_list:
                for comp in comp_list:
                    triplet_list.append((net, sta, comp))

    # loop over data
    for net, sta, comp in triplet_list:
        full_path = os.path.join(data_dir, net, sta, "%s.D" % comp)
        logging.debug("Full path : %s" % full_path)
        if os.path.exists(full_path):

            # construct the base filename for output and create the wf object
            filt_filename = os.path.join(data_dir,
                                         "%s.%s.%s.%s.%s" %
                                         (start_time.isoformat(),
                                          net, sta, comp, dataglob[1:]))
            logging.debug("Processing to create %s" % (filt_filename))
            wf = Waveform()

            # read and process the data
            try:
                # read and filter data
                wf.read_from_SDS(data_dir, net, sta, comp,
                                 starttime=start_time, endtime=end_time)
                wf.bp_filter(filter_c1, filter_c2, rmean=True, taper=True)
                if opdict['resample']:
                    wf.resample(opdict['fs'])
                wf.write_to_file_filled(filt_filename, format='MSEED',
                                        fill_value=0)

                # do kurtosis processing
                kurt_filename = os.path.join(data_dir,
                                             "%s.%s.%s.%s.%s" %
                                             (start_time.isoformat(),
                                              net, sta, comp, kurtglob[1:]))
                logging.debug("Processing to create %s" % (kurt_filename))
                wf.process_kurtosis(kurt_window, recursive=opdict['krec'],
                                    pre_taper=True, post_taper=True)
                wf.write_to_file_filled(kurt_filename, format='MSEED',
                                        fill_value=0)

                # calculate kurtosis gradient if requested
                if opdict['kderiv']:
                    kurt_grad_filename = os.path.join(data_dir,
                                                      "%s.%s.%s.%s.%s" %
                                                      (start_time.isoformat(),
                                                       net, sta, comp,
                                                       gradglob[1:]))
                    logging.debug("Processing to create %s" %
                                  (kurt_grad_filename))
                    wf.take_positive_derivative(pre_taper=True,
                                                post_taper=True)
                    wf.write_to_file_filled(kurt_grad_filename, format='MSEED',
                                            fill_value=0)

                # do convolution with gaussian if requested
                if opdict['gauss']:
                    thres = opdict['gthreshold']
                    mu = opdict['mu']
                    sigma = opdict['sigma']
                    gauss_filename = os.path.join(data_dir,
                                                  "%s.%s.%s.%s.%s" %
                                                  (start_time.isoformat(),
                                                   net, sta, comp,
                                                   gaussglob[1:]))
                    logging.debug("Processing to create %s" % (gauss_filename))
                    wf.process_gaussian(thres, mu, sigma)
                    wf.write_to_file_filled(gauss_filename, format='MSEED',
                                            fill_value=0)

            except UserWarning:
                logging.info('No data within time limits for %s %s %s' %
                             (net, sta, comp))
