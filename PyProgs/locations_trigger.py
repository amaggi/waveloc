#!/usr/bin/env python
# encoding: utf-8

import os
import glob
import h5py
import logging
import matplotlib.pyplot as plt
import numpy as np
from obspy.core import utcdatetime, read
from obspy.signal import trigger
from OP_waveforms import Waveform
from filters import smooth
from hdf5_grids import get_interpolated_time_ugrids
from ugrids import ugrid_closest_point_index


def plot_location_triggers(trace, trig_start, trig_end, trig_95_start,
                           trig_95_end, show=True):
    """
    Plot the location (in time) of the triggers.
    TODO : Flesh out this doc-string

    :param trace:
    :param trig_start:
    :param trig_end:
    :param trig_95_start:
    :param trig_95_end:
    :param show:
    """
    df = trace.stats.sampling_rate
    npts = trace.stats.npts
    t = np.arange(npts) / df
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax1.plot(t[trig_start-100:trig_end+100],
             trace.data[trig_start-100:trig_end+100], 'k')
    i, j = ax1.get_ylim()
    try:
        ax1.vlines(trig_start / df, i, j, color='r', lw=1, label="Trigger On")
        ax1.vlines(trig_end / df, i, j, color='b', lw=1, label="Trigger Off")
        ax1.vlines(trig_95_start / df, i, j, color='r', lw=2)
        ax1.vlines(trig_95_end / df, i, j, color='b', lw=2)
    except IndexError:
        pass
    fig.suptitle(trace.id)
    fig.canvas.draw()
    if show:
        plt.show()


def number_good_kurtosis_for_location(kurt_files, data_files, loc, time_dict,
                                      ugrid, snr_limit=10.0, snr_tr_limit=10.0,
                                      sn_time=10.0):
    """
    Analyses the filtered data and the kurtosis time-series to determine the
    number of stations whose traces have sufficiently high signal-to-noise
    ratio (SNR) to be useful for the location. Both time-series need to satisfy
    the conditions for a station to be counted as contributing to the location.

    :param kurt_files: Filenames for kurtosis files. Depending on the filenames
        given in this list, the function will analyse kurtosis,
        kurtosis-gradient or gaussian waveforms.
    :param data_files: Filenames for filtered data.
    :param loc: Location dictionary for the event to be analysed
    :param time_dict: Dictionary of travel-times for the location to be
        analysed
    :param snr_limit: SNR limit for kurtosis-type data
    :param snr_limit_tr: SNR limit for filtered data
    :param snr_time: Length of time in seconds before the event for computation
        of SNR.

    :rtype: integer
    :returns: Numer of stations that have contributed to the location.
    """

    x, y, z = ugrid
    o_time = loc['o_time']
    stack_x = loc['x_mean']
    stack_y = loc['y_mean']
    stack_z = loc['z_mean']

    n_good_kurt = 0
    wf = Waveform()

    for ifile in xrange(len(kurt_files)):
        kfilename = kurt_files[ifile]
        dfilename = data_files[ifile]

        st = read(kfilename, headonly=True)
        staname = st.traces[0].stats.station

        if staname in time_dict.keys():
            # approximate the traveltime by the travel-time of the closest
            # point in the irregular grid. This is a bad approximation, but
            # as it is only used for a signal to noise calculation, accuracy
            # is less important than speed
            ic, xc, yc, zc = ugrid_closest_point_index(x, y, z, stack_x,
                                                       stack_y, stack_z)
            traveltime = time_dict[staname][ic]

            start_time = o_time+traveltime-sn_time
            end_time = o_time+traveltime+sn_time
            try:
                wf.read_from_file(kfilename, starttime=start_time,
                                  endtime=end_time)
                snr = wf.get_snr(o_time+traveltime, start_time, end_time)

                wf.read_from_file(dfilename, starttime=start_time,
                                  endtime=end_time)
                snr_tr = wf.get_snr(o_time+traveltime, start_time, end_time)

                if snr > snr_limit and snr_tr > snr_tr_limit:
                    n_good_kurt = n_good_kurt + 1

            except UserWarning:
                logging.info('No data around %s for file %s.' %
                             (o_time.isoformat(), kfilename))

    return n_good_kurt


def trigger_locations_inner(max_val, max_x, max_y, max_z, left_trig,
                            right_trig, start_time, delta):
    """
    Inner loop of the location process.

    :param max_val: Time-series of stack-max.
    :param max_x: Time-series of the x-positions corresponding to stack-max.
    :param max_y: Time-series of the y-positions corresponding to stack-max.
    :param max_z: Time-series of the z-positions corresponding to stack-max.
    :param left_trig: Amplitude for trigger-on.
    :param right_trig: Amplitude for trigger-off.
    :param start_time: UTCDateTime of the first point in the time-series.
    :param delta: Sampling interval of the time-series

    :returns: List of locations. Each location is a dictionary containing all
        the relevant information.
    """

    locs = []
    trigs = trigger.triggerOnset(np.array(max_val), left_trig, right_trig)

    logging.debug('Found %d triggers.' % len(trigs))

    for trig in trigs:

        i_start = trig[0]
        i_end = trig[1]+1
        i_max_trig = np.argmax(max_val[i_start:i_end])+i_start
        max_trig = max_val[i_max_trig]
        max_trig_95 = 0.95*max_trig
        logging.debug('Max_trig = %.3f, max_trig_95 = %.3f' %
                      (max_trig, max_trig_95))
        trigs_95 = trigger.triggerOnset(max_val[i_start:i_end], max_trig_95,
                                        max_trig_95)
        for trig_95 in trigs_95:
            if i_max_trig >= trig_95[0]+i_start and \
               i_max_trig <= trig_95[1]+i_start:
                loc_dict = {}
                loc_dict['max_trig'] = max_trig
                i_start_95 = trig_95[0]+i_start
                i_end_95 = trig_95[1]+1+i_start
                loc_dict['x_mean'] = np.mean(max_x[i_start_95:i_end_95])
                loc_dict['x_sigma'] = np.std(max_x[i_start_95:i_end_95])
                loc_dict['y_mean'] = np.mean(max_y[i_start_95:i_end_95])
                loc_dict['y_sigma'] = np.std(max_y[i_start_95:i_end_95])
                loc_dict['z_mean'] = np.mean(max_z[i_start_95:i_end_95])
                loc_dict['z_sigma'] = np.std(max_z[i_start_95:i_end_95])

                loc_dict['o_time'] = start_time + i_max_trig*delta
                loc_dict['o_err_left'] = (i_max_trig-i_start_95)*delta
                loc_dict['o_err_right'] = (i_end_95-i_max_trig)*delta

                locs.append(loc_dict)

    return locs


def do_locations_trigger_setup_and_run(opdict):
    """
    Run locations trigger. Takes all options and paramters from
        WavelocOptions.opdict dictionary.

    :param opdict: Dictionary of waveloc options / parameters.
    """

    base_path = opdict['base_path']
    # parse command line
    data_dir = os.path.join(base_path, 'data', opdict['datadir'])
    kurt_files = glob.glob(os.path.join(data_dir, opdict['kurtglob']))
    data_files = glob.glob(os.path.join(data_dir, opdict['dataglob']))
    kurt_files.sort()
    data_files.sort()

    x, y, z, time_grids = get_interpolated_time_ugrids(opdict)

    logging.info("Starting log for combine_stacks.")

    out_path = os.path.join(base_path, 'out', opdict['outdir'])
    stack_path = os.path.join(out_path, 'stack')

    reloc = opdict['reloc']
    if reloc:
        loc_path = os.path.join(out_path, 'reloc')
        stack_files = glob.glob(os.path.join(stack_path,
                                             'reloc_stack_all*.hdf5'))
        stack_files.sort()
    else:
        loc_path = os.path.join(out_path, 'loc')
        stack_files = glob.glob(os.path.join(stack_path, 'stack_all*.hdf5'))
        stack_files.sort()

    n_stacks = len(stack_files)
    if n_stacks == 0:
        raise UserWarning('Empty list of stacks in %s' % (stack_path))

    loc_filename = os.path.join(loc_path, "locations.dat")
    logging.info("Path for stack files : %s" % stack_path)
    logging.info("Path for loc files : %s" % loc_path)
    logging.info("Location file : %s" % loc_filename)

    # DO DATA PREP ACCORDING TO RELOC OR NOT

    logging.info("\nDealing with continuous location, so merging stack files\
                 directly ...\n")

    # get basic info from first file
    f_stack = h5py.File(stack_files[0], 'r')
    max_val = f_stack['max_val']
    dt = max_val.attrs['dt']
    f_stack.close()

    # get start times  (get first and last times)
    start_times = []
    end_times = []
    for fname in stack_files:
        f_stack = h5py.File(fname, 'r')
        max_val = f_stack['max_val']
        start_times.append(utcdatetime.UTCDateTime(
                           max_val.attrs['start_time']))
        end_times.append(utcdatetime.UTCDateTime(
                         max_val.attrs['start_time'])+dt*len(max_val))
        f_stack.close()

    first_start_time = min(start_times)
    last_end_time = max(end_times)

    nt_full = int((last_end_time-first_start_time)/dt)+1

    # create - assume all stacks are of the same length and will be
    # concatenated end to end (this will give more than enough space)
    f = h5py.File(os.path.join(stack_path, 'combined_stack_all.hdf5'), 'w')
    cmax_val = f.create_dataset('max_val', (nt_full,), 'f', chunks=(nt_full,))
    cmax_x = f.create_dataset('max_x', (nt_full,), 'f', chunks=(nt_full,))
    cmax_y = f.create_dataset('max_y', (nt_full,), 'f', chunks=(nt_full,))
    cmax_z = f.create_dataset('max_z', (nt_full,), 'f', chunks=(nt_full,))

    # concatenate unsmoothed versions of max_val to avoid
    # problems at file starts and ends
    for i in range(n_stacks):
        f_stack = h5py.File(stack_files[i], 'r')
        max_val = f_stack['max_val']
        max_x = f_stack['max_x']
        max_y = f_stack['max_y']
        max_z = f_stack['max_z']

        # get time info for this stack
        nt = len(max_val)
        start_time = utcdatetime.UTCDateTime(max_val.attrs['start_time'])
        ibegin = np.int((start_time-first_start_time)/dt)

        # copy data over into the right place
        cmax_val[ibegin:ibegin+nt] = max_val[:]
        cmax_x[ibegin:ibegin+nt] = max_x[:]
        cmax_y[ibegin:ibegin+nt] = max_y[:]
        cmax_z[ibegin:ibegin+nt] = max_z[:]

        # close the stack
        f_stack.close()

    # create the smoothed version of the max stack
    cmax_val_smooth = f.create_dataset('max_val_smooth', (nt_full,), 'f',
                                       chunks=(nt_full,))
    cmax_val_smooth[:] = smooth(np.array(cmax_val), 51)

    for name in f:
        dset = f[name]
        dset.attrs['dt'] = dt
        dset.attrs['start_time'] = first_start_time.isoformat()

    # DO TRIGGERING AND LOCATION
    if opdict['auto_loclevel']:
        loclevel = opdict['snr_loclevel']*np.median(cmax_val_smooth)
        opdict['loclevel'] = loclevel
    else:
        loclevel = opdict['loclevel']
    left_trig = loclevel
    right_trig = loclevel

    loc_list = trigger_locations_inner(cmax_val_smooth[:], cmax_x, cmax_y,
                                       cmax_z, left_trig, right_trig,
                                       first_start_time, dt)
    logging.info('Found %d initial.' % (len(loc_list)))

    # close the stack file
    f.close()

    loc_file = open(loc_filename, 'w')
    write_header_options(loc_file, opdict)

    snr_limit = opdict['snr_limit']
    snr_tr_limit = opdict['snr_tr_limit']
    sn_time = opdict['sn_time']
    n_kurt_min = opdict['n_kurt_min']

    n_ok = 0
    locs = []
    ugrid = (x, y, z)
    for loc in loc_list:
        if number_good_kurtosis_for_location(kurt_files, data_files, loc,
                                             time_grids, ugrid, snr_limit,
                                             snr_tr_limit, sn_time) > \
                n_kurt_min:
            logging.info("Max = %.2f, %s - %.2fs + %.2f s, x=%.4f pm %.4f km,\
                         y=%.4f pm %.4f km, z=%.4f pm %.4f km" %
                         (loc['max_trig'], loc['o_time'].isoformat(),
                          loc['o_err_left'], loc['o_err_right'], loc['x_mean'],
                          loc['x_sigma'], loc['y_mean'], loc['y_sigma'],
                          loc['z_mean'], loc['z_sigma']))
            loc_file.write(u"Max = %.2f, %s - %.2f s + %.2f s, x= %.4f pm %.4f\
                           km, y= %.4f pm %.4f km, z= %.4f pm %.4f km\n" %
                           (loc['max_trig'], loc['o_time'].isoformat(),
                            loc['o_err_left'], loc['o_err_right'],
                            loc['x_mean'], loc['x_sigma'], loc['y_mean'],
                            loc['y_sigma'], loc['z_mean'], loc['z_sigma']))
            n_ok = n_ok+1
            locs.append(loc)
        else:
            logging.info("Not enough kurtosis picks for : Max = %.2f,\
                         %s - %.2fs + %.2fs, x=%.4f pm %.4f, y=%.4f pm %.4f,\
                         z=%.4f pm %.4f" % (loc['max_trig'],
                         loc['o_time'].isoformat(), loc['o_err_left'],
                         loc['o_err_right'], loc['x_mean'], loc['x_sigma'],
                         loc['y_mean'], loc['y_sigma'], loc['z_mean'],
                         loc['z_sigma']))
    loc_file.close()
    logging.info('Wrote %d locations to file %s.' % (n_ok, loc_filename))

    return locs


def read_locs_from_file(filename):
    """
    Read locations from file.

    :param filename: File to read.

    :returns: List of locations (each location is a dictionary)
    """

    from obspy.core import utcdatetime

    locs = []

    f = open(filename, 'r')
    lines = f.readlines()
    f.close()

    for line in lines:

        if not line.isspace() and line.split()[0][0] != '#':

            loc = {}

            loc['max_trig'] = np.float(line.split()[2].split(',')[0])
            loc['o_time'] = utcdatetime.UTCDateTime(line.split()[3])
            loc['o_err_left'] = np.float(line.split()[5])
            loc['o_err_right'] = np.float(line.split()[8])
            loc['x_mean'] = np.float(line.split()[11])
            loc['x_sigma'] = np.float(line.split()[13])
            loc['y_mean'] = np.float(line.split()[16])
            loc['y_sigma'] = np.float(line.split()[18])
            loc['z_mean'] = np.float(line.split()[21])
            loc['z_sigma'] = np.float(line.split()[23])

            if len(line.split()) > 25:
                loc['ml'] = np.float(line.split()[26])
                loc['ml_sigma'] = np.float(line.split()[28])

            locs.append(loc)

    return locs


def write_header_options(loc_file, opdict):
    """
    Write the information regarding how the locations were generated to an
    already opened location file. This information becomes a location "header".

    :param loc_file: File object to which the information will be written.
    :param opdict: Waveloc options dictionary
    """

    # Header of locations.dat
    loc_file.write(u'#FILTER : %.1f - %.1f Hz\n' % (opdict['c1'],
                                                    opdict['c2']))
    loc_file.write(u'#KURTOSIS = window: %.2f s, recurs: %s, grad: %s, \
                   gauss: %s\n' % (opdict['kwin'], opdict['krec'],
                                   opdict['kderiv'], opdict['gauss']))
    loc_file.write(u'#OPTIONS = reloc: %s\n' % opdict['reloc'])
    loc_file.write(u'#LOCATION = level: %d, window of analysis: %.2f s, \
                   kurtosis snr: %.2f, waveform snr: %.2f, number of \
                   stations: %d\n\n' % (opdict['loclevel'], opdict['sn_time'],
                                        opdict['snr_limit'],
                                        opdict['snr_tr_limit'],
                                        opdict['n_kurt_min']))


def read_header_from_file(filename, opdict):
    """
    Read header information from a location file into a WavelocOptions.opdict.

    :param filename: File to read.
    :param opdict: Waveloc options dictionary to write into.

    :returns: opdict
    """

    f = open(filename, 'r')
    lines = f.readlines()
    f.close()

    opdict['c1'] = np.float(lines[0].split()[2])
    opdict['c2'] = np.float(lines[0].split()[4])
    opdict['kwin'] = np.float(lines[1].split()[3])
    opdict['krec'] = lines[1].split()[6][:-1]
    opdict['kderiv'] = lines[1].split()[8][:-1]
    opdict['gauss'] = lines[1].split()[10]
    opdict['reloc'] = lines[2].split()[3]
    opdict['loclevel'] = np.int(lines[3].split()[3][:-1])
    opdict['snr_tr_limit'] = np.float(lines[3].split()[14][:-1])

    return opdict
