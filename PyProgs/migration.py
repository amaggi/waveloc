#!/usr/bin/env python
# encoding: utf-8

import os
import h5py
import logging
import glob
import numpy as np

from obspy.core import utcdatetime
from itertools import count, islice
from time import time

from OP_waveforms import read_data_compatible_with_time_dict
from hdf5_grids import get_interpolated_time_ugrids
from filters import smooth


def do_migration_setup_and_run(opdict):
    """
    Do setup and launch migration.

    :param opdict: WavelocOptions.opdict
    """

    base_path = opdict['base_path']
    runtime = opdict['time']
    reloc = opdict['reloc']

    # data
    data_dir = os.path.join(base_path, 'data', opdict['datadir'])
    if opdict['kderiv']:
        data_glob = opdict['gradglob']
        if opdict['gauss']:
            data_glob = opdict['gaussglob']
    else:
        data_glob = opdict['kurtglob']
    data_files = glob.glob(os.path.join(data_dir, data_glob))
    data_files.sort()
    if len(data_files) == 0:
        logging.error('No data files found for %s and %s' % (data_dir,
                                                             data_glob))
        raise UserWarning

    # grids
    x, y, z, time_grids = get_interpolated_time_ugrids(opdict)

    # start and end times
    starttime = opdict['starttime']
    endtime = opdict['endtime']
    data_length = opdict['data_length']
    data_overlap = opdict['data_overlap']

    initial_start_time = utcdatetime.UTCDateTime(starttime)
    initial_end_time = initial_start_time+data_length

    final_end_time = utcdatetime.UTCDateTime(endtime)

    time_shift_secs = data_length-data_overlap

    # ######## FOR EACH TIME SPAN - DO MIGRATION #############

    # start loop over time
    start_time = initial_start_time
    end_time = initial_end_time

    if runtime:
        t_ref = time()

    while (start_time < final_end_time):

        # read data
        logging.info("Reading data  : %s - %s." % (start_time.isoformat(),
                                                   end_time.isoformat()))
        data, delta = \
            read_data_compatible_with_time_dict(data_files, time_grids,
                                                start_time, end_time)

        if reloc:
            tr_glob = opdict['kurtglob']
            files = glob.glob(os.path.join(data_dir, tr_glob))
            traces, delta = \
                read_data_compatible_with_time_dict(files, time_grids,
                                                    start_time, end_time)
            sta_list = sorted(traces)
            for staname in sta_list:
                snr = np.max(traces[staname])/np.mean(np.abs(traces[staname]))
                if snr < opdict['reloc_snr']:
                    data[staname] = np.zeros(len(data[staname]))

        # re-read grid_info at each iteration to make sure it is a clean copy
        grid_info = (x, y, z)

        # do migration if have enough data (3 is bare minimum)
        if len(data.keys()) >= 3:
            logging.info("Migrating data : %s - %s." % (start_time.isoformat(),
                                                        end_time.isoformat()))
            do_migration_loop_continuous(opdict, data, delta, start_time,
                                         grid_info, time_grids)
        elif len(data.keys()) == 0:
            logging.warn('No data found between %s and %s.' %
                         (start_time.isoformat(), end_time.isoformat()))
        else:
            logging.warn('Insufficient data found between %s and %s.' %
                         (start_time.isoformat(), end_time.isoformat()))

        # Reset the start and end times to loop again
        start_time = start_time+time_shift_secs
        end_time = end_time+time_shift_secs

    if runtime:
        t = time()-t_ref
        logging.info("Time for migrating all time slices : %.2f s\n" % (t))


def do_migration_loop_continuous(opdict, data, delta, start_time, grid_info,
                                 time_grids, keep_grid=False,
                                 keep_stacks=True):
    """
    Do continuous migration loop. TODO : flesh out this doc-string.
    """
    logging.info("Processing time slice %s" % start_time.isoformat())

    options_time = opdict['time']
    use_ram = opdict['use_ram']
    output_dir = os.path.join(opdict['base_path'], 'out', opdict['outdir'])
    options_reloc = opdict['reloc']

    x, y, z = grid_info
    n_buf = len(x)
    min_npts = min([len(data[key]) for key in data.keys()])

    if options_time:
        t_ref = time()

    # open hdf5 file for stack_grid
    if options_reloc:
        grid_filename = os.path.join(output_dir, 'grid',
                                     'reloc_stack_grid_%s.hdf5' % start_time)
    else:
        grid_filename = os.path.join(output_dir, 'grid',
                                     'stack_grid_%s.hdf5' % start_time)
    logging.info('Creating grid file %s' % grid_filename)
    # if running on small memory machine write to hdf5 file
    if not use_ram:
        f = h5py.File(grid_filename, 'w')
        stack_grid = f.create_dataset('migrated_grid', (n_buf, min_npts), 'f',
                                      chunks=(1, min_npts))
        stack_grid[...] = 0.
        f.create_dataset('x', data=x)
        f.create_dataset('y', data=y)
        f.create_dataset('z', data=z)
    else:
        # if running on big memory machine work in ram
        stack_grid = np.empty((n_buf, min_npts), dtype='float32', order='F')
        stack_grid[...] = 0.

    # DO MIGRATION
    stack_shift_time = migrate_4D_stack(data, delta, time_grids, stack_grid,
                                        use_ram)
    stack_start_time = start_time-stack_shift_time
    n_buf, nt = stack_grid.shape

    if options_time:
        t = time()-t_ref
        logging.info("Time for migrating %d stacks, each of extent %d points :\
                     %.2f s\n" % (n_buf, nt, t))

    if keep_stacks:
        if options_time:
            t_ref = time()

        if options_reloc:
            stack_filename = os.path.join(output_dir, 'stack',
                                          'reloc_stack_all_%s.hdf5' %
                                          start_time)
        else:
            stack_filename = os.path.join(output_dir, 'stack',
                                          'stack_all_%s.hdf5' % start_time)
        logging.info('Extracting max_val etc. to %s' % stack_filename)
        f_stack = h5py.File(stack_filename, 'w')
        # extract maxima
        extract_max_values(stack_grid, grid_info, f_stack,  use_ram=use_ram)
        for name in f_stack:
            dset = f_stack[name]
            logging.debug('After extract_max_values : %s %f %f' %
                          (name, np.max(dset), np.sum(dset)))
            dset.attrs['start_time'] = stack_start_time.isoformat()
            dset.attrs['dt'] = delta

        f_stack.close()
        # keep the stack_filename
        if options_time:
            t = time()-t_ref
            logging.info("Time for extracting maxima : %.2f s\n" % (t))

    if keep_grid:
        if not use_ram:
            # working on small-ram machine
            sg = stack_grid
        else:
            # working on big-ram machine
            # create the hdf5 file
            f = h5py.File(grid_filename, 'w')
            sg = f.create_dataset('migrated_grid', data=stack_grid)
            f.create_dataset('x', data=x)
            f.create_dataset('y', data=y)
            f.create_dataset('z', data=z)
        # add useful attributes to the hdf5 dataset
        sg.attrs['n_buf'] = n_buf
        sg.attrs['nt'] = nt
        sg.attrs['dt'] = delta
        sg.attrs['start_time'] = stack_start_time.isoformat()

    # close the hdf5 file for the grid
    if (not use_ram) or keep_grid:
        f.close()
        return grid_filename, n_buf, nt, stack_shift_time

    # remove the grid file unless you want to keep it
    if (not keep_grid) and (not use_ram):
        logging.info('Removing grid file %s' % grid_filename)
        os.remove(grid_filename)


def migrate_4D_stack(data, delta, time_grids, stack_grid, use_ram=False):
    """
    Migrate the 4D stack. Central part of the waveloc process.
    TODO : flesh out this doc-string

    :param data:
    :param delta:
    :param time_grids:
    :param stack_grid:
    :param use_ram:
    """
    # save the list of data keys
    # note : keys of data are all included in keys of time_grid, but there may
    # be more times than data
    wf_ids = data.keys()
    n_wf_ids = len(wf_ids)

    n_buf, min_npts = stack_grid.shape
    logging.debug("Stack max dimension = %d x %d" % (n_buf, min_npts))

    # initialize the arrays we will need
    tmp_stack = np.empty(min_npts, dtype='float32')
    tmp_stack[:] = 0.
    i_times = np.zeros((n_wf_ids, n_buf), dtype='int')
    i_max_times = np.zeros(n_buf, dtype='int')
    i_min_times = np.zeros(n_buf, dtype='int')
    start_index = np.zeros(n_buf, dtype='int')
    start_indexes = np.zeros((n_wf_ids, n_buf), dtype='int')
    end_indexes = np.zeros((n_wf_ids, n_buf), dtype='int')
    n_lens = np.zeros((n_wf_ids, n_buf), dtype='int')

    # construct grid (n_buf x n_sta) grid of time_indexes for migration
    for i in islice(count(0), n_wf_ids):
        wf_id = wf_ids[i]
        i_times[i, :] = np.round(time_grids[wf_id][:]/delta)/1

    # find the min and max time indexes for point
    i_min_times = np.min(i_times, 0)
    i_max_times = np.max(i_times, 0)
    iextreme_min_time = np.min(i_min_times)
    iextreme_max_time = np.max(i_max_times)
    start_index = i_min_times-iextreme_min_time
    stack_shift_time = delta*iextreme_min_time

    # find start indexes, end indexes and lengths for each station and point
    start_indexes = i_times-i_min_times
    end_indexes = i_times-i_max_times+min_npts
    n_lens = end_indexes-start_indexes
    # keep the shortest length for each point
    n_len = np.min(n_lens, 0)
    # keep the shortest overall length
    shortest_n_len = np.min(n_len)

    # sill fix the length of the stack to the shortest possible length given
    # all the previous travel time information
    norm_stack_len = shortest_n_len-iextreme_max_time
    if norm_stack_len < 0:
        logging.error('Data length too short for coherent migration across\
                      network')

    # the actual migration loop
    # cannot seem to vectorize this any more... too bad !!

    for ib in islice(count(0), n_buf):

        # This is ugly, but necessary to avoid memory leak from inner loop
        _do_stack(ib, n_wf_ids, wf_ids, stack_grid, data, min_npts, n_lens,
                  start_indexes, end_indexes, start_index, norm_stack_len)

    # clean up what is no longer needed
    del i_times, i_min_times, i_max_times, start_indexes, end_indexes
    del n_lens, start_index

    # resize stack_grid
    if use_ram:
        stack_grid.resize((n_buf, norm_stack_len), refcheck=False)
    else:
        stack_grid.resize(norm_stack_len, axis=1)

    # end
    return stack_shift_time


def _do_stack(ib, n_wf_ids, wf_ids, stack_grid, data, min_npts, n_lens,
              start_indexes, end_indexes, start_index, norm_stack_len):

    tmp_stack = np.empty(min_npts, dtype='float32')
    tmp_stack[:] = 0.
    # stack shifted data from each station
    for i in islice(count(0), n_wf_ids):
        tmp_stack[0:n_lens[i, ib]] += \
            data[wf_ids[i]][start_indexes[i, ib]:end_indexes[i, ib]]

    # We need to homogenize, and get everything to start and end at the same
    # time
    stack_grid[ib, 0:norm_stack_len] = \
        tmp_stack[start_index[ib]:start_index[ib]+norm_stack_len]

    # cleanup
    del tmp_stack


def extract_max_values(stack_grid, search_info, f_stack, use_ram=False,
                       n_max=5e8):
    """
    Extracts maximum stack value from a 4D migrated grid. Also extracts x_max,
    y_max, z_max. TODO : flesh out this doc-string.

    :param stack_grid:
    :param search_info:
    :param f_stack:
    :param use_fam:
    :param n_max:
    """

    # get basic info
    x, y, z = search_info

    nb, nt = stack_grid.shape

    if use_ram:
        max_val = np.empty(nt, dtype='float32')
        max_val_smooth = np.empty(nt, dtype='float32')
        max_x = np.empty(nt, dtype='float32')
        max_y = np.empty(nt, dtype='float32')
        max_z = np.empty(nt, dtype='float32')
        # create temporary datasets
        max_ib = np.empty(nt, dtype=int)
    else:
        max_val = f_stack.create_dataset('max_val', (nt, ), 'f')
        max_val_smooth = f_stack.create_dataset('max_val_smooth', (nt, ), 'f')
        max_x = f_stack.create_dataset('max_x', (nt, ), 'f')
        max_y = f_stack.create_dataset('max_y', (nt, ), 'f')
        max_z = f_stack.create_dataset('max_z', (nt, ), 'f')
        # create temporary datasets
        max_ib = f_stack.create_dataset('max_ib', (nt, ), 'i')

    # extract values
    dt = int(n_max/nb)
    if nt <= dt:
        # do the extraction in one step
        max_ib[:] = np.argmax(stack_grid, 0)
        max_val[:] = np.max(stack_grid, 0)

    else:
        # do the extraction in steps
        n = nt/dt
        logging.info('Number of values exceeds %d. Doing extraction in %d\
                      steps' % (n_max, n))
        for i in islice(count(0), n):
            max_ib[i*dt:(i+1)*dt] = np.argmax(stack_grid[:, i*dt:(i+1)*dt], 0)
            max_val[i*dt:(i+1)*dt] = np.max(stack_grid[:, i*dt:(i+1)*dt], 0)
        max_ib[n*dt:nt] = np.argmax(stack_grid[:, n*dt:nt], 0)
        max_val[n*dt:nt] = np.max(stack_grid[:, n*dt:nt], 0)

    # find the corresponding x,y,z values
    np_max_ib = max_ib[:]   # copy data to a numpy array for indexing
    max_x[:] = x[np_max_ib]
    max_y[:] = y[np_max_ib]
    max_z[:] = z[np_max_ib]
    max_val_smooth[:] = smooth(np.array(max_val), 51)

    if use_ram:
        # add datasets to hdf5file
        max_val = f_stack.create_dataset('max_val', data=max_val)
        max_val_smooth = f_stack.create_dataset('max_val_smooth',
                                                data=max_val_smooth)
        max_x = f_stack.create_dataset('max_x', data=max_x)
        max_y = f_stack.create_dataset('max_y', data=max_y)
        max_z = f_stack.create_dataset('max_z', data=max_z)
    else:
        # clean up temporary datasets
        del f_stack['max_ib']
