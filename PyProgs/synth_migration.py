import os
import glob
import h5py
import logging
import numpy as np
from hdf5_grids import H5SingleGrid


def generateSyntheticDirac(opdict, time_grids=None, ugrid=True):
    """
    Generates a synthetic test and does the migration. All options are given in
    the WavelocOptions.opdict.

    :param opdict: Waveloc options / parameters
    :param time_grids:

    """

    # Creates the synthetic dataset for us to work with

    from NllGridLib import read_stations_file, read_hdr_file
    from migration import migrate_4D_stack, extract_max_values
    from hdf5_grids import get_interpolated_time_grids,\
        get_interpolated_time_ugrids

    load_time_grids = False
    if time_grids is None:
        load_time_grids = True

    #define length and sampling frequency of synthetic data
    s_amplitude = opdict['syn_amplitude']
    s_data_length = opdict['syn_datalength']
    s_sample_freq = opdict['syn_samplefreq']
    s_filename = opdict['syn_filename']

    s_npts = int(s_data_length*s_sample_freq)
    s_delta = 1/s_sample_freq
    s_kwidth = opdict['syn_kwidth']
    s_nkwidth = int(round(s_kwidth*s_sample_freq))

    # define origin time
    s_t0 = opdict['syn_otime']

    base_path = opdict['base_path']
    test_grid_file = os.path.join(base_path, 'out', opdict['outdir'], 'grid',
                                  s_filename)
    test_stack_file = os.path.join(base_path, 'out', opdict['outdir'], 'stack',
                                   'stack_all_'+s_filename)
    test_info_file = os.path.join(base_path, 'out', opdict['outdir'], 'grid',
                                  '%s.info' % s_filename)

    # get filenames for time-grids and search grids
    search_grid_filename = os.path.join(base_path, 'lib',
                                        opdict['search_grid'])
    stations_filename = os.path.join(base_path, 'lib', opdict['stations'])
    stations = read_stations_file(stations_filename)

    if 'sta_list' in opdict:
        sta_list = opdict['sta_list'].split(',')
    else:
        sta_list = stations.keys()

    # get parameters for noise etc
    syn_addnoise = opdict['syn_addnoise']

    #################################
    # start setting up synthetic data
    #################################

    grid_info = read_hdr_file(search_grid_filename)

    if load_time_grids:
        if ugrid:
            x, y, z, time_grids = get_interpolated_time_ugrids(opdict)
        else:
            x, y, z, time_grids = get_interpolated_time_grids(opdict)
    n_buf = len(x)

    #################################
    # create synthetic data
    #################################

    it = int(round(s_t0/s_delta))

    # retrieve travel times for chosen hypocenter
    # and station list
    syn_x = opdict['syn_x']
    syn_y = opdict['syn_y']
    syn_z = opdict['syn_z']

    # get full time grid names
    full_time_grids = glob.glob(os.path.join(base_path, 'lib',
                                opdict['time_grid']+'*.hdf5'))
    full_time_grids.sort()

    ttimes = {}
    for f_timegrid in full_time_grids:
        tgrid = H5SingleGrid(filename=f_timegrid)
        sta = tgrid.grid_info['station']
        if sta in sta_list:
            ttimes[sta] = tgrid.value_at_point(syn_x, syn_y, syn_z)
        del(tgrid)

    logging.debug('Travel-times for true hypocenter = %s' % ttimes)

    # construct data with these travel times
    data = {}
    for key, delay in ttimes.iteritems():
        if syn_addnoise:
            s_snr = opdict['syn_snr']
            s = np.random.rand(s_npts)*s_amplitude/s_snr
        else:
            s = np.zeros(s_npts)
        atime = s_t0+delay
        i_atime = np.int(atime/s_delta)
        if i_atime+s_nkwidth > len(s):
            logging.error('syn_datalength is too small compared with\
                           geographical size of network ')
        s[i_atime:i_atime+s_nkwidth] = s_amplitude-np.arange(s_nkwidth) * \
            (s_amplitude/float(s_nkwidth))
        data[key] = s

    # DO MIGRATION

    logging.info('Doing migration to %s' % test_grid_file)
    f = h5py.File(test_grid_file, 'w')
    stack_grid = f.create_dataset('stack_grid', (n_buf, s_npts), 'f',
                                  chunks=(1, s_npts))
    stack_shift_time = migrate_4D_stack(data, s_delta, time_grids, stack_grid)
    n_buf, nt = stack_grid.shape

    # add useful information to dataset
    stack_grid.attrs['dt'] = s_delta
    stack_grid.attrs['nt'] = nt
    stack_grid.attrs['n_buf'] = n_buf
    stack_grid.attrs['start_time'] = -stack_shift_time
    f.create_dataset('x', data=x)
    f.create_dataset('y', data=y)
    f.create_dataset('z', data=z)

    # extract max-stack
    logging.info('Extracting max_val etc. to %s' % test_stack_file)
    f_stack = h5py.File(test_stack_file, 'w')
    # extract maxima
    grid_info = (x, y, z)
    extract_max_values(stack_grid, grid_info, f_stack)
    for name in f_stack:
        dset = f_stack[name]
        logging.debug('After extract_max_values : %s %f %f' %
                      (name, np.max(dset), np.sum(dset)))
        dset.attrs['start_time'] = -stack_shift_time
        dset.attrs['dt'] = s_delta

    # close the stack and grid files
    f_stack.close()
    f.close()
    logging.info('Saved 4D grid to file %s' % test_grid_file)

    shifted_it = it+int(round(stack_shift_time/s_delta))

    # SETUP information to pass back
    test_info = {}
    test_info['dat_file'] = test_grid_file
    test_info['stack_file'] = test_stack_file
    test_info['grid_shape'] = n_buf, nt
    test_info['dt'] = s_delta
    test_info['true_it'] = shifted_it
    test_info['true_loc'] = syn_x, syn_y, syn_z
    test_info['start_time'] = -stack_shift_time

    logging.debug(test_info)
    f = open(test_info_file, 'w')
    f.write(str(test_info))

    return test_info
