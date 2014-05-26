#!/usr/bin/env python
# encoding: utf-8

import os
import h5py
import glob
import logging
import numpy as np

from locations_trigger import read_locs_from_file, \
    do_locations_trigger_setup_and_run
from NllGridLib import read_hdr_file
from hdf5_grids import get_interpolated_time_grids
from migration import do_migration_loop_continuous
from OP_waveforms import read_data_compatible_with_time_dict
from integrate4D import compute_expected_coordinates4D, \
    compute_expected_coordinates3D


def read_prob_locs_from_file(filename):
    """
    Read file containing probability determined locations.

    :param filename: File to be read.

    :returns: Dictionary of locations
    """

    from obspy.core import utcdatetime

    locs = []

    f = open(filename, 'r')
    lines = f.readlines()
    f.close()

    for line in lines:
        loc = {}

        loc['o_time'] = utcdatetime.UTCDateTime(line.split()[5])
        loc['o_err'] = np.float(line.split()[8])
        loc['x_mean'] = np.float(line.split()[11])
        loc['x_sigma'] = np.float(line.split()[13])
        loc['y_mean'] = np.float(line.split()[16])
        loc['y_sigma'] = np.float(line.split()[18])
        loc['z_mean'] = np.float(line.split()[21])
        loc['z_sigma'] = np.float(line.split()[23])

        locs.append(loc)

    return locs


def do_locations_prob_setup_and_run(opdict):
    """
    Setup and run probability-based locations on migration grids. Takes all
    parameters from WavelocOptions.opdict.

    :param opdict: Parameters and options for Waveloc.
    """

    # get / set info
    base_path = opdict['base_path']
    space_only = opdict['probloc_spaceonly']

    locfile = os.path.join(base_path, 'out', opdict['outdir'], 'loc',
                           'locations.dat')
    locfile_prob = os.path.join(base_path, 'out', opdict['outdir'], 'loc',
                                'locations_prob.dat')
    locfile_hdf5 = os.path.join(base_path, 'out', opdict['outdir'], 'loc',
                                'locations_prob.hdf5')
    f_prob = open(locfile_prob, 'w')

    # if locfile does not exist then make it by running trigger location
    if not os.path.exists(locfile):
        logging.info('No location found at %s.  Running trigger location \
            first...' % locfile)
        do_locations_trigger_setup_and_run(opdict)

    # data files
    data_dir = os.path.join(base_path, 'data', opdict['datadir'])
    data_glob = opdict['dataglob']
    kurt_glob = opdict['kurtglob']
    grad_glob = opdict['gradglob']
    data_files = glob.glob(os.path.join(data_dir, data_glob))
    kurt_files = glob.glob(os.path.join(data_dir, kurt_glob))
    grad_files = glob.glob(os.path.join(data_dir, grad_glob))
    data_files.sort()
    kurt_files.sort()
    grad_files.sort()

    # grids
    search_grid_filename = os.path.join(base_path, 'lib',
                                        opdict['search_grid'])

    # read time grid information
    time_grids = get_interpolated_time_grids(opdict)

    # read locations
    locs = read_locs_from_file(locfile)

    # prepare file for output of marginals
    f_marginals = h5py.File(locfile_hdf5, 'w')

    # iterate over locations
    for loc in locs:
        # create the appropriate grid on the fly

        # generate the grids
        o_time = loc['o_time']
        if space_only:
            start_time = o_time
            end_time = o_time
        else:
            start_time = o_time-3*loc['o_err_left']
            end_time = o_time+3*loc['o_err_right']

        # make a buffer for migration
        start_time_migration = start_time - 10.0
        end_time_migration = end_time + 10.0

        # re-read grid info to ensure clean copy
        grid_info = read_hdr_file(search_grid_filename)

        # read data
        grad_dict, delta = \
            read_data_compatible_with_time_dict(grad_files, time_grids,
                                                start_time_migration,
                                                end_time_migration)

        # do migration (all metadata on grid is added to grid_info)
        do_migration_loop_continuous(opdict, grad_dict, delta,
                                     start_time_migration, grid_info,
                                     time_grids, keep_grid=True)

        # integrate to get the marginal probability density distributions

        # get required info
        grid_starttime = grid_info['start_time']
        nx, ny, nz, nt = grid_info['grid_shape']
        dx, dy, dz, dt = grid_info['grid_spacing']
        x_orig, y_orig, z_orig = grid_info['grid_orig']

        # we are only interested in the time around the origin time of the
        # event
        it_left = np.int(np.round((start_time - grid_starttime)/dt))
        it_right = np.int(np.round((end_time - grid_starttime)/dt))
        it_true = np.int(np.round((o_time - grid_starttime)/dt))
        nt = (it_right-it_left)+1

        # set up integration axes (wrt reference)
        x = np.arange(nx)*dx
        y = np.arange(ny)*dy
        z = np.arange(nz)*dz
        if not space_only:
            t = np.arange(nt)*dt

        # open the grid file
        grid_filename = grid_info['dat_file']
        f = h5py.File(grid_filename, 'r')
        stack_grid = f['stack_grid']

        # extract the portion of interest (copy data)
        if space_only:
            stack_3D = np.empty((nx, ny, nz))
            stack_3D[:] = stack_grid[:, it_true].reshape(nx, ny, nz)
        else:
            stack_4D = np.empty((nx, ny, nz, nt))
            stack_4D[:] = stack_grid[:, it_left:it_right+1].reshape(nx, ny,
                                                                    nz, nt)

        # close the grid file
        f.close()

        # Get expected values (normalizes grid internally)
        if space_only:
            exp_x, exp_y, exp_z, cov_matrix, prob_dict = \
                compute_expected_coordinates3D(stack_3D, x, y, z,
                                               return_2Dgrids=True)
        else:
            exp_x, exp_y, exp_z, exp_t, cov_matrix, prob_dict = \
                compute_expected_coordinates4D(stack_4D, x, y, z, t,
                                               return_2Dgrids=True)

        # put reference location back
        exp_x = exp_x + x_orig
        exp_y = exp_y + y_orig
        exp_z = exp_z + z_orig
        if space_only:
            exp_t = o_time
        else:
            exp_t = start_time + exp_t

        # extract uncertainties from covariance matrix
        if space_only:
            sig_x, sig_y, sig_z = np.sqrt(np.diagonal(cov_matrix))
            sig_t = (loc['o_err_left']+loc['o_err_right'])/2.
        else:
            sig_x, sig_y, sig_z, sig_t = np.sqrt(np.diagonal(cov_matrix))

        # save the marginals to a hdf5 file in loc subdirectory (f_marginals)
        # each event becomes a group in this one file
        grp = f_marginals.create_group(exp_t.isoformat())
        grp.create_dataset('x', data=x+x_orig)
        grp.create_dataset('y', data=y+y_orig)
        grp.create_dataset('z', data=z+z_orig)
        grp.create_dataset('prob_x', data=prob_dict['prob_x0'])
        grp.create_dataset('prob_y', data=prob_dict['prob_x1'])
        grp.create_dataset('prob_z', data=prob_dict['prob_x2'])
        grp.create_dataset('prob_xy', data=prob_dict['prob_x0_x1'])
        grp.create_dataset('prob_xz', data=prob_dict['prob_x0_x2'])
        grp.create_dataset('prob_yz', data=prob_dict['prob_x1_x2'])
        if not space_only:
            grp.create_dataset('t', data=t-(o_time - start_time))
            grp.create_dataset('prob_t', data=prob_dict['prob_x3'])
            grp.create_dataset('prob_xt', data=prob_dict['prob_x0_x3'])
            grp.create_dataset('prob_yt', data=prob_dict['prob_x1_x3'])
            grp.create_dataset('prob_zt', data=prob_dict['prob_x2_x3'])

        # write the expected values to a plain text locations file
        f_prob.write("PROB DENSITY : T = %s s pm %.2f s, x= %.4f pm %.4f km, \
                     y= %.4f pm %.4f km, z= %.4f pm %.4f km\n" %
                     (exp_t.isoformat(), sig_t, exp_x, sig_x, exp_y, sig_y,
                      exp_z, sig_z))

    # close location files
    f_prob.close()
    f_marginals.close()


if __name__ == '__main__':

    from options import WavelocOptions
    logging.basicConfig(level=logging.INFO,
                        format='%(levelname)s : %(asctime)s : %(message)s')

    wo = WavelocOptions()
    args = wo.p.parse_args()

    wo.set_all_arguments(args)
    wo.verify_location_options()

    do_locations_prob_setup_and_run(wo.opdict)
