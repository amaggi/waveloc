#!/usr/bin/env python
# encoding: utf-8

"""
Routines to do double-difference relative location of events.
"""

import os
import numpy as np
import logging

from obspy.core import utcdatetime

from locations_trigger import read_locs_from_file, read_header_from_file, \
    write_header_options
from correlation import BinaryFile
from NllGridLib import read_stations_file, read_hdr_file
from hdf5_grids import get_interpolated_time_grids


def traveltimes(x, y, z, t_orig, stations, time_grids):
    """
    Compute theoretical traveltimes and arrival times for a point to all
    stations using pre-calculated time-grids

    :param x: x-coordinate of point
    :param y: y-coordinate of point
    :param z: z-coordinate of point
    :param t_orig: origin time of envent (used to calculate arrival times)
    :param stations: list of station names
    :param time_grids: dictionary of time-grids

    :type x: float
    :type y: float
    :type z: float
    :type t_orig: UTCDateTime
    :type stations: string

    :rtype: dictionary
    :returns: t_th, arr_times

    """

    t_th = {}
    arr_times = {}

    for staname in sorted(stations):
        if not staname in time_grids.keys():
            logging.info("%s station not in time_grids" % staname)
            continue
        t_th[staname] = []
        arr_times[staname] = []
        for i in range(len(x)):
            # travel-time
            t_th[staname].append(time_grids[staname].value_at_point(x[i], y[i],
                                 z[i]))
            # arrival-time
            arr_times[staname].append(utcdatetime.UTCDateTime(t_orig[i]) +
                                      t_th[staname][i])

    return t_th, arr_times


def partial_deriv(coord, ev, tth):
    """
    Computes partial derivatives

    :param coord: list of coordinates [x,y,z] of a given station
    :param ev: list of coordinates [x,y,z] of the event
    :param tth: theoretical traveltime from the event to the station

    :type coord: list
    :type ev: list
    :type tth: float

    :rtype: list
    :returns: dpx, dpy, dpz

    """

    norm = (coord[0]-ev[0])**2+(coord[1]-ev[1])**2+(coord[2]-ev[2])**2
    dpx = -(coord[0]-ev[0])*tth/norm
    dpy = -(coord[1]-ev[1])*tth/norm
    dpz = -(coord[2]-ev[2])*tth/norm

    return [dpx, dpy, dpz]


def fill_matrix(cluster, x, y, z, t_orig, stations, t_th, t_arr, coeff, delay,
                threshold):
    """
    Fills in matrix G (partial derivatives) data-vector d (double differences)
    and matrix W (weights)

    :param cluster: indices of events in the cluster
    :param x: x-coordinates of all events
    :param y: y-coordinates of all events
    :param z: z-coordinates of all events
    :param t_orig: origin times of all events
    :param stations: dictionary of station coordinates
    :param t_th: dictionary of theoretical traveltimes
    :param t_arr: dictionary of theoretical arrival times
    :param coeff: cross-correlation coefficients between all possible pairs of events
    :param delay: time delays measured between all possible pairs of events
    :param threshold: 

    :type cluster: list
    :type x: list
    :type y: list
    :type z: list
    :type t_orig: list
    :type stations: dictionary
    :type t_th: dictionary
    :type t_arr: dictionary
    :type coeff: dictionary
    :type delay: dictionary
    :type threshold: float
    """
    G, W, d = [], [], []
    N = len(cluster)
    nline, num = 0, 0

    for staname in sorted(stations):
        if not staname in delay.keys():
            continue
        if not staname in t_th.keys():
            continue
        coord = [stations[staname]['x'], stations[staname]['y'],
                 -stations[staname]['elev']]
        for n in range(N):
            ev1 = [x[n], y[n], z[n]]
            e1 = cluster[n]
            dp1 = partial_deriv(coord, ev1, t_th[staname][n])

            for nn in range(n+1, N):
                e2 = cluster[nn]
                if delay[staname][e1-1][e2-1] != 'NaN' and \
                   coeff[staname][e1-1][e2-1] >= threshold:
                    # fill G
                    G.append(np.zeros(4*N))
                    ev2 = [x[nn], y[nn], z[nn]]
                    dp2 = partial_deriv(coord, ev2, t_th[staname][nn])
                    dp2 = [-elt for elt in dp2]

                    G[nline][4*n:4*n+4] = dp1+[1]
                    G[nline][4*nn:4*nn+4] = dp2+[-1]
                    nline += 1

                    # fill d
                    obs = delay[staname][e1-1][e2-1]+(t_orig[n]-t_orig[nn])
                    theo = t_arr[staname][n]-t_arr[staname][nn]
                    d.append(obs-theo)

                    # fill W
                    W.append(coeff[staname][e1-1][e2-1])

        num += 1

    return G, d, W


def centroid_constraint(G, d, W):
    """
    Applies centroid constraint: sum(delta m)=0.

    :param G: matrix of partial derivatives
    :param d: vector of double-difference times
    :param W: weighting matrix

    :returns: Updated versions of G, d, W

    """

    ncol = np.size(G, 1)
    nline = np.size(G, 0)
    for i in range(4):
        G.append(np.zeros(ncol))
        G[nline][i::4] = 1
        W.append(0.5)
        nline += 1
    d.extend(np.zeros(4))
    W = np.diag(W)

    return np.matrix(G), np.transpose(np.matrix(d)), np.matrix(W)


def inversion(G, d, W):
    """
    Inverts G, d, W using LSQR to obtain parameter vector m

    :param G: matrix of partial derivatives
    :param d: vector of double-difference times
    :param W: weighting matrix

    :returns: Parameter vector m

    """

    Gt = np.transpose(G)
    Winv = W.getI()
    a = Gt*Winv*G
    ainv = a.getI()
    m = ainv*Gt*Winv*d

    return m


def coord_cluster(cluster, locs):
    """
    Extract the coordinates of the events of a given cluster

    :param cluster: indices of events in the cluster
    :param locs: list of the whole locations (each element of the list is a dictionary)

    :type cluster: list
    :type locs: list

    :rtype: list
    :returns: xini, yini, zini, zini_ph, to_ini

    """

    xini, yini, zini, zini_ph, to_ini = [], [], [], [], []
    for ind in cluster:
        xini.append(locs[ind-1]['x_mean'])
        yini.append(locs[ind-1]['y_mean'])
        zini.append(locs[ind-1]['z_mean'])
        zini_ph.append(-locs[ind-1]['z_mean'])  # positive z axis upwards
        to_ini.append(locs[ind-1]['o_time'])

    return xini, yini, zini, zini_ph, to_ini


def plot_events(cluster, locs, stations, x, y, z, i, threshold, nbmin, area,
                nbsta):
    """
    Plot old and new locations (uses mayavi)

    :param cluster: indices of events composing all clusters
    :param locs: list of the whole locations (each element of the list is a dictionary)
    :param stations: dictionary of stations
    :param x: new x-coordinates of events
    :param y: new y-coordinates of events
    :param z: new z-coordinates of events
    :param i: cluster index
    :param threshold: minimum value of cross-correlation coefficient used to form a cluster
    :param nbmin: minimum number of stations required to form a cluster
    :param area: coordinates of the study area
    :param nbsta: number of stations where the measured correlation coefficient was greater than the given threshold for all possible event pairs

    :type cluster: numpy array
    :type locs: list
    :type stations: dictionary
    :type x: list
    :type y: list
    :type z: list
    :type i: int
    :type threshold: float
    :type nbmin: int
    :type area: list
    :type nbsta: numpy array
    """
    from mayavi import mlab
    from CZ_color import CZ_W_2_color

    # Stations coordinates
    xsta, ysta, zsta = [], [], []
    for sta in sorted(stations):
        xsta.append(stations[sta]['x'])
        ysta.append(stations[sta]['y'])
        zsta.append(stations[sta]['elev'])

    z_ph = [-elt for elt in z]

    # Initial hypocentral parameters
    xini, yini, zini, zini_ph, to_ini = coord_cluster(cluster[i], locs)

    mlab.figure(i, bgcolor=(1, 1, 1), fgcolor=(0, 0, 0), size=(1000, 900))
    mlab.clf()
    # yellow : initial locations
    mlab.points3d(xini, yini, zini_ph, color=(1, 1, 0), scale_factor=0.2)
    mlab.points3d(xsta, ysta, zsta, color=(1, 0, 0), scale_factor=0.05,
                  mode='cube')
    # cyan : new locations
    mlab.points3d(x, y, z_ph, color=(0, 1, 1), scale_factor=0.2)
    mlab.axes(extent=area, color=(0, 0, 0))  # axe des z positif vers le haut
    mlab.outline(extent=area, color=(0, 0, 0))
    mlab.title("cluster=%s,  threshold=%s,  nbmin=%s" % (i, threshold, nbmin),
               height=0.1, size=0.35, color=(0, 0, 0))

    if len(cluster[i]) < 20:
        for ind_I in range(len(cluster[i])):
            for ind_J in range(ind_I+1, len(cluster[i])):
                ev_I = cluster[i][ind_I]-1
                ev_J = cluster[i][ind_J]-1
                W_IJ = nbsta[ev_I, ev_J]
                if W_IJ >= nbmin:
                    mlab.points3d(xini[ind_J], yini[ind_J], zini_ph[ind_J],
                                  scale_factor=0.1, color=(0, 0, 0))
                    mlab.points3d(xini[ind_I], yini[ind_I], zini_ph[ind_I],
                                  scale_factor=0.1, color=(0, 0, 0))
                    d = (xini[ind_J]-xini[ind_I], yini[ind_J]-yini[ind_I],
                         zini_ph[ind_J]-zini_ph[ind_I])
                    norm = np.sqrt(d[0]**2+d[1]**2+d[2]**2)
                    mlab.quiver3d(xini[ind_I], yini[ind_I], zini_ph[ind_I],
                                  d[0], d[1], d[2],
                                  color=tuple(CZ_W_2_color(W_IJ)),
                                  mode='2ddash', scale_factor=norm,
                                  scale_mode='scalar')

    mlab.show()


def do_double_diff(x, y, z, to, stations, coeff, delay, cluster, threshold,
                   t_th,  arr_times):
    """
    Do double difference location (inner routine) and return new coordinates. 

    :param x: x-coordinates of events of a given cluster
    :param y: y-coordinates of events of a given cluster
    :param z: z-coordinates of events of a given cluster
    :param to: origin times of events of a given cluster
    :param stations: dictionary of stations
    :param coeff: cross-correlation coefficients between all possible pairs of events
    :param delay: time delays measured between all possible pairs of events
    :param cluster: indices of events in the cluster
    :param threshold: minimum value of cross-correlation coefficient used to form a cluster
    :param t_th: theoretical traveltimes
    :param arr_times: theoretical arrival times

    :type x: list
    :type y: list
    :type z: list
    :type to: list
    :type stations: dictionary
    :type coeff: dictionary
    :type delay: dictionary
    :type cluster: list
    :type threshold: float
    :type t_th: dictionary
    :type arr_times: dictionary

    :rtype: list
    :returns: x, y, z, to
    """
    N = len(cluster)

    # Fill G, d and W
    G,  d,  W = fill_matrix(cluster, x, y, z, to, stations, t_th, arr_times,
                            coeff, delay, threshold)

    # Centroid constraint : add 4 lines to G, d and W
    G, d, W = centroid_constraint(G, d, W)

    # Inversion
    m = inversion(G, d, W)

    for i in range(N):
        x[i] = x[i]+m[4*i, 0]
        y[i] = y[i]+m[4*i+1, 0]
        z[i] = z[i]+m[4*i+2, 0]
        to[i] = utcdatetime.UTCDateTime(to[i])+m[4*i+3, 0]

    return x, y, z, to


def do_double_diff_setup_and_run(opdict):
    """
    Do double difference (outer routine). Takes options from a
    WavelocOptions.opdict dictionary.

    :param opdict: Dictionary of parameters and options
    """

    base_path = opdict['base_path']
    verbose = opdict['verbose']
    dd_loc = opdict['dd_loc']

    # Station
    stations_filename = os.path.join(base_path, 'lib', opdict['stations'])
    stations = read_stations_file(stations_filename)

    # Location file
    locdir = os.path.join(base_path, 'out', opdict['outdir'], 'loc')
    loc_filename = os.path.join(locdir, 'locations.dat')
    locs = read_locs_from_file(loc_filename)
    opdict = read_header_from_file(loc_filename, opdict)

    # ------------------------------------------------------------------------
    # search grid
    search_grid_filename = os.path.join(base_path, 'lib',
                                        opdict['search_grid'])
    # traveltimes grid
    grid_info = read_hdr_file(search_grid_filename)
    time_grids = get_interpolated_time_grids(opdict)

    # Extract the UTM coordinates of the area of study
    xstart = grid_info['x_orig']
    xend = xstart+grid_info['nx']*grid_info['dx']
    ystart = grid_info['y_orig']
    yend = ystart+grid_info['ny']*grid_info['dy']
    zend = -grid_info['z_orig']
    zstart = -(-zend+grid_info['nz']*grid_info['dz'])
    area = [xstart, xend, ystart, yend, zstart, zend]

    # ------------------------------------------------------------------------
    nbmin = int(opdict['nbsta'])
    threshold = float(opdict['clus'])

    # Correlation,  time delay and cluster files
    corr_file = os.path.join(locdir, opdict['xcorr_corr'])
    cfile = BinaryFile(corr_file)
    coeff = cfile.read_binary_file()

    delay_file = os.path.join(locdir, opdict['xcorr_delay'])
    dfile = BinaryFile(delay_file)
    delay = dfile.read_binary_file()

    cluster_file = os.path.join(locdir, 'cluster-%s-%s' % (str(threshold),
                                                           str(nbmin)))
    clfile = BinaryFile(cluster_file)
    cluster = clfile.read_binary_file()

    # ------------------------------------------------------------------------
    # Input parameters
    len_cluster_min = 2

    if dd_loc:
        new_loc_filename = os.path.join(locdir, 'relocations.dat')
        new_loc_file = open(new_loc_filename, 'w')
        write_header_options(new_loc_file, opdict)

    # ------------------------------------------------------------------------
    # Iterate over clusters
    for i in cluster.keys():
        print "CLUSTER %d:" % i, cluster[i], len(cluster[i])
        N = len(cluster[i])

        # Hypocentral parameters to be changed
        x, y, z, z_ph, to = coord_cluster(cluster[i], locs)

        # Replace bad locations by the centroid coordinates
        centroid_x = np.mean(x)
        centroid_y = np.mean(y)
        centroid_z = np.mean(z)

        for ii in range(len(cluster[i])):
            if np.abs(x[ii]-centroid_x) > .75:
                x[ii] = centroid_x
            if np.abs(y[ii]-centroid_y) > .75:
                y[ii] = centroid_y
            if np.abs(z[ii]-centroid_z) > .75:
                z[ii] = centroid_z

        if N > len_cluster_min:
            # Theroretical traveltimes and arrival times
            t_th, arr_times = traveltimes(x, y, z, to, stations, time_grids)
            # do double difference location
            x, y, z, to = do_double_diff(x, y, z, to, stations, coeff, delay,
                                         cluster[i], threshold, t_th,
                                         arr_times)

        if verbose:
            from clustering import compute_nbsta
            nbsta = compute_nbsta(len(locs), coeff, threshold)
            plot_events(cluster, locs, stations, x, y, z, i, threshold, nbmin,
                        area, nbsta)

        if dd_loc:
            ind = 0
        for j in cluster[i]:
            locs[j-1]['x_mean'] = x[ind]
            locs[j-1]['y_mean'] = y[ind]
            locs[j-1]['z_mean'] = z[ind]
            locs[j-1]['o_time'] = to[ind]
            locs[j-1]['x_sigma'] = 0
            locs[j-1]['y_sigma'] = 0
            locs[j-1]['z_sigma'] = 0
            locs[j-1]['o_err_right'] = 0
            locs[j-1]['o_err_left'] = 0
            ind += 1
            new_loc_file.write("Max = %.2f, %s - %.2f s + %.2f s, x= %.4f pm\
                %.4f km, y= %.4f pm %.4f km, z= %.4f pm %.4f km\n" %
                               (locs[j-1]['max_trig'],
                                locs[j-1]['o_time'].isoformat(),
                                locs[j-1]['o_err_left'],
                                locs[j-1]['o_err_right'], locs[j-1]['x_mean'],
                                locs[j-1]['x_sigma'], locs[j-1]['y_mean'],
                                locs[j-1]['y_sigma'], locs[j-1]['z_mean'],
                                locs[j-1]['z_sigma']))

    if dd_loc:
        new_loc_file.close()
