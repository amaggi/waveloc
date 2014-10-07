#!/usr/bin/env python
# encoding: utf-8
"""
Provides classes and functions for clustering of earthquakes based on the
depth-first algorithm.
"""

import os
import glob
import matplotlib.pyplot as plt
import numpy as np
from CZ_color import CZ_Clust_2_color, CZ_W_2_color
from OP_waveforms import Waveform
import logging
from correlation import BinaryFile
from locations_trigger import read_locs_from_file
from NllGridLib import read_stations_file


class Graph(object):
    """
    Class of Graph objects. Contains 4 attributes:

    ** Attributes **

    .. attribute:: flag

        Indicates if an event already belongs to a cluster (1) or not (0).

    .. attribute:: cluster_index

        Indicates the cluster number of all events. 0 means the event does not
        belong to any cluster.

    .. attribute:: neighbours

        Lists the indexes of the neighbour events for each event.
        GRAPH.neighbours is then a list of lists.

    """

    def __init__(self):
        """
        Initialises all attributes to empty lists []
        """
        self.flag = []
        self.cluster_index = []
        self.neighbours = []

    def set_flag(self, value):
        """
        Appends value to flag attribute.

        :param value: Value to append (0=False, 1=True).
        :type value: integer
        """
        self.flag.append(value)

    def set_cluster_index(self, value):
        """
        Appends value to cluster_index attribute

        :param value: Value to append.
        :type value: integer
        """
        self.cluster_index.append(value)

    def set_neighbours(self, value):
        """
        Appends value to neighbours attribute

        :param value: Value to append (list of integers)
        :type value: list
        """
        self.neighbours.append(value)


def DFS(GRAPH, summit_first, cluster_ind):
    """
    Implementation of the depth first search (DFS) algorithm. This is a
    recursive algorithm. See detailed description in ``do_clustering
    function``.

    :param GRAPH: the algorithm updates GRAPH.flag and GRAPH.cluster_index
    :param summit_first: event index by which the research of neighbours begins
    :param cluster_ind: cluster index

    :type GRAPH: Graph object
    :type summit_first: integer
    :type cluster_ind: integer

    """

    GRAPH.flag[summit_first] = 1
    GRAPH.cluster_index[summit_first] = cluster_ind

    for ind_summit_fils in range(len(GRAPH.neighbours[summit_first])):
        summit_fils = GRAPH.neighbours[summit_first][ind_summit_fils]
        if GRAPH.flag[summit_fils] == 0:
            cluster_ind = DFS(GRAPH, summit_fils, cluster_ind)

    return cluster_ind


def waveval(stack_time, t_before, t_after, dt, tdeb):
    """
    Finds the indexes corresponding to t_before and t_after

    :param stack_time: origin time of the event
    :param t_before: time taken before the origin time
    :param t_after: time taken after the origin time
    :param dt: waveform sample distance in seconds
    :param tdeb: start time of the waveform

    :type stack_time: UTCDateTime
    :type t_before: float
    :type t_after: float
    :type dt: float
    :type tdeb: UTCDateTime

    :rtype: integer
    :returns: indexes corresponding to t_before and t_after

    """

    tstart = stack_time-t_before-tdeb
    tend = stack_time+t_after-tdeb
    i_start = int(round(tstart*1./dt))
    i_end = int(round(tend*1./dt))

    return i_start, i_end


def plot_traces(CLUSTER, delay_file, coeff, locs, datadir,
                data_files, threshold):
    """
    Plots the waveforms of all possible event pairs within a cluster.
    On the same figure, displays the superimposed waveforms of the event pair
    for all stations. Also displays the correlation value.

    :param CLUSTER: dictionary containing the event indexes belonging to each
                    cluster
    :param delay_file: file name of the file containing the time delays
    :param coeff: cross-correlation values of all possible event pairs for all
                  stations
    :param locs: list of the whole Waveloc locations (each element of the list
                 is a dictionary)
    :param datadir: data directory path
    :param data_files: list of paths of data files
    :param threshold: correlation coefficient threshold

    :type CLUSTER: dictionary
    :type delay_file: string
    :type coeff: dictionary
    :type locs: list
    :type datadir: string
    :type data_files: list
    :type threshold: float
    """

    # Read the file containing the time delays
    a = BinaryFile(delay_file)
    delay = a.read_binary_file()

    t_before = 0.5
    t_after = 6.0

    tr = {}
    for data_file in data_files:
        wf = Waveform()
        wf.read_from_file(data_file)
        tr[wf.station] = wf.values
        dt = wf.delta
        tdeb = wf.starttime

    list_name = sorted(tr)

    for i in range(1, len(CLUSTER)+1):  # cluster index
        for j in range(len(CLUSTER[i])):  # first event index
            e1 = CLUSTER[i][j]
            for k in range(j+1, len(CLUSTER[i])):  # second event index
                e2 = CLUSTER[i][k]
                co = 0
                fig = plt.figure()
                fig.set_facecolor('white')
                for l in range(len(list_name)):
                    name = list_name[l]
                    if delay[name][e1-1][e2-1] != 'NaN':
                        stack_time_1 = locs[e1-1]['o_time']
                        i_start_1, i_end_1 = waveval(stack_time_1, t_before,
                                                     t_after, dt, tdeb)
                        val1 = tr[name][i_start_1-1:i_end_1]
                        stack_time_2 = locs[e2-1]['o_time']
                        i_start_2, i_end_2 = \
                            waveval(stack_time_2-delay[name][e1-1][e2-1],
                                    t_before, t_after, dt, tdeb)
                        val2 = tr[name][i_start_2-1:i_end_2]
                        t = np.linspace(0, t_after+t_before,
                                        (t_after+t_before)/dt+1)
                        ax = fig.add_subplot(len(list_name), 1, l+1)
                        ax.set_axis_off()
                        ax.plot(t, val1/max(val1), 'k')
                        ax.plot(t, val2/max(val2), 'y--')
                        c = 'k'
                        if coeff[name][e1-1][e2-1] >= threshold:
                            co = co+1
                            c = 'r'
                        ax.text(0.2, 0.5, "%s, %s, %s" %
                                (name, str(coeff[name][e1-1][e2-1]),
                                 str(delay[name][e1-1][e2-1])), color=c)
                fig.suptitle("Cluster : %s ; Event pair : (%s,%s) ; %d" %
                             (str(i), str(e1), str(e2), co))
                plt.show()


def compute_nbsta(event, coeff, threshold):
    """
    Computes the number of stations where the correlation value is greater than
    the threshold for every event pair

    :param event: total number of events located by Waveloc
    :param coeff: cross-correlation coefficients of all possible event pairs
                  for all stations
    :param threshold: correlation value threshold set to form a cluster

    :type event: integer
    :type coeff: dictionary
    :type threshold: float

    :rtype: numpy matrix containing integers
    :returns: matrix of number of stations

    """

    nbsta = []
    for i in xrange(event):
        liste = []
        for k in xrange(i):
            liste.append(0)
        for j in xrange(i, event):
            c = 0
            if i != j:
                for name in sorted(coeff):
                    if coeff[name] and coeff[name][i][j] >= threshold and \
                       coeff[name][i][j] != 'NaN':
                        c = c+1
                liste.append(c)
            else:
                liste.append(0)
        nbsta.append(liste)

    nbsta = np.matrix(nbsta)
    return nbsta


def do_clustering(event, nbsta, nbmin):
    """
    First step: all events are examined one by one. For each of them, the
    indexes of the events where there is a sufficient number of stations (i.e.
    the number exceeds or equals nbmin) are kept in GRAPH.neighbours. All
    events are flagged to 0 and belong to cluster 0 by default. The indexes of
    the events for which neighbours were found are written in the list
    ``summits``.

    From this list (summits), we extract the event which has the greatest
    number of neighbours.  The clustering process will start with that event.

    Then the DFS algorithm is applied: when an event is examined, it is flagged
    to 1 and its cluster number is also updated (GRAPH.cluster_index). The DFS
    algorithm is recursive and explores each possible path until it ends: it
    means that when the research starts with a given event, it will search for
    the the first neighbour, and then the first neighbour of this first
    neighbour and so on... When all the ``first neighbours`` are found, the
    research can concentrates on second neighbours, and so on until all the
    neighbours are found.

    Once the DFS algorithm has found all events linked to summit_first, we
    write the corresponding indexes into the dictionary CLUSTER (where the keys
    are the cluster indexes) and look for events with neighbours which are
    still not flagged to 1. The process keeps going until the whole events with
    neighbours belong to a cluster.

    The function finally returns the dictionary CLUSTER.

    :param event: total number of events in the Waveloc location file
    :param nbsta: 2-D matrix containing the number of stations where the
                  cross-correlation value is greater than a given threshold for
                  all possible event pairs
    :param nbmin: minimum number of stations where the cross-correlation value
                  should be greater than a given threshold to form a cluster

    :type event: integer
    :type nbsta: matrix
    :type nbmin: integer

    :rtype: dictionary
    :returns: indexes of events forming each cluster
    """

    # CODE CHRISTOPHE - CLUSTERING : DEPTH FIRST SEARCH
    neighbours_of_I_summit__horiz = []
    neighbours_of_I_summit__verti = []

    GRAPH = Graph()
    summits = []

    for I in range(event):
        neighbours_of_I_summit__verti = \
            (np.where(nbsta[:, I] >= nbmin)[0]).tolist()[0]
        neighbours_of_I_summit__horiz = \
            (np.where(nbsta[I, :] >= nbmin)[1]).tolist()[0]
        GRAPH.set_neighbours(neighbours_of_I_summit__verti +
                             neighbours_of_I_summit__horiz)
        GRAPH.set_flag(0)
        GRAPH.set_cluster_index(0)
        if neighbours_of_I_summit__verti + neighbours_of_I_summit__horiz:
            summits.append(I)

    if summits:
        nb_max_neighbours = 0
        for ind_summit in range(len(summits)):
            l = len(GRAPH.neighbours[summits[ind_summit]])
            if l > nb_max_neighbours:
                summit_first = summits[ind_summit]
                nb_max_neighbours = l

        ind_summit_first = 0
        cluster_ind = 0
        CLUSTER = {}
        while 1:
            event_index_flagged = []
            event_index_non_flagged_with_neighbours = []
            ind_summit_first = ind_summit_first+1
            cluster_ind = cluster_ind+1
            cluster_ind = DFS(GRAPH, summit_first, cluster_ind)
            for k in range(len(GRAPH.neighbours)):
                if GRAPH.flag[k] == 1 and \
                   GRAPH.cluster_index[k] == ind_summit_first:
                    event_index_flagged.append(k)
                elif GRAPH.flag[k] == 0 and len(GRAPH.neighbours[k]) != 0:
                    event_index_non_flagged_with_neighbours.append(k)

            # add 1 to each event number as the first one is number one (and
            # not zero)
            CLUSTER[cluster_ind] = list(event_index_flagged +
                                        np.ones(len(event_index_flagged),
                                                dtype=np.int))
            if len(event_index_non_flagged_with_neighbours) > 1:
                summit_first = event_index_non_flagged_with_neighbours[0]
            else:
                break

        return CLUSTER

    else:
        return {}


def plot_graphs(locs, stations, nbsta, CLUSTER, nbmin, threshold):
    """
    Displays two figures.  On the first plot, all events are represented. Those
    belonging to a cluster are color-coded and labelled with the cluster index.
    On the second plot, all events are also represented, but those belonging to
    a cluster are colored in black. The links between events of a same cluster
    are also plotted and color-coded in function of the number of stations
    where the correlation value exceeds the threshold.

    Uses mlab from mayavi

    :param locs: list of the whole Waveloc locations (each element of the list
                 is a dictionary)
    :param stations: dictionary of stations
    :param nbsta: 2-D matrix containing the number of stations where the
                  cross-correlation value is greater than a given threshold for
                  all possible event pairs
    :param CLUSTER: dictionary containing the indexes of the events for each
                    cluster
    :param nbmin: minimum number of stations where the cross-correlation value
                  should be greater than a given threshold to form a cluster
    :param threshold: correlation value threshold set to form a cluster

    :type locs: dictionary
    :type stations: dictionary
    :type nbsta: matrix
    :type CLUSTER: dictionary
    :type nbmin: integer
    :type threshold: float
    """
    from mayavi import mlab

    # Event coordinates
    stack_x, stack_y, stack_z = [], [], []
    for loc in locs:
        stack_x.append(loc['x_mean'])
        stack_y.append(loc['y_mean'])
        stack_z.append(-loc['z_mean'])

    # Extract coordinates
    xsta, ysta, zsta = [], [], []
    for sta in sorted(stations):
        xsta.append(stations[sta]['x'])
        ysta.append(stations[sta]['y'])
        zsta.append(stations[sta]['elev'])

    # 3D PLOT USING MAYAVI
    logging.info("Plotting...")
    mlab.figure(1, bgcolor=(1, 1, 1), fgcolor=(0, 0, 0), size=(1000, 900))
    mlab.points3d(xsta, ysta, zsta, color=(1, 0, 0), scale_factor=0.05,
                  mode='cube')
    mlab.axes(extent=[362, 370, 7647, 7653, -0.5, 2.5], color=(0, 0, 0))
    mlab.outline(extent=[362, 370, 7647, 7653, -0.5, 2.5], color=(0, 0, 0))
    mlab.points3d(stack_x, stack_y, stack_z, scale_factor=0.1,
                  color=(0.8, 0.8, 0.8))
    mlab.title("threshold=%s,  nbmin=%s" % (threshold, nbmin), height=0.1,
               size=0.35, color=(0, 0, 0))
    for i_ev in range(len(nbsta)):
        for i_c in range(1, len(CLUSTER)+1):
            if i_ev+1 in CLUSTER[i_c]:
                mlab.points3d(stack_x[i_ev], stack_y[i_ev], stack_z[i_ev],
                              scale_factor=0.1,
                              color=tuple(CZ_Clust_2_color(100 *
                                          (len(CLUSTER)-i_c) / len(CLUSTER))))
                mlab.text3d(stack_x[i_ev], stack_y[i_ev], stack_z[i_ev],
                            str(i_c), color=(0, 0, 0), scale=0.1)
    logging.info("Done!")

    logging.info("Plotting...")
    mlab.figure(2, bgcolor=(1, 1, 1), fgcolor=(0, 0, 0), size=(1000, 900))
    mlab.points3d(xsta, ysta, zsta, color=(1, 0, 0), scale_factor=0.05,
                  mode='cube')
    mlab.axes(extent=[362, 370, 7647, 7653, -0.5, 2.5], color=(0, 0, 0))
    mlab.outline(extent=[362, 370, 7647, 7653, -0.5, 2.5], color=(0, 0, 0))
    mlab.points3d(stack_x, stack_y, stack_z, scale_factor=0.1,
                  color=(0.8, 0.8, 0.8))
    mlab.title("threshold=%s,  nbmin=%s" % (threshold, nbmin), height=0.1,
               size=0.35, color=(0, 0, 0))
    print nbsta
    for ind_I in range(len(nbsta)):
        for ind_J in range(ind_I+1, len(nbsta)):
            W_IJ = nbsta[ind_I, ind_J]
            if W_IJ >= nbmin:
                mlab.points3d(stack_x[ind_J], stack_y[ind_J], stack_z[ind_J],
                              scale_factor=0.1, color=(0, 0, 0))
                mlab.points3d(stack_x[ind_I], stack_y[ind_I], stack_z[ind_I],
                              scale_factor=0.1, color=(0, 0, 0))
                d = (stack_x[ind_J]-stack_x[ind_I],
                     stack_y[ind_J]-stack_y[ind_I],
                     stack_z[ind_J]-stack_z[ind_I])
                norm = np.sqrt(d[0]**2+d[1]**2+d[2]**2)
                mlab.quiver3d(stack_x[ind_I], stack_y[ind_I], stack_z[ind_I],
                              d[0], d[1], d[2],
                              color=tuple(CZ_W_2_color(W_IJ)), mode='2ddash',
                              scale_factor=norm, scale_mode='scalar')
    logging.info("Done!")
    mlab.show()


def do_clustering_setup_and_run(opdict):
    """
    Does clustering by applying the depth first search algorithm and saves the
    result (= a dictionary containing the event indexes forming each cluster)
    in a binary file.  Needs to define the correlation value threshold and the
    minimum number of stations where this threshold should be reached to form a
    cluster (should be done in the options dictionary)

    :param opdict: Dictionary of waveloc options

    """

    base_path = opdict['base_path']
    verbose = opdict['verbose']

    # stations
    stations_filename = os.path.join(base_path, 'lib', opdict['stations'])

    # data
    data_dir = os.path.join(base_path, 'data', opdict['datadir'])
    data_glob = opdict['dataglob']
    data_files = glob.glob(os.path.join(data_dir, data_glob))
    data_files.sort()

    # location file
    locdir = os.path.join(base_path, 'out', opdict['outdir'], 'loc')
    loc_filename = os.path.join(locdir, 'locations.dat')

    # file containing correlation values
    coeff_file = os.path.join(locdir, opdict['xcorr_corr'])
    # Read correlation values
    b = BinaryFile(coeff_file)
    coeff = b.read_binary_file()

    # INPUT PARAMETERS
    nbmin = int(opdict['nbsta'])
    if nbmin > len(coeff.keys()):
        raise Exception('the minimum number of stations cannot be > to the\
                         number of stations !!')
    event = len(coeff.values()[0])
    tplot = float(opdict['clus'])  # threshold for which we save and plot
    cluster_file = "%s/cluster-%s-%s" % (locdir, str(tplot), str(nbmin))

    corr = [opdict['clus']]
    for threshold in corr:
        threshold = float(threshold)
        nbsta = compute_nbsta(event, coeff, threshold)

        CLUSTER = do_clustering(event, nbsta, nbmin)

        if threshold == tplot:

            print "----------------------------------------------"
            print "THRESHOLD : ", threshold, " # STATIONS : ", nbmin
            print "# CLUSTERS : ", len(CLUSTER)
            print CLUSTER

            c = BinaryFile(cluster_file)
            c.write_binary_file(CLUSTER)
            print "Written in %s" % cluster_file

            if verbose:  # PLOT
                # Read location file
                locs = read_locs_from_file(loc_filename)
                # Read station file
                stations = read_stations_file(stations_filename)

                # Look at the waveforms
                # plot_traces(CLUSTER, delay_file, coeff, locs,
                #            data_dir, data_files, threshold)

                # Plot graphs
                plot_graphs(locs, stations, nbsta, CLUSTER, nbmin, threshold)
