import os
import h5py

from locations_trigger import read_locs_from_file
from locations_prob import read_prob_locs_from_file
from plot_mpl import plotProbLoc


def do_probloc_plotting_setup_and_run(opdict):
    """
    Plot the results of a wavloc run (migration and location using probability
    density). All options and parameters are taken from an opdict.

    :param opdict: WavlocOptions.opdict that contains the options / parameters.
    """

    # get / set info
    base_path = opdict['base_path']
    space_only = opdict['probloc_spaceonly']

    locfile = os.path.join(base_path, 'out', opdict['outdir'], 'loc',
                           'locations.dat')
    problocfile = os.path.join(base_path, 'out', opdict['outdir'], 'loc',
                               'locations_prob.dat')
    problocgrid = os.path.join(base_path, 'out', opdict['outdir'], 'loc',
                               'locations_prob.hdf5')

    figdir = os.path.join(base_path, 'out', opdict['outdir'], 'fig')

    # read locations
    locs = read_locs_from_file(locfile)
    prob_locs = read_prob_locs_from_file(problocfile)

    # open hdf5 file
    f = h5py.File(problocgrid, 'r')

    # for each loc
    for i in xrange(len(locs)):
        loc = locs[i]
        prob_loc = prob_locs[i]

        # hdf5 group name is prob loc origin time as string
        grp = f[prob_loc['o_time'].isoformat()]

        # do plotting
        plotProbLoc(grp, prob_loc, loc, figdir,  space_only)

    f.close()
