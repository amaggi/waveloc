import os
import h5py
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from ugrids import ugrid_closest_point_index, select_points_closeto_plane


def plotWavelocResults(plotopt):
    """
    Creates plot of waveloc results
    """

    # get filenames
    grid_filename = plotopt.getGridFilename()
    stack_filename = plotopt.getStackFilename()
    fig_filename = plotopt.getFigFilename()

    # get ugrid type
    ugrid_type = plotopt.opdict['ugrid_type']

    # get x, y, z
    x, y, z = plotopt.getXYZ()

    # get indexes of location
    dt = plotopt.opdict['dt']
    xc = plotopt.opdict['x_loc']
    yc = plotopt.opdict['y_loc']
    zc = plotopt.opdict['z_loc']
    tc = plotopt.opdict['t_loc_rel']
    ic = ugrid_closest_point_index(x, y, z, xc, yc, zc)
    it = int(tc/dt)

    # open grid_file
    f = h5py.File(grid_filename, 'r')
    grid4D = f['migrated_grid']
    grid1D = grid4D[:,it]
    f.close()

    # get max, min values for norm
    min_val = _round_sig(np.min(grid1D))
    max_val = _round_sig(np.max(grid1D))
    mean_val = _round_sig((max_val+min_val)/2.)
    norm = mpl.colors.Normalize(vmin=min_val,
                                vmax=max_val)

    # get xy, xz, yz cuts
    if ugrid_type == 'FULL':
        # is a regular grid - can treat it as such.
        nx = len(np.unique(x))
        ny = len(np.unique(y))
        nz = len(np.unique(z))
        grid3D = grid1D.reshape(nx, ny, nz)

        # cut the grid
        ix, iy, iz = np.unravel_index(ic, (nx, ny, nz))
        xy_cut = grid3D[:, :, iz]
        xz_cut = grid3D[:, iy, :]
        yz_cut = grid3D[ix, :, :]

    else:
        raise NotImplemented ('Plotting of user grids not implemented yet')

    ###########
    # do plot #
    ###########

    cmap = mpl.cm.hot_r

    plt.clf()
    fig = plt.figure()

    fig.suptitle('x = %.2fkm  y = %.2fkm  z = %.2fkm' % (xc, yc, zc))

    # plot xy plane
    p = plt.subplot(221)
    pos = list(p.get_position().bounds)
    fig.text(pos[0]-0.08, pos[1]+pos[3], '(a)', fontsize=12)
    if ugrid_type == 'FULL':
        plt.imshow(xy_cut.T, origin='lower', interpolation='none',
                   extent=[0, np.max(x)-np.min(x), 0, np.max(y)-np.min(y)],
                   cmap=cmap, norm=norm)
    else:
        raise NotImplemented ('Plotting of user grids not implemented yet')
    p.tick_params(labelsize=10)
    p.xaxis.set_ticks_position('bottom')
    plt.xlabel('x (km wrt ref)', size=10)
    plt.ylabel('y (km wrt ref)', size=10)

    # plot xz plane
    p = plt.subplot(425)
    pos = list(p.get_position().bounds)
    fig.text(pos[0]-0.08, pos[1]+pos[3], '(d)', fontsize=12)
    plt.imshow(xz_cut.T, origin='upper', interpolation='none',
               extent=[0, np.max(x)-np.min(x), 0, np.max(z)-np.min(z)],
               cmap=cmap, norm=norm)
    p.tick_params(labelsize=10)
    p.xaxis.set_ticks_position('top')
    p.xaxis.set_ticklabels('')
    plt.ylabel('z (km up)', size=10)

    # plot yz plane
    p = plt.subplot(427)
    pos = list(p.get_position().bounds)
    fig.text(pos[0]-0.08, pos[1]+pos[3], '(f)', fontsize=12)
    plt.imshow(yz_cut.T, origin='upper', interpolation='none',
               extent=[0, np.max(y)-np.min(y), 0, np.max(z)-np.min(z)],
               cmap=cmap, norm=norm)
    p.xaxis.set_ticks_position('bottom')
    p.tick_params(labelsize=10)
    plt.xlabel('y (km wrt ref)', size=10)
    plt.ylabel('z (km up)', size=10)

    # add independent colorbar
    ax1 = fig.add_axes([0.40, 0.03, 0.2, 0.015])
    ax1.tick_params(labelsize=8)
    ax1.xaxis.set_ticks_position('bottom')
    mpl.colorbar.ColorbarBase(ax1, cmap=cmap, norm=norm,
                              orientation='horizontal',
                              ticks=[min_val, mean_val, max_val])
    pos = list(ax1.get_position().bounds)
    fig.text(pos[0]+pos[2]/2., pos[1]+pos[3]+0.01, 'Stack max', fontsize=8,
             horizontalalignment='center', verticalalignment='bottom')

    plt.savefig(fig_filename)


def _round_sig(x, sig=2):
    from numpy import int, log10, floor, abs
    if abs(x > 0.):
        result = round(x, sig-int(floor(log10(abs(x))))-1)
    else:
        result = 0.
    return result
