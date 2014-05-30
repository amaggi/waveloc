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
    it = plotopt.opdict['dt']
    xc = plotopt.opdict['x_loc']
    yc = plotopt.opdict['y_loc']
    zc = plotopt.opdict['z_loc']
    ic = ugrid_closest_point_index(x, y, z, xc, yc, zc)

    # open grid_file
    f = h5py.File(grid_filename, 'r')
    grid4D = f['migrated_grid']
    grid3D = grid4D[:,it]
    f.close()

    # get max, min values for norm
    norm = mpl.colors.Normalize(vmin=np.min(grid3D), vmax=np.max(grid3D))

    # get xy, xz, yz cuts
    if ugrid_type == 'FULL':
        x_vector = np.unique(x)
        y_vector = np.unique(y)
        z_vector = np.unique(z)
        x_vector.sort()
        y_vector.sort()
        z_vector.sort()
        dx = x_vector[1]-x_vector[0]
        dy = y_vector[1]-y_vector[0]
        dz = z_vector[1]-z_vector[0]
        nx = len(x_vector)
        ny = len(y_vector)
        nz = len(z_vector)

        indexes = select_points_closeto_plane(x, y, z, 'z', zc, dz)
        xy_cut = grid3D[indexes].reshape(nx, ny)

        indexes = select_points_closeto_plane(x, y, z, 'y', yc, dy)
        xz_cut = grid3D[indexes].reshape(nx, nz)

        indexes = select_points_closeto_plane(x, y, z, 'x', xc, dx)
        yz_cut = grid3D[indexes].reshape(ny, nz)
    else:
        raise NotImplemented ('Plotting of user grids not implemented yet')

    ###########
    # do plot #
    ###########

    cmap = mpl.cm.hot_r

    plt.clf()
    fit = plt.figure()

    fig.suptitle('x = %.2fkm  y = %.2fkm  z = %.2fkm' % (xc, yc, zc))

    # plot xy plane
    p = plt.subplot(221)
    pos = list(p.get_position().bounds)
    fig.text(pos[0]-0.08, pos[1]+pos[3], '(a)', fontsize=12)
    if ugrid_type == 'FULL':
        plt.imshow(xy_cut.T, origin='lower', interpolation='none',
                   extent=[np.min(x), np.max(x), np.min(y), np.max(y)],
                   cmap=cmap, norm=norm)
    else:
        raise NotImplemented ('Plotting of user grids not implemented yet')
    p.tick_params(labelsize=10)
    p.xaxis.set_ticks_position('bottom')
    plt.xlabel('x (km wrt ref)', size=10)
    plt.ylabel('y (km wrt ref)', size=10)

    # add independent colorbar
    ax1 = fig.add_axes([0.40, 0.03, 0.2, 0.015])
    ax1.tick_params(labelsize=8)
    ax1.xaxis.set_ticks_position('bottom')
    mpl.colorbar.ColorbarBase(ax1, cmap=cmap, norm=norm,
                              orientation='horizontal',
                              ticks=[0, int(np.max(max_val)/2),
                                     int(np.max(max_val))])
    pos = list(ax1.get_position().bounds)
    fig.text(pos[0]+pos[2]/2., pos[1]+pos[3]+0.01, 'Stack max', fontsize=8,
             horizontalalignment='center', verticalalignment='bottom')

    plt.savefig(fig_filename)
