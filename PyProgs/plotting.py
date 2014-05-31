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
    nt = plotopt.opdict['nt']
    xc = plotopt.opdict['x_loc']
    yc = plotopt.opdict['y_loc']
    zc = plotopt.opdict['z_loc']
    tc = plotopt.opdict['t_loc_rel']
    ic = ugrid_closest_point_index(x, y, z, xc, yc, zc)
    it = int(tc/dt)
    t = np.arange(0, dt*nt, dt)

    # open grid_file
    f = h5py.File(grid_filename, 'r')
    grid4D = f['migrated_grid']
    grid1D = grid4D[:,it]
    f.close()

    # get max, min values for norm
    val_min = _round_sig(np.min(grid1D))
    val_max = _round_sig(np.max(grid1D))
    val_mean = _round_sig((val_min+val_max)/2.)
    norm = mpl.colors.Normalize(vmin=val_min,
                                vmax=val_max)

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

    # get max_val etc
    max_val = np.empty(nt, dtype='float')
    max_x = np.empty(nt, dtype='float')
    max_y = np.empty(nt, dtype='float')
    max_z = np.empty(nt, dtype='float')
    f = h5py.File(stack_filename, 'r')
    max_val = f['max_val_smooth'][:]
    max_x = f['max_x'][:]
    max_y = f['max_y'][:]
    max_z = f['max_z'][:]

    f.close()
    ###########
    # do plot #
    ###########

    import pdb; pdb.set_trace()

    cmap = mpl.cm.hot_r

    plt.clf()
    fig = plt.figure()

    fig.suptitle('x = %.2fkm  y = %.2fkm  z = %.2fkm' % (xc, yc, -zc))

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
               extent=[0, np.max(x)-np.min(x), -np.max(z), -np.min(z)],
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
               extent=[0, np.max(y)-np.min(y), -np.max(z), -np.min(z)],
               cmap=cmap, norm=norm)
    p.xaxis.set_ticks_position('bottom')
    p.tick_params(labelsize=10)
    plt.xlabel('y (km wrt ref)', size=10)
    plt.ylabel('z (km up)', size=10)

    # choose portion of time series to plot
    otime_window = plotopt.opdict['plot_otime_window']
    illim = max(it-int(otime_window/dt), 0)
    irlim = min(it+int(otime_window/dt), nt-1)
    llim = illim*dt
    rlim = irlim*dt

    # plot max value
    p = plt.subplot(422, frameon=False)
    pos = list(p.get_position().bounds)
    fig.text(pos[0]-.05, pos[1]+pos[3], '(b)', fontsize=12)
    p.tick_params(labelsize=10)
    plt.plot(t, max_val, 'k')
    p.xaxis.set_ticks_position('bottom')
    p.xaxis.set_ticklabels('')
    plt.ylabel('Stack max', size=10)
    p.yaxis.set_ticks_position('right')
    p.set_xlim(llim, rlim)
    p.set_ylim(np.min(max_val[illim:irlim]), max(max_val))
    plt.vlines(tc, np.min(max_val[illim:irlim]), max(max_val), 'r',
               linewidth=2)
    if 't_err' in plotopt.opdict:
        t_left, t_right = plotopt.opdict['t_err']
        plt.axvspan(tc-t_left, tc+t_right, facecolor='r', alpha=0.2)

    # plot max x
    p = plt.subplot(424, frameon=False)
    pos = list(p.get_position().bounds)
    fig.text(pos[0]-.05, pos[1]+pos[3], '(c)', fontsize=12)
    p.tick_params(labelsize=10)
    plt.scatter(t[illim:irlim], max_x[illim:irlim], s=40,
                c=max_val[illim:irlim], marker='.', linewidths=(0, ),
                clip_on=False, cmap=cmap, norm=norm)
    p.xaxis.set_ticks_position('bottom')
    p.xaxis.set_ticklabels('')
    plt.ylabel('x (km)', size=10)
    p.yaxis.set_ticks_position('right')
    p.set_xlim(llim, rlim)
    plt.hlines(xc, llim, rlim, 'r', linewidth=2)
    plt.vlines(tc, min(max_x), max(max_x), 'r', linewidth=2)
    if 'x_err' in plotopt.opdict:
        x_low, x_high = plotopt.opdict['x_err']
        plt.axhspan(xc-x_low, xc+x_high, facecolor='r', alpha=0.2)
    if 't_err' in plotopt.opdict:
        t_left, t_right = plotopt.opdict['t_err']
        plt.axvspan(tc-t_left, tc+t_right, facecolor='r', alpha=0.2)

    # plot max y
    p = plt.subplot(426, frameon=False)
    pos = list(p.get_position().bounds)
    fig.text(pos[0]-.05, pos[1]+pos[3], '(e)', fontsize=12)
    p.tick_params(labelsize=10)
    plt.scatter(t[illim:irlim], max_y[illim:irlim], s=40,
                c=max_val[illim:irlim], marker='.', linewidths=(0, ),
                clip_on=False, cmap=cmap, norm=norm)
    p.xaxis.set_ticks_position('bottom')
    p.xaxis.set_ticklabels('')
    plt.ylabel('y (km)', size=10)
    p.yaxis.set_ticks_position('right')
    p.set_xlim(llim, rlim)
    plt.hlines(yc, llim, rlim, 'r', linewidth=2)
    plt.vlines(tc, min(max_y), max(max_y), 'r', linewidth=2)
    if 'y_err' in plotopt.opdict:
        y_low, y_high = plotopt.opdict['y_err']
        plt.axhspan(yc-y_low, yc+y_high, facecolor='r', alpha=0.2)
    if 't_err' in plotopt.opdict:
        t_left, t_right = plotopt.opdict['t_err']
        plt.axvspan(tc-t_left, tc+t_right, facecolor='r', alpha=0.2)

    # plot max z
    # note : for calculation purposes z is posivive down
    # for plotting, z should be positive up
    p = plt.subplot(428, frameon=False)
    pos = list(p.get_position().bounds)
    fig.text(pos[0]-.05, pos[1]+pos[3], '(g)', fontsize=12)
    p.tick_params(labelsize=10)
    plt.scatter(t[illim:irlim], -max_z[illim:irlim], s=40,
                c=max_val[illim:irlim], marker='.', linewidths=(0, ),
                clip_on=False, cmap=cmap, norm=norm)
    plt.xlabel('Time (s)', size=10)
    p.xaxis.set_ticks_position('bottom')
    plt.ylabel('z (km up)', size=10)
    p.yaxis.set_ticks_position('right')
    p.set_xlim(llim, rlim)
    plt.hlines(-zc, llim, rlim, 'r', linewidth=2)
    plt.vlines(tc, -max(max_z), -min(max_z), 'r', linewidth=2)
    if 'z_err' in plotopt.opdict:
        z_low, z_high = plotopt.opdict['z_err']
        plt.axhspan(-zc-z_low, -zc+z_high, facecolor='r', alpha=0.2)
    if 't_err' in plotopt.opdict:
        t_left, t_right = plotopt.opdict['t_err']
        plt.axvspan(tc-t_left, tc+t_right, facecolor='r', alpha=0.2)

    # add independent colorbar
    ax1 = fig.add_axes([0.40, 0.03, 0.2, 0.015])
    ax1.tick_params(labelsize=8)
    ax1.xaxis.set_ticks_position('bottom')
    mpl.colorbar.ColorbarBase(ax1, cmap=cmap, norm=norm,
                              orientation='horizontal',
                              ticks=[val_min, val_mean, val_max])
    pos = list(ax1.get_position().bounds)
    fig.text(pos[0]+pos[2]/2., pos[1]+pos[3]+0.01, 'Stack max', fontsize=8,
             horizontalalignment='center', verticalalignment='bottom')

    plt.savefig(fig_filename)
    print(fig_filename)


def _round_sig(x, sig=2):
    from numpy import int, log10, floor, abs
    if abs(x > 0.):
        result = round(x, sig-int(floor(log10(abs(x))))-1)
    else:
        result = 0.
    return result
