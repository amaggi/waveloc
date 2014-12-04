import os
import h5py
import glob
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from obspy.core import UTCDateTime
from ugrids import ugrid_closest_point_index, select_points_closeto_plane
from hdf5_grids import get_interpolated_time_ugrids
from locations_trigger import read_locs_from_file
from OP_waveforms import read_data_compatible_with_time_dict
from migration import do_migration_loop_continuous
from plot_options import PlotOptions


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
    grid1D = grid4D[:, it]
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
        raise NotImplemented('Plotting of user grids not implemented yet')

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
        raise NotImplemented('Plotting of user grids not implemented yet')
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
    otime_window = plotopt.opdict['otime_window']
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


def plotLocationWaveforms(plotopt):
    """
    Plots the waveforms re-aligned after a waveloc location.
    Uses information in the plotopt object
    """

    # get filenames
    fig_filename = plotopt.getWfmFigFilename()

    # get location parameters
    dt = plotopt.opdict['dt']
    loc = plotopt.opdict['loc']
    otime = loc['o_time']
    otime_left = -loc['o_err_left']
    otime_right = loc['o_err_right']
    stack_wfm = plotopt.opdict['stack_wfm']
    start_time = plotopt.opdict['start_time']
    data_dict = plotopt.opdict['data_dict']
    mig_dict = plotopt.opdict['mig_dict']

    # get the time range
    t = np.arange(len(stack_wfm))*dt - (otime - start_time)
    #t = np.arange(len(stack_wfm))*dt + start_time

    # get and sort the stations
    stations = data_dict.keys()
    stations.sort()
    n_traces = len(stations)+1

    # start the plot (set up the two panels)
    plt.clf()
    fig = plt.figure()
    ax = fig.add_subplot(n_traces, 2, 1, title='Data')
    ax.text(-0.25, 2.0, "(a)", transform=ax.transAxes)
    ax.set_axis_off()
    ax = fig.add_subplot(n_traces, 2, 2, title='Characteristic function')
    ax.text(-0.12, 2.0, "(b)", transform=ax.transAxes)
    ax.set_axis_off()

    # iterate over stations
    i = 0
    for sta in stations:
        # plot the data in the first column
        if i != len(stations)-1:
            ax = fig.add_subplot(n_traces, 2, 2*i+1)
            ax.set_axis_off()
        else:
            ax = fig.add_subplot(n_traces, 2, 2*i+1, xlabel='time (s)')
            ax.xaxis.set_ticks_position('none')
            ax.yaxis.set_ticks([])
            ax.spines['top'].set_color('none')
            ax.spines['left'].set_color('none')
            ax.spines['right'].set_color('none')
        ax.plot(t, data_dict[sta], 'b')
        ax.axvspan(otime_left, otime_right, facecolor='r', alpha=0.2)
        # add the station name
        pos = list(ax.get_position().bounds)
        fig.text(pos[0]-0.01, pos[1]+pos[3]/2., sta, fontsize=10,
                 horizontalalignment='right',
                 verticalalignment='center')

        # plot the characteristic function in the second column
        if i != len(stations)-1:
            ax = fig.add_subplot(n_traces, 2, 2*i+2)
            ax.set_axis_off()
        else:
            ax = fig.add_subplot(n_traces, 2, 2*i+2, xlabel='time (s)')
            ax.xaxis.set_ticks_position('none')
            ax.yaxis.set_ticks([])
            ax.spines['top'].set_color('none')
            ax.spines['left'].set_color('none')
            ax.spines['right'].set_color('none')
        ax.plot(t, mig_dict[sta], 'b')
        ax.axvspan(otime_left, otime_right, facecolor='r',  alpha=0.2)
        # add the maximum kurtosis value
        pos = list(ax.get_position().bounds)
        fig.text(pos[0]+pos[2]+0.05, pos[1], '%.1f' % np.max(mig_dict[sta]),
                 fontsize=10, horizontalalignment='right')
        i = i+1

    fig.suptitle(otime.isoformat(), x=0.5, y=0.05)

    plt.savefig(fig_filename)
    print(fig_filename)
    plt.clf()


def do_plotting_setup_and_run(opdict, plot_wfm=True, plot_grid=True):
    """
    Plot the results of a wavloc run (migration and location). All options and
    parameters are taken from an opdict.

    :param opdict: WavlocOptions.opdict that contains the options / parameters.
    :param plot_wfm: If ``True`` plots waveforms after location (filtered data
        and kurtosis).
    :param plot_grid: If ``True``plots the migration grid.

    :type plot_wfm: boolean
    :type plot_grid: boolean
    """

    # get / set info
    base_path = opdict['base_path']

    locfile = os.path.join(base_path, 'out', opdict['outdir'], 'loc',
                           'locations.dat')
    stackfile = os.path.join(base_path, 'out', opdict['outdir'], 'stack',
                             'combined_stack_all.hdf5')

    data_dir = os.path.join(base_path, 'data', opdict['datadir'])

    data_glob = opdict['dataglob']
    data_files = glob.glob(os.path.join(data_dir, data_glob))
    data_files.sort()

    kurt_glob = opdict['kurtglob']
    kurt_files = glob.glob(os.path.join(data_dir, kurt_glob))
    kurt_files.sort()
    mig_files = kurt_files

    if opdict['kderiv']:
        grad_glob = opdict['gradglob']
        grad_files = glob.glob(os.path.join(data_dir, grad_glob))
        grad_files.sort()
        mig_files = grad_files

        if opdict['gauss']:
            gauss_glob = opdict['gaussglob']
            gauss_files = glob.glob(os.path.join(data_dir, gauss_glob))
            gauss_files.sort()
            mig_files = gauss_files

    # figdir = os.path.join(base_path, 'out', opdict['outdir'], 'fig')

    # read time grid information
    x, y, z, time_grids = get_interpolated_time_ugrids(opdict)
    grid_info = (x, y, z)

    # read locations
    locs = read_locs_from_file(locfile)

    # open stack file
    f_stack = h5py.File(stackfile, 'r')
    max_val = f_stack['max_val_smooth']
    stack_start_time = UTCDateTime(max_val.attrs['start_time'])

    for loc in locs:
        # generate the grids
        o_time = loc['o_time']
        start_time = o_time-opdict['plot_tbefore']
        end_time = o_time+opdict['plot_tafter']

        # get location coordinates
        xm = loc['x_mean']
        ym = loc['y_mean']
        zm = loc['z_mean']

        # get the corresponding travel-times for time-shifting
        ttimes = {}
        for sta in time_grids.keys():
            ic = ugrid_closest_point_index(x, y, z, xm, ym, zm)
            ttimes[sta] = time_grids[sta][ic]

        tshift_migration = max(ttimes.values())
        start_time_migration = start_time-tshift_migration
        end_time_migration = end_time+tshift_migration

        if plot_grid:

            # read data
            mig_dict, delta = \
                read_data_compatible_with_time_dict(mig_files, time_grids,
                                                    start_time_migration,
                                                    end_time_migration)
            # do migration
            grid_filename, n_buf, nt, stack_shift_time = \
                do_migration_loop_continuous(opdict, mig_dict, delta,
                                             start_time_migration, grid_info,
                                             time_grids, keep_grid=True)
            # plot
            plotopt = PlotOptions(opdict)
            grid_filename = os.path.basename(grid_filename)
            stack_filename = 'stack_all_' + grid_filename.split('_')[-1]
            plotopt.opdict['grid_filename'] = grid_filename
            plotopt.opdict['stack_filename'] = stack_filename
            plotopt.opdict['loc'] = loc

            plotopt.opdict['dt'] = delta
            plotopt.opdict['n_buf'] = n_buf
            plotopt.opdict['nt'] = nt
            plotopt.opdict['t_loc_rel'] = \
                o_time - start_time_migration + stack_shift_time
            plotopt.opdict['x_loc'] = xm
            plotopt.opdict['y_loc'] = ym
            plotopt.opdict['z_loc'] = zm
            plotopt.opdict['start_time'] = \
                start_time_migration - stack_shift_time
            plotWavelocResults(plotopt)

        if plot_wfm:

            # read data
            data_dict, delta = \
                read_data_compatible_with_time_dict(data_files, time_grids,
                                                    start_time_migration,
                                                    end_time_migration)
            mig_dict, delta = \
                read_data_compatible_with_time_dict(mig_files, time_grids,
                                                    start_time_migration,
                                                    end_time_migration)
            # cut desired portion out of data
            for sta in data_dict.keys():
                tmp = data_dict[sta]

                # alignment on origin time
                istart = np.int(np.round((start_time + ttimes[sta] -
                                          start_time_migration) / delta))
                iend = istart + np.int(np.round((opdict['plot_tbefore'] +
                                                 opdict['plot_tafter']) /
                                                delta))

                # sanity check in case event is close to start or end of data
                if istart < 0:
                    istart = 0
                if iend > len(tmp):
                    iend = len(tmp)
                data_dict[sta] = tmp[istart:iend]
                # do slice
                tmp = mig_dict[sta]
                mig_dict[sta] = tmp[istart:iend]

            # retrieve relevant portion of stack max
            istart = np.int(np.round((o_time - opdict['plot_tbefore'] -
                                      stack_start_time) / delta))
            iend = istart + np.int(np.round((opdict['plot_tbefore'] +
                                             opdict['plot_tafter']) / delta))
            # sanity check in case event is close to start or end of data
            if istart < 0:
                start_time = start_time + np.abs(istart)*delta
                istart = 0
            if iend > len(max_val):
                iend = len(max_val)
            # do slice
            stack_wfm = max_val[istart:iend]

            # plot
            plotopt = PlotOptions(opdict)
            plotopt.opdict['dt'] = delta
            plotopt.opdict['data_dict'] = data_dict
            plotopt.opdict['mig_dict'] = mig_dict
            plotopt.opdict['stack_wfm'] = stack_wfm
            plotopt.opdict['loc'] = loc
            plotopt.opdict['start_time'] = start_time
            plotLocationWaveforms(plotopt)

    f_stack.close()
