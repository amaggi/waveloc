import os
import h5py
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from obspy.core import utcdatetime


def plotLocationWaveforms(loc, start_time, dt, data_dict, grad_dict, stack_wfm,
                          fig_dir):
    """
    Creates plot for located waveforms.  Assumes data and grad are ready
    for plotting. TODO : Flesh out this doc-string.

    :param loc:
    :param start_time:
    :param dt:
    :param data_dict:
    :param grad_dict:
    :param stack_wfm:
    :param fig_dir:
    """

    otime = loc['o_time']
    otime_left = -loc['o_err_left']
    otime_right = loc['o_err_right']

    plot_filename = os.path.join(fig_dir, 'loc_%s.pdf' % (otime.isoformat()))

    t = np.arange(len(stack_wfm))*dt - (otime - start_time)

    stations = data_dict.keys()
    stations.sort()
    n_traces = len(stations)+1

    plt.clf()
    fig = plt.figure()
    ax = fig.add_subplot(n_traces, 2, 1, title='Data')
    ax.text(-0.25, 2.0, "(a)", transform=ax.transAxes)
    ax.set_axis_off()
    ax = fig.add_subplot(n_traces, 2, 2, title='Characteristic function')
    ax.text(-0.12, 2.0, "(b)", transform=ax.transAxes)
    ax.set_axis_off()

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
        # plot the kurtosis gradient in the second column
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
        ax.plot(t, grad_dict[sta], 'b')
        ax.axvspan(otime_left, otime_right, facecolor='r',  alpha=0.2)
        # add the maximum kurtosis value
        pos = list(ax.get_position().bounds)
        fig.text(pos[0]+pos[2]+0.05, pos[1], '%.1f' % np.max(grad_dict[sta]),
                 fontsize=10, horizontalalignment='right')
        i = i+1

    fig.suptitle(otime.isoformat(), x=0.5, y=0.05)

    plt.savefig(plot_filename)
    plt.clf()


def plotLocationGrid(loc, grid_filename, fig_dir, otime_window):
    """
    Plots location grid. TODO : Flesh out this doc-string.

    :param loc:
    :param grid_filename:
    :param fig_dir:
    :param otime_window:

    """

    # get useful information from grid file
    f = h5py.File(grid_filename, 'r')
    sg = f['stack_grid']
    n_buf, nt = sg.shape
    stack_starttime = utcdatetime.UTCDateTime(sg.attrs['start_time'])
    dt = sg.attrs['dt']
    f.close()

    # get location info
    o_time = loc['o_time']
    x_mean = loc['x_mean']
    y_mean = loc['y_mean']
    z_mean = loc['z_mean']
    o_err_left = loc['o_err_left']
    o_err_right = loc['o_err_right']
    x_sigma = loc['x_sigma']
    y_sigma = loc['y_sigma']
    z_sigma = loc['z_sigma']

    #get indexes correponding to location
    it_true = np.int(np.round((o_time-stack_starttime)/dt))
    if it_true > nt-1:
        msg = 'Origin time after last time for plot by %3f seconds. \n\
            Increase plot_tafter.' % ((it_true - nt + 1)*dt)
        raise UserWarning(msg)
    if it_true < 0:
        msg = 'Origin time before first time for plot by %3f seconds. \n\
            Increase plot_tbefore.' % (-1*it_true*dt)
        raise UserWarning(msg)

    # get indexes corresponding to location uncertainties
    # times are wrt stack_starttime
    t_left = o_time - o_err_left - stack_starttime
    t_right = o_time + o_err_right - stack_starttime
    # coordinates are absolute
    x_low = x_mean - x_sigma
    y_low = y_mean - y_sigma
    z_low = z_mean - z_sigma
    x_high = x_mean + x_sigma
    y_high = y_mean + y_sigma
    z_high = z_mean + z_sigma

    plot_info = {}
    plot_info['dt'] = dt
    plot_info['start_time'] = stack_starttime
    plot_info['o_time'] = o_time
    plot_info['t_err'] = (t_left, t_right)
    plot_info['x_err'] = (x_low, x_high)
    plot_info['y_err'] = (y_low, y_high)
    plot_info['z_err'] = (z_low, z_high)
    plot_info['true_values'] = (x_mean, y_mean, z_mean, o_time-stack_starttime)
    plot_info['grid_file'] = grid_filename

    plotDiracTest(plot_info, fig_dir, otime_window)


def plotDiracTest(test_info, fig_dir, otime_window):
    """
    Creates plot for synthetic test.
    TODO : flesh out this doc-string

    :param test_info:
    :param fig_dir:
    :param otime_window:
    """

    # set up plot using info from test_info
    nx, ny, nz, nt = test_info['grid_shape']
    dx, dy, dz, dt = test_info['grid_spacing']
    x_orig, y_orig, z_orig = test_info['grid_orig']
    ix_true, iy_true, iz_true, it_true = test_info['true_indexes']
    if 'true_values' in test_info:
        x_true, y_true, z_true, t_true = test_info['true_values']
    stack_start_time = test_info['start_time']
    grid_filename = test_info['dat_file']
    stack_filename = test_info['stack_file']
    if 'o_time' in test_info:
        fig_filename = os.path.join(fig_dir, "grid_%s.pdf" %
                                    test_info['o_time'].isoformat())
    else:
        fig_filename = \
            os.path.join(fig_dir, "%s.pdf" %
                         os.path.splitext(os.path.basename(grid_filename))[0])

    # read the stack file
    f = h5py.File(grid_filename, 'r')
    stack_grid = f['stack_grid']

    stack_3D = stack_grid[:, it_true].reshape(nx, ny, nz)

    col = plt.cm.hot_r

    # cut through the true location at the true time
    xy_cut = stack_3D[:, :, iz_true]
    xz_cut = stack_3D[:, iy_true, :]
    yz_cut = stack_3D[ix_true, :, :]

    # extract the max stacks
    f_stack = h5py.File(stack_filename, 'r')
    # if have a smoothed version, use it for the plots
    if 'max_val_smooth' in f_stack:
        max_val = f_stack['max_val_smooth']
    else:
        max_val = f_stack['max_val']
    max_x = f_stack['max_x']
    max_y = f_stack['max_y']
    max_z = f_stack['max_z']

    # set up the 4 axes
    x = np.arange(nx)*dx
    y = np.arange(ny)*dy
    z = (np.arange(nz)*dz+z_orig)*(-1)
    # setup of t-axis depends on type of stack_start_time
    if type(stack_start_time) == float:
        t = np.arange(nt)*dt+stack_start_time
    else:
        t = np.arange(nt)*dt

    # do plot
    plt.clf()
    fig = plt.figure()

    if 'true_values' in test_info and 'o_time' in test_info:
        fig.suptitle('%s   x = %.2fkm  y = %.2fkm  z = %.2fkm' %
                     (test_info['o_time'].isoformat(), x_true, y_true, z_true))

    # plot xy plane
    if dx > 0 and dy > 0:
        p = plt.subplot(2, 2, 1)
        pos = list(p.get_position().bounds)
        fig.text(pos[0]-0.08, pos[1]+pos[3], '(a)', fontsize=12)
        plt.imshow(xy_cut.T, origin='lower', interpolation='none',
                   extent=[np.min(x), np.max(x), np.min(y), np.max(y)],
                   cmap=col)
        p.tick_params(labelsize=10)
        p.xaxis.set_ticks_position('bottom')
        plt.xlabel('x (km wrt ref)', size=10)
        plt.ylabel('y (km wrt ref)', size=10)

    #plot xz plane
    if dx > 0 and dz > 0:
        p = plt.subplot(4, 2, 5)
        pos = list(p.get_position().bounds)
        fig.text(pos[0]-0.08, pos[1]+pos[3], '(d)', fontsize=12)
        plt.imshow(xz_cut.T, origin='upper', interpolation='none',
                   extent=[np.min(x), np.max(x), np.min(z), np.max(z)],
                   cmap=col)
        p.tick_params(labelsize=10)
        p.xaxis.set_ticks_position('top')
        p.xaxis.set_ticklabels('')
        plt.ylabel('z (km up)', size=10)

    # plot yz plane
    if dy > 0 and dz > 0:
        p = plt.subplot(4, 2, 7)
        pos = list(p.get_position().bounds)
        fig.text(pos[0]-0.08, pos[1]+pos[3], '(f)', fontsize=12)
        plt.imshow(yz_cut.T, origin='upper', interpolation='none',
                   extent=[np.min(y), np.max(y), np.min(z), np.max(z)],
                   cmap=col)
        p.xaxis.set_ticks_position('bottom')
        p.tick_params(labelsize=10)
        plt.xlabel('y (km wrt ref)', size=10)
        plt.ylabel('z (km up)', size=10)

    # choose portion of time series to plot
    if 'true_values' in test_info:
        llim = max(t_true-otime_window, t[0])
        rlim = min(t_true+otime_window, t[-1])
    else:
        llim = max(t[it_true]-otime_window, t[0])
        rlim = min(t[it_true]+otime_window, t[-1])

    illim = int((llim-t[0])/dt)
    irlim = int((rlim-t[0])/dt)

    # plot max value
    p = plt.subplot(4, 2, 2,  frameon=False)
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
    plt.vlines(t[it_true], np.min(max_val[illim:irlim]), max(max_val), 'r',
               linewidth=2)
    if 't_err' in test_info:
        t_left, t_right = test_info['t_err']
        plt.axvspan(t_left, t_right, facecolor='r', alpha=0.2)

    # put the origin back in for the last plots
    x = np.arange(nx)*dx+x_orig
    y = np.arange(ny)*dy+y_orig
    z = np.arange(nz)*dz+z_orig

    # plot max x
    if dx > 0:
        p = plt.subplot(4, 2, 4,  frameon=False)
        pos = list(p.get_position().bounds)
        fig.text(pos[0]-.05, pos[1]+pos[3], '(c)', fontsize=12)
        p.tick_params(labelsize=10)
        plt.scatter(t[illim:irlim], max_x[illim:irlim], s=40,
                    c=max_val[illim:irlim], marker='.', linewidths=(0, ),
                    clip_on=False, cmap=col)
        p.xaxis.set_ticks_position('bottom')
        p.xaxis.set_ticklabels('')
        plt.ylabel('x (km)', size=10)
        p.yaxis.set_ticks_position('right')
        p.set_xlim(llim, rlim)
        if 'true_values' in test_info:
            if not 't_err' in test_info:
                plt.hlines(x_true, llim, rlim, 'r', linewidth=2)
                plt.vlines(t_true, min(max_x), max(max_x), 'r', linewidth=2)
        else:
            plt.hlines(x[ix_true], llim, rlim, 'r', linewidth=2)
            plt.vlines(t[it_true], min(max_x), max(max_x), 'r', linewidth=2)
        if 'x_err' in test_info:
            x_low, x_high = test_info['x_err']
            plt.axhspan(x_low, x_high, facecolor='r', alpha=0.2)
        if 't_err' in test_info:
            t_left, t_right = test_info['t_err']
            plt.axvspan(t_left, t_right, facecolor='r', alpha=0.2)

    # plot max y
    if dy > 0:
        p = plt.subplot(4, 2, 6,  frameon=False)
        pos = list(p.get_position().bounds)
        fig.text(pos[0]-.05, pos[1]+pos[3], '(e)', fontsize=12)
        p.tick_params(labelsize=10)
        plt.scatter(t[illim:irlim], max_y[illim:irlim], s=40,
                    c=max_val[illim:irlim], marker='.', linewidths=(0, ),
                    clip_on=False, cmap=col)
        p.xaxis.set_ticks_position('bottom')
        p.xaxis.set_ticklabels('')
        plt.ylabel('y (km)', size=10)
        p.yaxis.set_ticks_position('right')
        p.set_xlim(llim, rlim)
        if 'true_values' in test_info:
            if not 't_err' in test_info:
                plt.hlines(y_true, llim, rlim, 'r', linewidth=2)
                plt.vlines(t_true, min(max_y), max(max_y), 'r', linewidth=2)
        else:
            plt.hlines(y[iy_true], llim, rlim, 'r', linewidth=2)
            plt.vlines(t[it_true], min(max_y), max(max_y), 'r', linewidth=2)
        if 'y_err' in test_info:
            y_low, y_high = test_info['y_err']
            plt.axhspan(y_low, y_high, facecolor='r', alpha=0.2)
        if 't_err' in test_info:
            t_left, t_right = test_info['t_err']
            plt.axvspan(t_left, t_right, facecolor='r', alpha=0.2)

    # plot max z
    if dz > 0:
        p = plt.subplot(4, 2, 8, frameon=False)
        pos = list(p.get_position().bounds)
        fig.text(pos[0]-.05, pos[1]+pos[3], '(g)', fontsize=12)
        p.tick_params(labelsize=10)
        plt.scatter(t[illim:irlim], max_z[illim:irlim], s=40,
                    c=max_val[illim:irlim], marker='.', linewidths=(0, ),
                    clip_on=False, cmap=col)
        plt.xlabel('Time (s)', size=10)
        p.xaxis.set_ticks_position('bottom')
        plt.ylabel('z (km down)', size=10)
        p.yaxis.set_ticks_position('right')
        p.set_xlim(llim, rlim)
        if 'true_values' in test_info:
            if not 't_err' in test_info:
                plt.hlines(z_true, llim, rlim, 'r', linewidth=2)
                plt.vlines(t_true, min(max_z), max(max_z), 'r', linewidth=2)
        else:
            plt.hlines(z[iz_true], llim, rlim, 'r', linewidth=2)
            plt.vlines(t[it_true], min(max_z), max(max_z), 'r', linewidth=2)
        if 'z_err' in test_info:
            z_low, z_high = test_info['z_err']
            plt.axhspan(z_low, z_high, facecolor='r', alpha=0.2)
        if 't_err' in test_info:
            t_left, t_right = test_info['t_err']
            plt.axvspan(t_left, t_right, facecolor='r', alpha=0.2)

    # add independent colorbar
    ax1 = fig.add_axes([0.40, 0.03, 0.2, 0.015])
    ax1.tick_params(labelsize=8)
    ax1.xaxis.set_ticks_position('bottom')
    cmap = mpl.cm.hot_r
    norm = mpl.colors.Normalize(vmin=np.min(max_val), vmax=np.max(max_val))
    mpl.colorbar.ColorbarBase(ax1, cmap=cmap, norm=norm,
                              orientation='horizontal',
                              ticks=[0, int(np.max(max_val)/2),
                                     int(np.max(max_val))])
    pos = list(ax1.get_position().bounds)
    fig.text(pos[0]+pos[2]/2., pos[1]+pos[3]+0.01, 'Stack max', fontsize=8,
             horizontalalignment='center', verticalalignment='bottom')

    plt.savefig(fig_filename)

    f.close()
    f_stack.close()


def gaussian(x, mu, sigma):
    return (1./(sigma*np.sqrt(2.*np.pi))) * \
        np.exp(-1*(x-mu)*(x-mu)/(2*sigma*sigma))


def plotProbLoc(marginals, prob_loc, loc, fig_dir, space_only):
    """
    Creates plot for probabilistic location
    TODO : flesh out this doc-string

    :param marginals:
    :param prob_loc:
    :param loc:
    :param fig_dir:
    :param space_only:
    """

    # get basic parameters

    x = marginals['x'][:]
    y = marginals['y'][:]
    z = marginals['z'][:]*(-1)
    z.sort()
    prob_x = marginals['prob_x'][:]
    prob_y = marginals['prob_y'][:]
    prob_z = marginals['prob_z'][:]
    prob_z_rev = prob_z[::-1]     # reverse array for plot (z is up in plot)
    prob_xy = marginals['prob_xy'][:, :]
    prob_xz = marginals['prob_xz'][:, :]
    prob_yz = marginals['prob_yz'][:, :]

    x_low = loc['x_mean'] - loc['x_sigma']
    y_low = loc['y_mean'] - loc['y_sigma']
    z_low = -loc['z_mean'] - loc['z_sigma']
    x_high = loc['x_mean'] + loc['x_sigma']
    y_high = loc['y_mean'] + loc['y_sigma']
    z_high = -loc['z_mean'] + loc['z_sigma']

    # set filename
    o_time = prob_loc['o_time']
    fig_filename = os.path.join(fig_dir, "probloc_%s.pdf" % o_time.isoformat())

    plt.figure(1, figsize=(18, 6))

    # XY plot
    left, width = 0.05, 0.2
    left_h = left + width
    bottom, height = 0.1, 0.6
    bottom_h = bottom + height

    # XY plot
    rect_2D = [left, bottom, width, height]
    rect_top = [left, bottom_h, width, 0.2]
    rect_rig = [left_h, bottom, 0.075, height]

    ax_2D = plt.axes(rect_2D)
    ax_2D.tick_params(labelsize=10)
    ax_top = plt.axes(rect_top, sharex=ax_2D)
    ax_top.tick_params(labelsize=10)
    ax_top.xaxis.set_ticks_position('top')
    ax_top.yaxis.set_ticks(())
    ax_rig = plt.axes(rect_rig, sharey=ax_2D)
    ax_rig.tick_params(labelsize=10)
    ax_rig.yaxis.set_ticks_position('right')
    ax_rig.xaxis.set_ticks(())

    ax_2D.imshow(prob_xy.T, origin='lower', interpolation='none',
                 extent=[np.min(x), np.max(x), np.min(y), np.max(y)])
    ax_top.plot(x, prob_x, 'b')
    ax_top.axvspan(x_low, x_high, facecolor='r', alpha=0.2)
    ax_rig.plot(prob_y, y, 'b')
    ax_rig.axhspan(y_low, y_high, facecolor='r', alpha=0.2)

    #XZ plot
    left, width = left_h+0.05+0.075, 0.2
    left_h = left + width

    rect_2D = [left, bottom, width, height]
    rect_top = [left, bottom_h, width, 0.2]
    rect_rig = [left_h, bottom, 0.075, height]

    ax_2D_2 = plt.axes(rect_2D)
    ax_2D_2.tick_params(labelsize=10)
    ax_top_2 = plt.axes(rect_top, sharex=ax_2D_2)
    ax_top_2.tick_params(labelsize=10)
    ax_top_2.xaxis.set_ticks_position('top')
    ax_top_2.yaxis.set_ticks(())
    ax_rig_2 = plt.axes(rect_rig, sharey=ax_2D_2)
    ax_rig_2.tick_params(labelsize=10)
    ax_rig_2.yaxis.set_ticks_position('right')
    ax_rig_2.xaxis.set_ticks(())

    ax_2D_2.imshow(prob_xz.T, origin='upper', interpolation='none',
                   extent=[np.min(x), np.max(x), np.min(z), np.max(z)])
    ax_top_2.plot(x, prob_x, 'b')
    ax_top_2.axvspan(x_low, x_high, facecolor='r', alpha=0.2)
    ax_rig_2.plot(prob_z_rev, z, 'b')
    ax_rig_2.axhspan(z_low, z_high, facecolor='r', alpha=0.2)

    #XZ plot
    left, width = left_h+0.05+0.075, 0.2
    left_h = left + width

    rect_2D = [left, bottom, width, height]
    rect_top = [left, bottom_h, width, 0.2]
    rect_rig = [left_h, bottom, 0.075, height]

    ax_2D_3 = plt.axes(rect_2D)
    ax_2D_3.tick_params(labelsize=10)
    ax_top_3 = plt.axes(rect_top, sharex=ax_2D_3)
    ax_top_3.tick_params(labelsize=10)
    ax_top_3.xaxis.set_ticks_position('top')
    ax_top_3.yaxis.set_ticks(())
    ax_rig_3 = plt.axes(rect_rig, sharey=ax_2D_3)
    ax_rig_3.tick_params(labelsize=10)
    ax_rig_3.yaxis.set_ticks_position('right')
    ax_rig_3.xaxis.set_ticks(())

    ax_2D_3.imshow(prob_yz.T, origin='upper', interpolation='none',
                   extent=[np.min(y), np.max(y), np.min(z), np.max(z)])
    ax_top_3.plot(y, prob_y, 'b')
    ax_top_3.axvspan(y_low, y_high, facecolor='r', alpha=0.2)
    ax_rig_3.plot(prob_z_rev, z, 'b')
    ax_rig_3.axhspan(z_low, z_high, facecolor='r', alpha=0.2)

    plt.savefig(fig_filename)
    plt.clf()
