import os
from ugrids import read_ugrid, nll2reg_ugrid


class PlotOptions(object):
    """
    The PlotOptions class contains a single attribute **opdict** containing all
    the options and parameters to control the plotting of waveloc results.
    """

    def __init__(self, wo_opdict, syn=False):

        self.opdict = {}

        # paths
        for key in 'base_path', 'outdir':
            self.opdict[key] = wo_opdict[key]

        # if not synthetic test, also need the datadir
        if not syn:
            for key in 'datadir', :
                self.opdict[key] = wo_opdict[key]

        # the waveforms that are needed for the plot
        if wo_opdict['gauss']:
            self.opdict['ktype'] = 'gauss'
        elif wo_opdict['kderiv']:
            self.opdict['ktype'] = 'kgrad'
        else:
            self.opdict['ktype'] = 'kurt'

        # ugrid stuff
        for key in 'ugrid_type', 'ugrid_file', 'search_grid':
            if key in wo_opdict:
                self.opdict[key] = wo_opdict[key]

        # other stuff
        for key in 'otime_window', :
            if key in wo_opdict:
                self.opdict[key] = wo_opdict[key]

    def getXYZ(self):

        ugrid_type = self.opdict['ugrid_type']
        base_path = self.opdict['base_path']

        if ugrid_type == 'USER':
            ugrid_file = self.opdict['ugrid_file']
            filename = os.path.join(base_path, 'lib', ugrid_file)
            x, y, z = read_ugrid(filename)
        elif ugrid_type == 'FULL':
            search_grid = self.opdict['search_grid']
            filename = os.path.join(base_path, 'lib', search_grid)
            x, y, z = nll2reg_ugrid(filename)

        return x, y, z

    def getFigDir(self):
        return os.path.join(self.opdict['base_path'], 'out',
                            self.opdict['outdir'], 'fig')

    def getGridDir(self):
        return os.path.join(self.opdict['base_path'], 'out',
                            self.opdict['outdir'], 'grid')

    def getStackDir(self):
        return os.path.join(self.opdict['base_path'], 'out',
                            self.opdict['outdir'], 'stack')

    def getGridFilename(self):
        dir = self.getGridDir()
        return os.path.join(dir, self.opdict['grid_filename'])

    def getStackFilename(self):
        dir = self.getStackDir()
        return os.path.join(dir, self.opdict['stack_filename'])

    def getFigFilename(self):
        grid_filename = self.getGridFilename()
        basename = os.path.basename(grid_filename)
        fig_filename = "%s_grid.pdf" % os.path.splitext(basename)[0]
        dir = self.getFigDir()
        return os.path.join(dir, fig_filename)

    def getWfmFigFilename(self):
        grid_filename = self.getGridFilename()
        basename = os.path.basename(grid_filename)
        fig_filename = "%s_wfm.pdf" % os.path.splitext(basename)[0]
        dir = self.getFigDir()
        return os.path.join(dir, fig_filename)
