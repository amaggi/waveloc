import os
import glob
import logging


class WavelocOptions(object):
    """
    The WavelocOptions class contains a single attribute **opdict** containing
    all the options and parameters to control the behaviour of waveloc.
    It also contains many methods to verify the consistency of the options
    before launching the various steps of the waveloc procedure.
    """

    def __init__(self):

        self.opdict = {}

        # set some default values

        # general profiling / debugging behaviour
        self.opdict['time'] = False
        self.opdict['verbose'] = False
        self.opdict['use_ram'] = False

        # data processing
        self.opdict['resample'] = False
        self.opdict['krec'] = False
        self.opdict['kderiv'] = False
        self.opdict['gauss'] = False

        # migration
        self.opdict['ugrid_type'] = 'FULL'
        self.opdict['load_ttimes_buf'] = True
        self.opdict['reloc'] = False
        self.opdict['reloc_snr'] = 12.

        # location
        self.opdict['auto_loclevel'] = False
        self.opdict['loclevel'] = 50.
        self.opdict['snr_loclevel'] = 10.
        self.opdict['snr_limit'] = 10.
        self.opdict['snr_tr_limit'] = 10.
        self.opdict['sn_time'] = 10.
        self.opdict['n_kurt_min'] = 4

        # prob density location
        self.opdict['probloc_spaceonly'] = True

        # synthetic
        self.opdict['syn_addnoise'] = False
        self.opdict['syn_amplitude'] = 1.
        self.opdict['syn_kwidth'] = 0.1

        # cross-correlation
        self.opdict['xcorr_threshold'] = 0.7
        self.opdict['xcorr_before'] = 0.5
        self.opdict['xcorr_after'] = 6.0

        # clustering
        self.opdict['nbsta'] = 3
        self.opdict['clus'] = 0.8

        # double-difference
        self.opdict['dd_loc'] = False

        # kurtogram
        self.opdict['new_kurtfile'] = False

    def set_test_options(self):
        """
        Sets up a standard version of opdict. Used notably in the example
        cases.
        """
        self.opdict['time'] = True
        self.opdict['verbose'] = True

        self.opdict['test_datadir'] = 'test_data'
        self.opdict['datadir'] = 'TEST'
        self.opdict['outdir'] = 'TEST'

        self.opdict['net_list'] = 'YA'
        self.opdict['sta_list'] = "FJS,FLR,FOR,HDL,RVL,SNE,UV01,UV02,UV03,\
UV04,UV05,UV06,UV07,UV08,UV09,UV10,UV11,UV12,UV13,UV14,UV15"
        self.opdict['comp_list'] = "HHZ"

        self.opdict['starttime'] = "2010-10-14T00:14:00.0Z"
        self.opdict['endtime'] = "2010-10-14T00:18:00.0Z"

        self.opdict['time_grid'] = 'Slow_len.100m.P'
        self.opdict['search_grid'] = 'test_grid.search.hdr'
        self.opdict['stations'] = 'coord_stations_test'

        self.opdict['resample'] = False
        self.opdict['fs'] = None

        self.opdict['c1'] = 4.0
        self.opdict['c2'] = 10.0

        self.opdict['kwin'] = 4
        self.opdict['krec'] = False
        self.opdict['kderiv'] = True
        self.opdict['gauss'] = False
        self.opdict['gthreshold'] = 0.1
        self.opdict['mu'] = 0
        self.opdict['sigma'] = 0.1

        self.opdict['data_length'] = 600
        self.opdict['data_overlap'] = 20

        self.opdict['dataglob'] = '*filt.mseed'
        self.opdict['kurtglob'] = '*kurt.mseed'
        self.opdict['gradglob'] = '*grad.mseed'
        self.opdict['gaussglob'] = '*gauss.mseed'

        self.opdict['load_ttimes_buf'] = True

        self.opdict['reloc'] = False
        self.opdict['reloc_snr'] = 12.

        self.opdict['auto_loclevel'] = False
        self.opdict['loclevel'] = 50.0
        self.opdict['snr_limit'] = 10.0
        self.opdict['snr_tr_limit'] = 10.0
        self.opdict['sn_time'] = 10.0
        self.opdict['n_kurt_min'] = 4

        self.opdict['syn_addnoise'] = False

        self.opdict['new_kurtfile'] = False

        self.opdict['xcorr_threshold'] = 0.7
        self.opdict['xcorr_before'] = 0.5
        self.opdict['xcorr_after'] = 6.0
        self.opdict['xcorr_corr'] = 'corr'
        self.opdict['xcorr_delay'] = 'delay'

        self.opdict['clus'] = 0.8
        self.opdict['nbsta'] = 3

        self.opdict['dd_loc'] = True

    def verify_base_path(self):
        """
        Verifies that the base_path is set. If the 'base_path' option is not
        directly set in the opdict, its value is read from the environment
        variable $WAVELOC_PATH.
        """

        # if the option base_path is not set, then check the environment
        # variable
        # if the environment variable is not set, quit with error message
        if 'base_path' not in self.opdict:
            logging.info('No base_path set in options, getting base_path from \
                          $WAVELOC_PATH')
            base_path = os.getenv('WAVELOC_PATH')
            if not os.path.isdir(base_path):
                raise UserWarning('Environment variable WAVELOC_PATH not set \
                  correctly.')
            self.opdict['base_path'] = base_path

        base_path = self.opdict['base_path']
        lib_path = os.path.join(base_path, 'lib')
        if not os.path.isdir(lib_path):
            raise UserWarning('Directory %s does not exist.' % lib_path)

    def _verify_lib_path(self):
        self.verify_base_path()
        base_path = self.opdict['base_path']
        lib_path = os.path.join(base_path, 'lib')
        if not os.path.isdir(lib_path):
            raise UserWarning('Directory %s does not exist.' % lib_path)

    def _verify_datadir(self):
        self.verify_base_path()
        base_path = self.opdict['base_path']
        if 'datadir' not in self.opdict:
            raise UserWarning('datadir option not set')

        datadir = os.path.join(base_path, 'data', self.opdict['datadir'])
        if not os.path.isdir(datadir):
            raise UserWarning('Directory %s does not exist.' % datadir)

    def _verify_outdir(self):
        self.verify_base_path()
        base_path = self.opdict['base_path']
        if 'outdir' not in self.opdict:
            raise UserWarning('outdir option not set')

        outdir = os.path.join(base_path, 'out', self.opdict['outdir'])
        if not os.path.isdir(outdir):
            os.makedirs(outdir)
        if not os.path.isdir(os.path.join(outdir, 'fig')):
            os.makedirs(os.path.join(outdir, 'fig'))
        if not os.path.isdir(os.path.join(outdir, 'grid')):
            os.makedirs(os.path.join(outdir, 'grid'))
        if not os.path.isdir(os.path.join(outdir, 'loc')):
            os.makedirs(os.path.join(outdir, 'loc'))
        if not os.path.isdir(os.path.join(outdir, 'stack')):
            os.makedirs(os.path.join(outdir, 'stack'))
        if not os.path.isdir(os.path.join(outdir, 'time_grids')):
            os.makedirs(os.path.join(outdir, 'time_grids'))
        if self.opdict['reloc'] and not os.path.isdir(os.path.join(outdir,
                                                                   'reloc')):
            os.makedirs(os.path.join(outdir, 'reloc'))

    def _verify_net_list(self):
        if 'net_list' not in self.opdict:
            raise UserWarning('net_list option not set')

    def _verify_sta_list(self):
        if not 'sta_list'in self.opdict:
            raise UserWarning('sta_list option not set')

    def _verify_comp_list(self):
        if 'comp_list' not in self.opdict:
            raise UserWarning('comp_list option not set')

    def _verify_channel_file(self):
        if 'channel_file' not in self.opdict:
            raise UserWarning('channel_file option not set')
        self._verify_lib_path()
        base_path = self.opdict['base_path']
        filename = os.path.join(base_path, 'lib', self.opdict['channel_file'])
        if not os.path.isfile(filename):
            raise UserWarning('Cannot find %s' % filename)

    def _verify_starttime(self):
        if 'starttime' not in self.opdict:
            raise UserWarning('starttime option not set')

    def _verify_endtime(self):
        if 'endtime' not in self.opdict:
            raise UserWarning('endtime option not set')

    def _verify_resample(self):
        if 'resample' not in self.opdict:
            raise UserWarning('resample option not set')

    def _verify_fs(self):
        self._verify_resample()
        resample = self.opdict['resample']
        if resample:
            if 'fs' not in self.opdict:
                raise UserWarning('fs option not set')

    def _verify_c1(self):
        if 'c1' not in self.opdict:
            raise UserWarning('c1 option not set')

    def _verify_c2(self):
        if 'c2' not in self.opdict:
            raise UserWarning('c2 option not set')

    def _verify_kwin(self):
        if 'kwin' not in self.opdict:
            raise UserWarning('kwin option not set')

    def _verify_gthreshold(self):
        if 'gthreshold' not in self.opdict:
            raise UserWarning('gthreshold option not set')

    def _verify_mu(self):
        if 'mu' not in self.opdict:
            raise UserWarning('mu option not set')

    def _verify_sigma(self):
        if 'sigma' not in self.opdict:
            raise UserWarning('sigma option not set')

    def _verify_dataless(self):
        if 'dataless' not in self.opdict:
            raise UserWarning('dataless option not set')
        self._verify_lib_path()
        base_path = self.opdict['base_path']
        lib_path = os.path.join(base_path, 'lib')
        dataless_names = glob.glob(os.path.join(lib_path,
                                                self.opdict['dataless']))
        if len(dataless_names) == 0:
            raise UserWarning('No dataless files found: %s' % dataless_names)

    def _verify_dataglob(self):
        if 'dataglob' not in self.opdict:
            raise UserWarning('dataglob option not set')
        self._verify_datadir()
        base_path = self.opdict['base_path']
        datadir = os.path.join(base_path, 'data', self.opdict['datadir'])
        data_names = glob.glob(os.path.join(datadir, self.opdict['dataglob']))
        if len(data_names) == 0:
            raise UserWarning('No data files found : %s' % data_names)

    def _verify_kurtglob(self):
        if 'kurtglob' not in self.opdict:
            raise UserWarning('kurtglob option not set')
        self._verify_datadir()
        base_path = self.opdict['base_path']
        datadir = os.path.join(base_path, 'data', self.opdict['datadir'])
        kurt_names = glob.glob(os.path.join(datadir, self.opdict['kurtglob']))
        if len(kurt_names) == 0:
            raise UserWarning('No kurtosis files found : %s' % kurt_names)

    def _verify_gradglob(self):
        if 'gradglob' not in self.opdict:
            raise UserWarning('gradglob option not set')
        self._verify_datadir()
        base_path = self.opdict['base_path']
        datadir = os.path.join(base_path, 'data', self.opdict['datadir'])
        grad_names = glob.glob(os.path.join(datadir, self.opdict['gradglob']))
        if len(grad_names) == 0:
            raise UserWarning('No kurtosis gradient files found : %s' %
                              grad_names)

    def _verify_gaussglob(self):
        if 'gaussglob' not in self.opdict:
            raise UserWarning('gaussglob option not set')
        self._verify_datadir()
        base_path = self.opdict['base_path']
        datadir = os.path.join(base_path, 'data', self.opdict['datadir'])
        gauss_names = glob.glob(os.path.join(datadir,
                                self.opdict['gaussglob']))
        if len(gauss_names) == 0:
            raise UserWarning('No gaussian files found : %s' % gauss_names)

    def _verify_ugrid_type(self):
        if 'ugrid_type' not in self.opdict:
            raise UserWarning('ugrid_type option not set')
        utype = self.opdict['ugrid_type']
        if utype not in ['FULL', 'USER']:
            raise UserWarning('Unknown ugrid_type %s' % utype)
        if utype == 'USER':
            self._verify_ugrid_file()
        else:
            self._verify_search_grid()

    def _verify_time_grid(self):
        if 'time_grid' not in self.opdict:
            raise UserWarning('time_grid option not set')
        self._verify_lib_path()
        base_path = self.opdict['base_path']
        time_grid = os.path.join(base_path, 'lib', self.opdict['time_grid'])
        tg_glob = time_grid+'*'
        tg_files = glob.glob(tg_glob)
        if len(tg_files) == 0:
            raise UserWarning('No time grid files found %s' % tg_glob)

    def _verify_data_length(self):
        if 'data_length' not in self.opdict:
            raise UserWarning('data_length option not set')

    def _verify_data_overlap(self):
        if 'data_overlap' not in self.opdict:
            raise UserWarning('data_overlap option not set')

    def _verify_snr_limit(self):
        if 'snr_limit' not in self.opdict:
            raise UserWarning('snr_limit option not set')

    def _verify_snr_tr_limit(self):
        if 'snr_tr_limit' not in self.opdict:
            raise UserWarning('snr_tr_limit option not set')

    def _verify_sn_time(self):
        if 'sn_time' not in self.opdict:
            raise UserWarning('sn_time option not set')

    def _verify_n_kurt_min(self):
        if 'n_kurt_min' not in self.opdict:
            raise UserWarning('n_kurt_min option not set')

    def _verify_stations(self):
        if 'stations' not in self.opdict:
            raise UserWarning('stations option not set')
        self._verify_lib_path()
        base_path = self.opdict['base_path']
        stations = os.path.join(base_path, 'lib', self.opdict['stations'])
        if not os.path.isfile(stations):
            raise UserWarning('Cannot find %s' % stations)

    def _verify_search_grid(self):
        if 'search_grid' not in self.opdict:
            raise UserWarning('search_grid option not set')
        self._verify_lib_path()
        base_path = self.opdict['base_path']
        search_grid = os.path.join(base_path, 'lib',
                                   self.opdict['search_grid'])
        if not os.path.isfile(search_grid):
            raise UserWarning('Cannot find %s' % search_grid)

    def _verify_ugrid_file(self):
        if 'ugrid_file' not in self.opdict:
            raise UserWarning('ugrid_file option not set')
        self._verify_lib_path()
        base_path = self.opdict['base_path']
        search_grid = os.path.join(base_path, 'lib',
                                   self.opdict['ugrid_file'])
        if not os.path.isfile(search_grid):
            raise UserWarning('Cannot find %s' % search_grid)

    def _verify_reloc(self):
        if 'reloc' not in self.opdict:
            raise UserWarning('reloc option not set')

    def _verify_reloc_snr(self):
        if 'reloc_snr' not in self.opdict:
            raise UserWarning('reloc_snr option not set')

    def _verify_auto_loclevel(self):
        if 'auto_loclevel' not in self.opdict:
            raise UserWarning('auto_loclevel option not set')

    def _verify_snr_loclevel(self):
        self._verify_auto_loclevel()
        auto_loclevel = self.opdict['auto_loclevel']
        if auto_loclevel:
            if 'snr_loclevel' not in self.opdict:
                raise UserWarning('snr_loclevel option not set')

    def _verify_loclevel(self):
        self._verify_auto_loclevel()
        auto_loclevel = self.opdict['auto_loclevel']
        if not auto_loclevel:
            if 'loclevel' not in self.opdict:
                raise UserWarning('loclevel option not set')

    def _verify_probloc_spaceonly(self):
        if 'probloc_spaceonly' not in self.opdict:
            raise UserWarning('probloc_spaceonly option not set')

    def _verify_xcorr_threshold(self):
        if 'xcorr_threshold' not in self.opdict:
            raise UserWarning('xcorr_threshold option not set')

    def _verify_newkurtfile(self):
        if 'new_kurtfile' not in self.opdict:
            raise UserWarning('new_kurtfile option not set')

    def _verify_xcorr_before(self):
        if 'xcorr_before' not in self.opdict:
            raise UserWarning('xcorr_before option not set')

    def _verify_xcorr_after(self):
        if 'xcorr_after' not in self.opdict:
            raise UserWarning('xcorr_after option not set')

    def _verify_xcorr_corr(self):
        if 'xcorr_corr' not in self.opdict:
            raise UserWarning('xcorr_corr option not set')

    def _verify_xcorr_delay(self):
        if 'xcorr_delay' not in self.opdict:
            raise UserWarning('xcorr_delay option not set')

    def _verify_nbsta(self):
        if 'nbsta' not in self.opdict:
            raise UserWarning('nbsta option not set')

    def _verify_clus(self):
        if 'clus' not in self.opdict:
            raise UserWarning('clus option not set')

    def _verify_dd_loc(self):
        if 'dd_loc' not in self.opdict:
            raise UserWarning('dd_loc option not set')

    def _verify_syn_addnoise(self):
        if 'syn_addnoise' not in self.opdict:
            raise UserWarning('syn_addnoise option not set')

    def _verify_syn_snr(self):
        self._verify_syn_addnoise()
        syn_addnoise = self.opdict['syn_addnoise']
        if syn_addnoise:
            if 'syn_snr' not in self.opdict:
                raise UserWarning('syn_snr option not set')

    def _verify_syn_amplitude(self):
        if 'syn_amplitude' not in self.opdict:
            raise UserWarning('syn_amplitude option not set')

    def _verify_syn_datalength(self):
        if 'syn_datalength' not in self.opdict:
            raise UserWarning('syn_datalength option not set')

    def _verify_syn_samplefreq(self):
        if 'syn_samplefreq' not in self.opdict:
            raise UserWarning('syn_samplefreq option not set')

    def _verify_syn_kwidth(self):
        if 'syn_kwidth' not in self.opdict:
            raise UserWarning('syn_kwidth option not set')

    def _verify_syn_otime(self):
        if 'syn_otime' not in self.opdict:
            raise UserWarning('syn_otime option not set')

    def _verify_syn_x(self):
        if 'syn_x' not in self.opdict:
            raise UserWarning('syn_x option not set')

    def _verify_syn_y(self):
        if 'syn_y' not in self.opdict:
            raise UserWarning('syn_y option not set')

    def _verify_syn_z(self):
        if 'syn_z' not in self.opdict:
            raise UserWarning('syn_z option not set')

    def _verify_syn_filename(self):
        if 'syn_filename' not in self.opdict:
            raise UserWarning('syn_filename option not set')

    def _verify_plot_tbefore(self):
        if 'plot_tbefore' not in self.opdict:
            raise UserWarning('plot_tbefore option not set')

    def _verify_plot_tafter(self):
        if 'plot_tafter' not in self.opdict:
            raise UserWarning('plot_tafter option not set')

    def _verify_otime_window(self):
        if 'otime_window' not in self.opdict:
            raise UserWarning('otime_window option not set')

    def _verify_channel_net_sta_comp(self):
        # if have channel_file option, check that
        if 'channel_file' in self.opdict:
            self._verify_channel_file()
        # else check net sta comp lists are set
        else:
            self._verify_net_list()
            self._verify_sta_list()
            self._verify_comp_list()

    def verify_SDS_processing_options(self):
        """
        Verify presence of all options necessary for SDS_processing.
        """

        self.verify_base_path()
        self._verify_datadir()

        self._verify_channel_net_sta_comp()

        self._verify_starttime()
        self._verify_endtime()

        self._verify_fs()
        self._verify_c1()
        self._verify_c2()
        self._verify_kwin()

        if self.opdict['gauss']:
            self._verify_gthreshold()
            self._verify_mu()
            self._verify_sigma()

    def verify_migration_options(self):
        """
        Verify presence of all options necessary for migration.
        """

        self.verify_base_path()
        self._verify_lib_path()
        self._verify_datadir()
        self._verify_outdir()

        self._verify_channel_net_sta_comp()

        self._verify_kurtglob()
        if self.opdict['kderiv']:
            self._verify_gradglob()
            if self.opdict['gauss']:
                self._verify_gaussglob()

        self._verify_starttime()
        self._verify_endtime()
        self._verify_data_length()
        self._verify_data_overlap()

        self._verify_stations()
        self._verify_ugrid_type()
        self._verify_search_grid()
        self._verify_time_grid()

        self._verify_reloc()
        if self.opdict['reloc']:
            self._verify_reloc_snr()

    def verify_location_options(self):
        """
        Verify presence of all options necessary for location.
        """

        self.verify_base_path()
        self._verify_lib_path()
        self._verify_datadir()
        self._verify_outdir()

        self._verify_kurtglob()

        self._verify_loclevel()
        self._verify_snr_loclevel()
        self._verify_snr_limit()
        self._verify_snr_tr_limit()
        self._verify_sn_time()
        self._verify_n_kurt_min()
        self._verify_probloc_spaceonly()

        self._verify_search_grid()
        self._verify_time_grid()

        self._verify_reloc()

    def verify_kurtogram_options(self):
        """
        Verify presence of all options necessary for kurtogram.
        """

        self.verify_base_path()
        self._verify_datadir()
        self._verify_outdir()

        self._verify_dataglob()
        self._verify_kurtglob()

        base_path = self.opdict['base_path']
        locdir = os.path.join(base_path, 'out', self.opdict['outdir'], 'loc')
        locfile = os.path.join(locdir, 'locations.dat')
        if not os.path.isfile(locfile):
            raise UserWarning('Cannot find %s' % locfile)

        self._verify_newkurtfile()

    def verify_magnitude_options(self):
        """
        Verify presence of all options necessary for magnitude.
        """

        self.verify_base_path()
        self._verify_lib_path()
        self._verify_datadir()
        self._verify_outdir()

        self._verify_channel_net_sta_comp()
        self._verify_dataless()

        base_path = self.opdict['base_path']
        locdir = os.path.join(base_path, 'out', self.opdict['outdir'], 'loc')
        locfile = os.path.join(locdir, 'locations.dat')
        if not os.path.isfile(locfile):
            raise UserWarning('Cannot find %s' % locfile)

    def verify_correlation_options(self):
        """
        Verify presence of all options necessary for correlation.
        """

        self.verify_base_path()
        self._verify_lib_path()
        self._verify_datadir()
        self._verify_outdir()

        self._verify_dataglob()
        self._verify_xcorr_threshold()
        self._verify_xcorr_before()
        self._verify_xcorr_after()
        self._verify_xcorr_corr()
        self._verify_xcorr_delay()

        base_path = self.opdict['base_path']
        locdir = os.path.join(base_path, 'out', self.opdict['outdir'], 'loc')
        locfile = os.path.join(locdir, 'locations.dat')
        if not os.path.isfile(locfile):
            raise UserWarning('Cannot find %s' % locfile)

    def verify_cluster_options(self):
        """
        Verify presence of all options necessary for clustering.
        """

        self.verify_base_path()
        self._verify_lib_path()
        self._verify_datadir()
        self._verify_outdir()

        base_path = self.opdict['base_path']
        locdir = os.path.join(base_path, 'out', self.opdict['outdir'], 'loc')

        self._verify_dataglob()

        self._verify_stations()
        self._verify_xcorr_corr()
        self._verify_xcorr_delay()

        coeff_file = os.path.join(locdir, self.opdict['xcorr_corr'])
        if not os.path.isfile(coeff_file):
            raise UserWarning('Cannot find %s' % coeff_file)

        delay_file = os.path.join(locdir, self.opdict['xcorr_delay'])
        if not os.path.isfile(delay_file):
            raise UserWarning('Cannot find %s' % delay_file)

        self._verify_nbsta()
        self._verify_clus()

    def verify_doublediff_options(self):
        """
        Verify presence of all options necessary for double difference
        location.
        """

        self.verify_base_path()
        self._verify_lib_path()
        self._verify_outdir()
        base_path = self.opdict['base_path']

        self._verify_time_grid()
        self._verify_search_grid()

        locdir = os.path.join(base_path, 'out', self.opdict['outdir'], 'loc')

        self._verify_stations()
        self._verify_xcorr_corr()
        self._verify_xcorr_delay()

        coeff_file = os.path.join(locdir, self.opdict['xcorr_corr'])
        if not os.path.isfile(coeff_file):
            raise UserWarning('Cannot find %s' % coeff_file)

        delay_file = os.path.join(locdir, self.opdict['xcorr_delay'])
        if not os.path.isfile(delay_file):
            raise UserWarning('Cannot find %s' % delay_file)

        self._verify_nbsta()
        self._verify_clus()

        self._verify_dd_loc()

    def verify_synthetic_options(self):
        """
        Verify presence of all options necessary for synthetic migration runs.
        """

        self.verify_base_path()
        self._verify_lib_path()
        self._verify_outdir()

        self._verify_time_grid()

        self._verify_stations()
        self._verify_syn_addnoise()
        self._verify_syn_snr()

        self._verify_syn_amplitude()
        self._verify_syn_datalength()
        self._verify_syn_samplefreq()
        self._verify_syn_kwidth()
        self._verify_syn_otime()
        self._verify_syn_x()
        self._verify_syn_y()
        self._verify_syn_z()
        self._verify_syn_filename()

    def verify_plotting_options(self):
        """
        Verify presence of all options necessary for plotting waveloc results.
        """

        self.verify_base_path()
        self._verify_lib_path()
        self._verify_datadir()
        self._verify_outdir()

        base_path = self.opdict['base_path']
        locdir = os.path.join(base_path, 'out', self.opdict['outdir'], 'loc')

        locfile = os.path.join(locdir, 'locations.dat')
        if not os.path.isfile(locfile):
            raise UserWarning('Locations file %s does not exist.' % locfile)

        self._verify_dataglob()
        self._verify_kurtglob()
        self._verify_gradglob()

        self._verify_plot_tbefore()
        self._verify_plot_tafter()
        self._verify_otime_window()

        self._verify_search_grid()
        self._verify_time_grid()

        self._verify_stations()

    def verify_probloc_plotting_options(self):
        """
        Verify presence of all options necessary for plotting waveloc using the
        probability density location method.
        """

        self.verify_base_path()
        self._verify_lib_path()
        self._verify_datadir()
        self._verify_outdir()

        base_path = self.opdict['base_path']
        locdir = os.path.join(base_path, 'out', self.opdict['outdir'], 'loc')

        locfile = os.path.join(locdir, 'locations.dat')
        if not os.path.isfile(locfile):
            raise UserWarning('Locations file %s does not exist.' % locfile)

        locfile = os.path.join(locdir, 'locations_prob.dat')
        if not os.path.isfile(locfile):
            raise UserWarning('Locations file %s does not exist.' % locfile)

        locfile = os.path.join(locdir, 'locations_prob.hdf5')
        if not os.path.isfile(locfile):
            raise UserWarning('Locations file %s does not exist.' % locfile)
