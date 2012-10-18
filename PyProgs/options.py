import os, glob, argparse, logging

class WavelocOptions(object):
  """
  Describe the WavelocOptions class, and all the options
  """

  def __init__(self):

    self.opdict={}


    # set some default values

    # general profiling / debugging behaviour
    self.opdict['time']=False
    self.opdict['verbose']=False

    # data processing
    self.opdict['resample']=False
    self.opdict['krec']=False
    self.opdict['kderiv']=False

    # migration
    self.opdict['load_ttimes_buf']=True
    self.opdict['reloc']=False
    self.opdict['reloc_snr']=12.

    # location
    self.opdict['auto_loclevel']=False
    self.opdict['loclevel']=50.
    self.opdict['snr_loclevel']=10.
    self.opdict['snr_limit']=10.
    self.opdict['snr_tr_limit']=10.
    self.opdict['sn_time']=10.
    self.opdict['n_kurt_min']=4

    # synthetic
    self.opdict['syn_addnoise']=False
    self.opdict['syn_amplitude']=1.
    self.opdict['syn_kwidth']=0.1

    # cross-correlation
    self.opdict['xcorr_threshold']=0.7
    self.opdict['xcorr_before']=0.5
    self.opdict['xcorr_after']=6.0

    # clustering
    self.opdict['nbsta']=3
    self.opdict['clus']=0.8

    # double-difference
    self.opdict['dd_loc']=False

  
    # For now, continue to support command-line arguments
    # TODO : get rid of these evenutally for ease of maintenance
    self.p = argparse.ArgumentParser()

    self.p.add_argument('--time', '-t', action='store_true',
            default=self.opdict['time'], 
            help='print timing information to stout')
    self.p.add_argument('--verbose', '-v', action='store_true',
            default=self.opdict['verbose'], 
            help='print debugging information to stout')
  
    self.p.add_argument('--datadir',action='store', 
            help="subdirectory of base_path/data")
    self.p.add_argument('--outdir', action='store', 
            help='subdirectory of base_path/out for stocking output files')

    self.p.add_argument('--net_list', action='store', 
            help="list of network codes (e.g. \"BE,G\") ")
    self.p.add_argument('--sta_list', action='store',
            help="list of station names (e.g. \"STA1,STA2\") ")
    self.p.add_argument('--comp_list',action='store',
            help="list of component names (e.g. \"HHZ,LHZ\") ")

    self.p.add_argument('--resample',action='store_true',
            default=self.opdict['resample'], help="resample data")
    self.p.add_argument('--fs',      action='store', type=float,
            help="resample frequency")

    self.p.add_argument('--c1',action='store',type=float,  
            help="low frequency corner of band pass filter ")
    self.p.add_argument('--c2',action='store',type=float,  
            help="high frequency corner of band pass filter ")
    self.p.add_argument('--kwin',action='store',type=float, 
            help="length of kurtosis window (seconds)")
    self.p.add_argument('--krec',action='store_true',
            default=self.opdict['krec'], 
            help="use recursive kurtosis calculation (faster but less precise)")
    self.p.add_argument('--kderiv',action='store_true',
            default=self.opdict['kderiv'], help="use derivative of kurtosis")

    self.p.add_argument('--dataglob',action='store',help="data glob")
    self.p.add_argument('--kurtglob',action='store',help="kurtosis glob")
    self.p.add_argument('--gradglob',action='store',help="gradient glob")

    self.p.add_argument('--starttime', action='store', 
            help="start time for data e.g. 2010-10-14T00:00:00.0Z")
    self.p.add_argument('--endtime',  action='store', 
            help="end time for data e.g. 2010-10-14T10:00:00.0Z")
    self.p.add_argument('--data_length', action='store', type=float,
            help="length in seconds for data segments to analyse (e.g. 630)")
    self.p.add_argument('--data_overlap', action='store', type=float,
            help="length in seconds for overlapping data segments (e.g. 30)")

    self.p.add_argument('--stations',action='store',
            help='station list (found in base_path/lib)') 
    self.p.add_argument('--search_grid',action='store',
            help="search grid (found in base_path/lib)")
    self.p.add_argument('--time_grid',  action='store',
            help="time grid basename (found in base_path/lib)")
    self.p.add_argument('--load_ttimes_buf',action='store_true',
            default=self.opdict['load_ttimes_buf'], help = 
            'load pre-calculated travel-times for the search grid from file')

    self.p.add_argument('--reloc', action='store_true',
            default=self.opdict['reloc'], help='select data for migration')
    self.p.add_argument('--reloc_snr', action='store',
            default=self.opdict['reloc_snr'], help='signal to noise ratio level over which kurtosis are migrated ')

    self.p.add_argument('--auto_loclevel', action='store',
            default=self.opdict['auto_loclevel'], type=float,
            help='automatically set trigger stack level for locations ')
    self.p.add_argument('--loclevel', action='store',
            default=self.opdict['loclevel'],   type=float,
            help='trigger stack level for locations (e.g. 50) ')
    self.p.add_argument('--snr_loclevel', action='store',
            default=self.opdict['snr_loclevel'],   type=float, help= 
            'SNR for automatically setting trigger stack level for locations')
    self.p.add_argument('--snr_limit',action='store',
            default=self.opdict['snr_limit'], type=float,
            help="signal_to_noise level for kurtosis acceptance")
    self.p.add_argument('--snr_tr_limit',action='store',
            default=self.opdict['snr_tr_limit'], type=float,
            help="signal_to_noise level for data acceptance")
    self.p.add_argument('--sn_time',action='store',
            default=self.opdict['sn_time'], type=float, help="time over which \
            to calculate the signal_to_noise ratio for kurtosis acceptance")
    self.p.add_argument('--n_kurt_min',action='store',
            default=self.opdict['n_kurt_min'], type=int,  
            help="min number of good kurtosis traces for a location")

    self.p.add_argument('--syn_addnoise',action='store_true',
            default=self.opdict['syn_addnoise'], 
            help="add noise to synthetic tests")
    self.p.add_argument('--syn_snr',action='store',type=float, 
            help="Signal to noise ratio for synthetic tests")
    self.p.add_argument('--syn_amplitude',action='store',type=float,
            default=self.opdict['syn_amplitude'], 
            help="amplitude of kurtosis gradient peak on synthetic waveforms")
    self.p.add_argument('--syn_datalength',action='store',type=float,
            help="length of synthetic waveforms")
    self.p.add_argument('--syn_samplefreq',action='store',type=float,
            help="sample frequency (Hz) of synthetic waveforms")
    self.p.add_argument('--syn_kwidth',action='store',type=float,
            default=self.opdict['syn_kwidth'], 
            help="width of kurtosis gradient pulse on synthetic waveforms")
    self.p.add_argument('--syn_otime',action='store',type=float, help=
            "origin time for synthetic waveforms (wrt start of waveforms)")
    self.p.add_argument('--syn_ix',action='store',type=int, 
            help="x grid index for syntetic hypocenter")
    self.p.add_argument('--syn_iy',action='store',type=int, 
            help="y grid index for syntetic hypocenter")
    self.p.add_argument('--syn_iz',action='store',type=int, 
            help="z grid index for syntetic hypocenter")
    self.p.add_argument('--syn_filename',action='store', 
            help="filename for synthetic grid")

    self.p.add_argument('--plot_tbefore',action='store',type=float, 
            help="time before origin time for plots")
    self.p.add_argument('--plot_tafter',action='store',type=float, 
            help="time after origin time for plots")

    self.p.add_argument('--xcorr_threshold',action='store',
            default=self.opdict['xcorr_threshold'], type=float, 
            help="correlation value over which the correlation is computed \
                    again in the Fourier domain")
    self.p.add_argument('--xcorr_before',action='store',
            default=self.opdict['xcorr_before'], type=float, help=
            "cross-correlation window: time interval before the origin time")
    self.p.add_argument('--xcorr_after',action='store', default=self.opdict['xcorr_after'],
            type=float, help="cross-correlation window: time interval after \
                    the origin time")
    self.p.add_argument('-xcorr_-corr', action='store', 
            help="name of the file containing all correlation values")
    self.p.add_argument('--xcorr_delay', action='store', 
            help="name of the file containing all time delays")

    self.p.add_argument('--nbsta',action='store', default=self.opdict['nbsta'],
            type=int, help="number of stations over which an event pair is \
            considered provided that its correlation coefficient is greater \
            than a given threshold")
    self.p.add_argument('--clus',action='store', default=self.opdict['clus'],
            type=float, 
            help="correlation value over which an event pair is considered")

    self.p.add_argument('--dd_loc',action='store_true',help="create a new location file with double difference locations")

  def set_all_arguments(self,args):
    self.opdict['time']=args.time
    self.opdict['verbose']=args.verbose

    self.opdict['datadir']=args.datadir
    self.opdict['outdir']=args.outdir

    self.opdict['net_list']=args.net_list
    self.opdict['sta_list']=args.sta_list
    self.opdict['comp_list']=args.comp_list

    self.opdict['resample']=args.resample
    self.opdict['fs']=args.fs

    self.opdict['c1']=args.c1
    self.opdict['c2']=args.c2
    self.opdict['kwin']=args.kwin
    self.opdict['krec']=args.krec
    self.opdict['kderiv']=args.kderiv

    self.opdict['dataglob']=args.dataglob
    self.opdict['kurtglob']=args.kurtglob
    self.opdict['gradglob']=args.gradglob

    self.opdict['starttime']=args.starttime
    self.opdict['endtime']=args.endtime
    self.opdict['data_length']=args.data_length
    self.opdict['data_overlap']=args.data_overlap

    self.opdict['stations']=args.stations
    self.opdict['search_grid']=args.search_grid
    self.opdict['time_grid']=args.time_grid
    self.opdict['load_ttimes_buf']=args.load_ttimes_buf

    self.opdict['reloc']=args.reloc
    self.opdict['reloc_snr']=args.reloc_snr

    self.opdict['auto_loclevel']=args.auto_loclevel
    self.opdict['loclevel']=args.loclevel
    self.opdict['snr_loclevel']=args.snr_loclevel
    self.opdict['snr_limit']=args.snr_limit
    self.opdict['snr_tr_limit']=args.snr_tr_limit
    self.opdict['sn_time']=args.sn_time
    self.opdict['n_kurt_min']=args.n_kurt_min

    self.opdict['syn_addnoise']=args.syn_addnoise
    self.opdict['syn_snr']=args.syn_snr
    self.opdict['syn_amplitude']=args.syn_amplitude
    self.opdict['syn_datalength']=args.syn_datalength
    self.opdict['syn_samplefreq']=args.syn_samplefreq
    self.opdict['syn_kwidth']=args.syn_kwidth
    self.opdict['syn_otime']=args.syn_otime
    self.opdict['syn_ix']=args.syn_ix
    self.opdict['syn_iy']=args.syn_iy
    self.opdict['syn_iz']=args.syn_iz

    self.opdict['plot_tbefore']=args.plot_tbefore
    self.opdict['plot_tafter']=args.plot_tafter

    self.opdict['xcorr_threshold']=args.xcorr_threshold
    self.opdict['xcorr_before']=args.xcorr_before
    self.opdict['xcorr_after']=args.xcorr_after
    self.opdict['xcorr_corr']=args.xcorr_corr
    self.opdict['xcorr_delay']=args.xcorr_delay

    self.opdict['clus']=args.clus
    self.opdict['nbsta']=args.nbsta

    self.opdict['dd_loc']=args.dd_loc


  def set_test_options(self):
    self.opdict['time']=True
    self.opdict['verbose']=False

    self.opdict['test_datadir']='test_data'
    self.opdict['datadir']='TEST'
    self.opdict['outdir']='TEST'

    self.opdict['net_list']='YA'
    self.opdict['sta_list']="FJS,FLR,FOR,HDL,RVL,SNE,UV01,UV02,UV03,UV04,UV05,UV06,UV07,UV08,UV09,UV10,UV11,UV12,UV13,UV14,UV15"
    self.opdict['comp_list']="HHZ"

    self.opdict['starttime']="2010-10-14T00:14:00.0Z"
    self.opdict['endtime']="2010-10-14T00:18:00.0Z"

    self.opdict['time_grid']='Slow_len.100m.P'
    self.opdict['search_grid']='test_grid.search.hdr'
    self.opdict['stations']='coord_stations_test'

    self.opdict['resample']=False
    self.opdict['fs']=None

    self.opdict['c1']=4.0
    self.opdict['c2']=10.0

    self.opdict['kwin']=4
    self.opdict['krec']=False
    self.opdict['kderiv']=True

    self.opdict['data_length']=600
    self.opdict['data_overlap']=20

    self.opdict['dataglob']='*filt.mseed'
    self.opdict['kurtglob']='*kurt.mseed'
    self.opdict['gradglob']='*grad.mseed'

    self.opdict['load_ttimes_buf']=True

    self.opdict['reloc']=False
    self.opdict['reloc_snr']=12.

    self.opdict['auto_loclevel']=False
    self.opdict['loclevel']=50.0
    self.opdict['snr_limit']=10.0
    self.opdict['snr_tr_limit']=10.0
    self.opdict['sn_time']=10.0
    self.opdict['n_kurt_min']=4

    self.opdict['syn_addnoise']=False

    self.opdict['xcorr_threshold']=0.7
    self.opdict['xcorr_before']=0.5
    self.opdict['xcorr_after']=6.0
    self.opdict['xcorr_corr']='corr'
    self.opdict['xcorr_delay']='delay'

    self.opdict['clus']=0.8
    self.opdict['nbsta']=3

    self.opdict['dd_loc']=True


  def verify_base_path(self):
    """
    Verifies that the base_path is set
    """

    # if the option base_path is not set, then check the environment variable
    # if the environment variable is not set, quit with error message
    if not self.opdict.has_key('base_path'):
      logging.info('No base_path set in options, getting base_path from \
              $WAVELOC_PATH')
      base_path=os.getenv('WAVELOC_PATH')
      if not os.path.isdir(base_path): 
          raise UserWarning('Environment variable WAVELOC_PATH not set \
                  correctly.')
      self.opdict['base_path']=base_path
    
    base_path=self.opdict['base_path']
    lib_path=os.path.join(base_path,'lib')
    if not os.path.isdir(lib_path): 
        raise UserWarning('Directory %s does not exist.'%lib_path)

  def _verify_lib_path(self):
    self.verify_base_path()
    base_path=self.opdict['base_path']
    lib_path=os.path.join(base_path,'lib')
    if not os.path.isdir(lib_path): 
        raise UserWarning('Directory %s does not exist.'%lib_path)

  def _verify_datadir(self):
    self.verify_base_path()
    base_path=self.opdict['base_path']
    if not self.opdict.has_key('datadir'):
        raise UserWarning('datadir option not set')

    datadir=os.path.join(base_path,'data',self.opdict['datadir'])
    if not os.path.isdir(datadir):  
        raise UserWarning('Directory %s does not exist.'%datadir)

  def _verify_outdir(self):
    self.verify_base_path()
    base_path=self.opdict['base_path']
    if not self.opdict.has_key('outdir'):
        raise UserWarning('outdir option not set')

    outdir=os.path.join(base_path,'out',self.opdict['outdir'])
    if not os.path.isdir(outdir):  
      os.makedirs(outdir)
    if not os.path.isdir(os.path.join(outdir,'fig')):  
      os.makedirs(os.path.join(outdir,'fig'))
    if not os.path.isdir(os.path.join(outdir,'grid')):  
      os.makedirs(os.path.join(outdir,'grid'))
    if not os.path.isdir(os.path.join(outdir,'loc')):  
      os.makedirs(os.path.join(outdir,'loc'))
    if not os.path.isdir(os.path.join(outdir,'stack')):  
      os.makedirs(os.path.join(outdir,'stack'))
    if not os.path.isdir(os.path.join(outdir,'time_grids')):  
      os.makedirs(os.path.join(outdir,'time_grids'))
    if self.opdict['reloc'] and not os.path.isdir(os.path.join(outdir,'reloc')):
      os.makedirs(os.path.join(outdir,'reloc'))


  def _verify_net_list(self):
    if not self.opdict.has_key('net_list'):
        raise UserWarning('net_list option not set')
 
  def _verify_sta_list(self):
    if not self.opdict.has_key('sta_list'):
        raise UserWarning('sta_list option not set')

  def _verify_comp_list(self):
    if not self.opdict.has_key('comp_list'):
        raise UserWarning('comp_list option not set')

  def _verify_starttime(self):
    if not self.opdict.has_key('starttime'):
        raise UserWarning('starttime option not set')

  def _verify_endtime(self):
    if not self.opdict.has_key('endtime'):
        raise UserWarning('endtime option not set')

  def _verify_resample(self):
    if not self.opdict.has_key('resample'):
        raise UserWarning('resample option not set')

  def _verify_fs(self):
    self._verify_resample()
    resample=self.opdict['resample']
    if resample:
      if not self.opdict.has_key('fs'):
        raise UserWarning('fs option not set')

  def _verify_c1(self):
    if not self.opdict.has_key('c1'):
        raise UserWarning('c1 option not set')

  def _verify_c2(self):
    if not self.opdict.has_key('c2'):
        raise UserWarning('c2 option not set')

  def _verify_kwin(self):
    if not self.opdict.has_key('kwin'):
        raise UserWarning('kwin option not set')

  def _verify_dataglob(self):
    if not self.opdict.has_key('dataglob'):
        raise UserWarning('dataglob option not set')
    self._verify_datadir()
    base_path=self.opdict['base_path']
    datadir=os.path.join(base_path,'data',self.opdict['datadir'])
    data_names=glob.glob(os.path.join(datadir,self.opdict['dataglob']))
    if len(data_names)==0: 
        raise UserWarning('No data files found : %s'%data_names)

  def _verify_kurtglob(self):
    if not self.opdict.has_key('kurtglob'):
        raise UserWarning('kurtglob option not set')
    self._verify_datadir()
    base_path=self.opdict['base_path']
    datadir=os.path.join(base_path,'data',self.opdict['datadir'])
    kurt_names=glob.glob(os.path.join(datadir,self.opdict['kurtglob']))
    if len(kurt_names)==0: 
        raise UserWarning('No kurtosis files found : %s'%kurt_names)

  def _verify_gradglob(self):
    if not self.opdict.has_key('gradglob'):
        raise UserWarning('gradglob option not set')
    self._verify_datadir()
    base_path=self.opdict['base_path']
    datadir=os.path.join(base_path,'data',self.opdict['datadir'])
    grad_names=glob.glob(os.path.join(datadir,self.opdict['gradglob']))
    if len(grad_names)==0: 
        raise UserWarning('No kurtosis gradient files found : %s'%grad_names)

  def _verify_time_grid(self):
    if not self.opdict.has_key('time_grid'):
        raise UserWarning('time_grid option not set')
    self._verify_lib_path()
    base_path=self.opdict['base_path']
    time_grid=os.path.join(base_path,'lib',self.opdict['time_grid'])
    tg_glob=time_grid+'*'
    tg_files=glob.glob(tg_glob)
    if len(tg_files) == 0 : 
        raise UserWarning('No time grid files found %s'%tg_glob)


  def _verify_data_length(self):
    if not self.opdict.has_key('data_length'):
        raise UserWarning('data_length option not set')

  def _verify_data_overlap(self):
    if not self.opdict.has_key('data_overlap'):
        raise UserWarning('data_overlap option not set')

  def _verify_snr_limit(self):
    if not self.opdict.has_key('snr_limit'):
        raise UserWarning('snr_limit option not set')

  def _verify_snr_tr_limit(self):
    if not self.opdict.has_key('snr_tr_limit'):
        raise UserWarning('snr_tr_limit option not set')

  def _verify_sn_time(self):
    if not self.opdict.has_key('sn_time'):
        raise UserWarning('sn_time option not set')

  def _verify_n_kurt_min(self):
    if not self.opdict.has_key('n_kurt_min'):
        raise UserWarning('n_kurt_min option not set')

  def _verify_stations(self):
    if not self.opdict.has_key('stations'):
        raise UserWarning('stations option not set')
    self._verify_lib_path()
    base_path=self.opdict['base_path']
    stations=os.path.join(base_path,'lib',self.opdict['stations'])
    if not os.path.isfile(stations) : 
        raise UserWarning('Cannot find %s'%stations)

  def _verify_search_grid(self):
    if not self.opdict.has_key('search_grid'):
        raise UserWarning('search_grid option not set')
    self._verify_lib_path()
    base_path=self.opdict['base_path']
    search_grid=os.path.join(base_path,'lib',self.opdict['search_grid'])
    if not os.path.isfile(search_grid) : 
        raise UserWarning('Cannot find %s'%search_grid)

  def _verify_reloc(self):
    if not self.opdict.has_key('reloc'):
        raise UserWarning('reloc option not set')

  def _verify_reloc_snr(self):
    if not self.opdict.has_key('reloc_snr'):
        raise UserWarning('reloc_snr option not set')

  def _verify_auto_loclevel(self):
    if not self.opdict.has_key('auto_loclevel'):
        raise UserWarning('auto_loclevel option not set')

  def _verify_snr_loclevel(self):
    self._verify_auto_loclevel()
    auto_loclevel=self.opdict['auto_loclevel']
    if auto_loclevel:
      if not self.opdict.has_key('snr_loclevel'):
        raise UserWarning('snr_loclevel option not set')

  def _verify_loclevel(self):
    self._verify_auto_loclevel()
    auto_loclevel=self.opdict['auto_loclevel']
    if not auto_loclevel:
      if not self.opdict.has_key('loclevel'):
        raise UserWarning('loclevel option not set')

  def _verify_xcorr_threshold(self):
    if not self.opdict.has_key('xcorr_threshold'):
        raise UserWarning('xcorr_threshold option not set')

  def _verify_xcorr_before(self):
    if not self.opdict.has_key('xcorr_before'):
        raise UserWarning('xcorr_before option not set')

  def _verify_xcorr_after(self):
    if not self.opdict.has_key('xcorr_after'):
        raise UserWarning('xcorr_after option not set')

  def _verify_xcorr_corr(self):
    if not self.opdict.has_key('xcorr_corr'):
        raise UserWarning('xcorr_corr option not set')

  def _verify_xcorr_delay(self):
    if not self.opdict.has_key('xcorr_delay'):
        raise UserWarning('xcorr_delay option not set')

  def _verify_nbsta(self):
    if not self.opdict.has_key('nbsta'):
        raise UserWarning('nbsta option not set')

  def _verify_clus(self):
    if not self.opdict.has_key('clus'):
        raise UserWarning('clus option not set')
  
  def _verify_dd_loc(self):
    if not self.opfict.has_key('dd_loc'):
        raise UserWarning('dd_loc option not set') 

  def _verify_syn_addnoise(self):
    if not self.opdict.has_key('syn_addnoise'):
        raise UserWarning('syn_addnoise option not set')

  def _verify_syn_snr(self):
    self._verify_syn_addnoise()
    syn_addnoise=self.opdict['syn_addnoise']
    if syn_addnoise:
      if not self.opdict.has_key('syn_snr'):
        raise UserWarning('syn_snr option not set')

  def _verify_syn_amplitude(self):
    if not self.opdict.has_key('syn_amplitude'):
        raise UserWarning('syn_amplitude option not set')

  def _verify_syn_datalength(self):
    if not self.opdict.has_key('syn_datalength'):
        raise UserWarning('syn_datalength option not set')

  def _verify_syn_samplefreq(self):
    if not self.opdict.has_key('syn_samplefreq'):
        raise UserWarning('syn_samplefreq option not set')

  def _verify_syn_kwidth(self):
    if not self.opdict.has_key('syn_kwidth'):
        raise UserWarning('syn_kwidth option not set')

  def _verify_syn_otime(self):
    if not self.opdict.has_key('syn_otime'):
        raise UserWarning('syn_otime option not set')

  def _verify_syn_ix(self):
    if not self.opdict.has_key('syn_ix'):
        raise UserWarning('syn_ix option not set')

  def _verify_syn_iy(self):
    if not self.opdict.has_key('syn_iy'):
        raise UserWarning('syn_iy option not set')

  def _verify_syn_iz(self):
    if not self.opdict.has_key('syn_iz'):
        raise UserWarning('syn_iz option not set')

  def _verify_syn_filename(self):
    if not self.opdict.has_key('syn_filename'):
        raise UserWarning('syn_filename option not set')

  def _verify_plot_tbefore(self):
    if not self.opdict.has_key('plot_tbefore'):
        raise UserWarning('plot_tbefore option not set')

  def _verify_plot_tafter(self):
    if not self.opdict.has_key('plot_tafter'):
        raise UserWarning('plot_tafter option not set')


  def verify_SDS_processing_options(self):

    self.verify_base_path()
    self._verify_datadir()

    self._verify_net_list()
    self._verify_sta_list()
    self._verify_comp_list()

    self._verify_starttime()
    self._verify_endtime()

    self._verify_fs()
    self._verify_c1()
    self._verify_c2()
    self._verify_kwin()


  def verify_migration_options(self):

    self.verify_base_path()
    self._verify_lib_path()
    self._verify_datadir()
    self._verify_outdir()

    self._verify_gradglob()
    self._verify_starttime()
    self._verify_endtime()
    self._verify_data_length()
    self._verify_data_overlap()

    self._verify_stations()
    self._verify_search_grid()
    self._verify_time_grid()

    self._verify_reloc()
    if opdict['reloc']:
      self._verify_reloc_snr()

  def verify_location_options(self):

    self.verify_base_path()
    self._verify_lib_path()
    self._verify_datadir()
    self._verify_outdir()

    self._verify_kurtglob()
    self._verify_gradglob()


    self._verify_loclevel()
    self._verify_snr_loclevel()
    self._verify_snr_limit()
    self._verify_snr_tr_limit()
    self._verify_sn_time()
    self._verify_n_kurt_min()

    self._verify_search_grid()
    self._verify_time_grid()

    self._verify_reloc()

  def verify_correlation_options(self):
    
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


  def verify_cluster_options(self):
    
    self.verify_base_path()
    self._verify_lib_path()
    self._verify_datadir()
    self._verify_outdir()

    base_path=self.opdict['base_path']
    locdir=os.path.join(base_path,'out',self.opdict['outdir'],'loc')

    self._verify_dataglob()

    self._verify_stations()
    self._verify_xcorr_corr()
    self._verify_xcorr_delay()

    coeff_file=os.path.join(locdir,self.opdict['xcorr_corr'])
    if not os.path.isfile(coeff_file):
        raise UserWarning('Cannot find %s'%coeff_file)
    
    delay_file=os.path.join(locdir,self.opdict['xcorr_delay'])
    if not os.path.isfile(delay_file):
        raise UserWarning('Cannot find %s'%delay_file)

    self._verify_nbsta()
    self._verify_clus()


  def verify_doublediff_options(self):

    self.verify_base_path()
    self._verify_lib_path()
    self._verify_outdir()
    base_path=self.opdict['base_path']

    self._verify_time_grid()
    self._verify_search_grid()

    locdir=os.path.join(base_path,'out',self.opdict['outdir'],'loc')

    self._verify_stations()
    self._verify_corr()
    self._verify_delay()

    coeff_file=os.path.join(locdir,self.opdict['corr'])
    if not os.path.isfile(coeff_file):
        raise UserWarning('Cannot find %s'%coeff_file)
   
    delay_file=os.path.join(locdir,self.opdict['delay'])
    if not os.path.isfile(delay_file):
        raise UserWarning('Cannot find %s'%delay_file)

    self._verify_nbsta()
    self._verify_clus()

    self._verify_dd_loc()


  def verify_synthetic_options(self):

    self.verify_base_path()
    self._verify_lib_path()
    self._verify_outdir()
    base_path=self.opdict['base_path']

    self._verify_time_grid()

    self._verify_stations()
    self._verify_syn_addnoise()
    self._verify_syn_snr()

    self._verify_syn_amplitude()
    self._verify_syn_datalength()
    self._verify_syn_samplefreq()
    self._verify_syn_kwidth()
    self._verify_syn_otime()
    self._verify_syn_ix()
    self._verify_syn_iy()
    self._verify_syn_iz()
    self._verify_syn_filename()
          


  def verify_plotting_options(self):

    self.verify_base_path()
    self._verify_lib_path()
    self._verify_datadir()
    self._verify_outdir()

    base_path=self.opdict['base_path']
    locdir=os.path.join(base_path,'out',self.opdict['outdir'],'loc')

    locfile=os.path.join(locdir,'locations.dat')
    if not os.path.isfile(locfile): 
        raise UserWarning('Locations file %s does not exist.'%locfile)

    self._verify_dataglob()
    self._verify_kurtglob()
    self._verify_gradglob()

    self._verify_plot_tbefore()
    self._verify_plot_tafter()

    self._verify_search_grid()
    self._verify_time_grid()


    self._verify_stations()
