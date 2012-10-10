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

    # location
    self.opdict['reloc']=False
    self.opdict['auto_loclevel']=False
    self.opdict['loclevel']=50.
    self.opdict['snr_loclevel']=10.
    self.opdict['snr_limit']=10.
    self.opdict['sn_time']=10.
    self.opdict['n_kurt_min']=4

    # synthetic
    self.opdict['syn_addnoise']=False
    self.opdict['syn_amplitude']=1.
    self.opdict['syn_kwidth']=0.1

    # cross-correlation
    self.opdict['threshold']=0.7
    self.opdict['before']=0.5
    self.opdict['after']=6.0

    # clustering
    self.opdict['nbsta']=3
    self.opdict['clus']=0.8

  
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
            default=self.opdict['reloc'], help='apply to relocated events')
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

    self.p.add_argument('--threshold',action='store',
            default=self.opdict['threshold'], type=float, 
            help="correlation value over which the correlation is computed \
                    again in the Fourier domain")
    self.p.add_argument('--before',action='store',
            default=self.opdict['before'], type=float, help=
            "cross-correlation window: time interval before the origin time")
    self.p.add_argument('--after',action='store', default=self.opdict['after'],
            type=float, help="cross-correlation window: time interval after \
                    the origin time")
    self.p.add_argument('--corr', action='store', 
            help="name of the file containing all correlation values")
    self.p.add_argument('--delay', action='store', 
            help="name of the file containing all time delays")

    self.p.add_argument('--nbsta',action='store', default=self.opdict['nbsta'],
            type=int, help="number of stations over which an event pair is \
            considered provided that its correlation coefficient is greater \
            than a given threshold")
    self.p.add_argument('--clus',action='store', default=self.opdict['clus'],
            type=float, 
            help="correlation value over which an event pair is considered")


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
    self.opdict['auto_loclevel']=args.auto_loclevel
    self.opdict['loclevel']=args.loclevel
    self.opdict['snr_loclevel']=args.snr_loclevel
    self.opdict['snr_limit']=args.snr_limit
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

    self.opdict['threshold']=args.threshold
    self.opdict['before']=args.before
    self.opdict['after']=args.after
    self.opdict['corr']=args.corr
    self.opdict['delay']=args.delay

    self.opdict['clus']=args.clus
    self.opdict['nbsta']=args.nbsta

    self.opdict['new_file']=args.new_file
    self.opdict['refine']=args.refine


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
    self.opdict['auto_loclevel']=False
    self.opdict['loclevel']=50.0
    self.opdict['snr_limit']=10.0
    self.opdict['sn_time']=10.0
    self.opdict['n_kurt_min']=4

    self.opdict['syn_addnoise']=False

    self.opdict['threshold']=0.7
    self.opdict['before']=0.5
    self.opdict['after']=6.0
    self.opdict['corr']='corr'
    self.opdict['delay']='delay'

    self.opdict['clus']=0.8
    self.opdict['nbsta']=3


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

  def _verify_data_length(self):
    if not self.opdict.has_key('data_length'):
        raise UserWarning('data_length option not set')

  def _verify_data_overlap(self):
    if not self.opdict.has_key('data_overlap'):
        raise UserWarning('data_overlap option not set')


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

    base_path=self.opdict['base_path']

    self._verify_gradglob()
    self._verify_starttime()
    self._verify_endtime()
    self._verify_data_length()
    self._verify_data_overlap()

    if not self.opdict.has_key('stations') or self.opdict['stations']==None:   raise UserWarning('Empty stations coordinate file') 
    stations=os.path.join(base_path,'lib',self.opdict['stations'])
    if not os.path.isfile(stations) : raise UserWarning('Cannot find %s'%stations)

    if not self.opdict.has_key('search_grid') or self.opdict['search_grid']==None:   raise UserWarning('Empty search grid filename') 
    search_grid=os.path.join(base_path,'lib',self.opdict['search_grid'])
    if not os.path.isfile(search_grid) : raise UserWarning('Cannot find %s'%search_grid)

    if self.opdict['time_grid']==None:   raise UserWarning('Empty time grid base filename') 
    time_grid=os.path.join(base_path,'lib',self.opdict['time_grid'])
    tg_glob=time_grid+'*'
    tg_files=glob.glob(tg_glob)
    if len(tg_files) == 0 : raise UserWarning('No time grid files found %s'%tg_glob)

  def verify_location_options(self):

    self.verify_base_path()
    self._verify_lib_path()
    self._verify_datadir()
    self._verify_outdir()

    base_path=self.opdict['base_path']
    datadir=os.path.join(base_path,'data',self.opdict['datadir'])
    out_path=os.path.join(base_path,'out',self.opdict['outdir'])

    self._verify_kurtglob()
    self._verify_gradglob()


    if self.opdict['auto_loclevel']: 
      if self.opdict['snr_loclevel']==None :   raise UserWarning('Empty snr for automatic location threshold') 
    else :
      if self.opdict['loclevel']==None :   raise UserWarning('Empty location threshold') 
    if self.opdict['snr_limit']==None:   raise UserWarning('Empty threshold for signal to noise ratio') 
    if self.opdict['sn_time']==None:   raise UserWarning('Empty time span for signal to noise ratio computation') 
    if self.opdict['n_kurt_min']==None:   raise UserWarning('Empty minimum number of good kurtosis for location') 

    if self.opdict['search_grid']==None:   raise UserWarning('Empty search grid filename') 
    search_grid=os.path.join(base_path,'lib',self.opdict['search_grid'])
    if not os.path.isfile(search_grid) : raise UserWarning('Cannot find %s'%search_grid)

    if self.opdict['time_grid']==None:   raise UserWarning('Empty time grid base filename') 
    time_grid=os.path.join(base_path,'lib',self.opdict['time_grid'])
    tg_glob=time_grid+'*'
    tg_files=glob.glob(tg_glob)
    if len(tg_files) == 0 : raise UserWarning('No time grid files found %s'%tg_glob)

  def verify_correlation_options(self):
    
    self.verify_base_path()
    self._verify_lib_path()
    self._verify_datadir()
    self._verify_outdir()

    base_path=self.opdict['base_path']
    datadir=os.path.join(base_path,'data',self.opdict['datadir'])


    self._verify_dataglob()

    if self.opdict['threshold']==None:  raise UserWarning('Empty correlation threshold for refinement in Fourier domain')
    if self.opdict['before']==None:  raise UserWarning('Empty lower limit for correlation time window')
    if self.opdict['after']==None:  raise UserWarning('Empty upper limit for correlation time window')

    if self.opdict['corr']==None:  raise UserWarning('Empty correlation file name')
    coeff_file=os.path.join(locdir,self.opdict['corr'])
    
    if self.opdict['delay']==None:  raise UserWarning('Empty time delays file name')
    delay_file=os.path.join(locdir,self.opdict['delay'])

  def verify_cluster_options(self):
    
    self.verify_base_path()
    self._verify_lib_path()
    self._verify_datadir()
    self._verify_outdir()

    base_path=self.opdict['base_path']
    datadir=os.path.join(base_path,'data',self.opdict['datadir'])

    self._verify_dataglob()


    if self.opdict['stations']==None:   raise UserWarning('Empty stations coordinate file') 
    stations=os.path.join(base_path,'lib',self.opdict['stations'])
    if not os.path.isfile(stations) : raise UserWarning('Cannot find %s'%stations)

    if self.opdict['corr']==None:  raise UserWarning('Empty correlation file')
    coeff_file=os.path.join(locdir,self.opdict['corr'])
    if not os.path.isfile(coeff_file):  raise UserWarning('Cannot find %s'%coeff_file)
    
    if self.opdict['delay']==None:  raise UserWarning('Empty time delays file')
    delay_file=os.path.join(locdir,self.opdict['delay'])
    if not os.path.isfile(delay_file):  raise UserWarning('Cannot find %s'%delay_file)

    if self.opdict['nbsta']==None:  raise UserWarning('Empty minimum number of stations')
    if self.opdict['clus']==None:  raise UserWarning('Empty correlation threshold for clustering')


  def verify_synthetic_options(self):

    self.verify_base_path()
    self._verify_lib_path()
    self._verify_outdir()
    base_path=self.opdict['base_path']

    if not self.opdict.has_key('time_grid') or self.opdict['time_grid']==None:   raise UserWarning('Empty time grid base filename') 
    time_grid=os.path.join(base_path,'lib',self.opdict['time_grid'])
    tg_glob=time_grid+'*'
    tg_files=glob.glob(tg_glob)
    if len(tg_files) == 0 : raise UserWarning('No time grid files found %s'%tg_glob)

    if self.opdict['stations']==None:   raise UserWarning('Empty stations coordinate file') 
    stations=os.path.join(base_path,'lib',self.opdict['stations'])
    if not os.path.isfile(stations) : raise UserWarning('Cannot find %s'%stations)

    if self.opdict['syn_addnoise'] :
      if self.opdict['syn_snr']==None:	raise UserWarning('No SNR set for synthetic test')
          
    if self.opdict['syn_amplitude']==None:	raise UserWarning('No synthetic amplitudue set')
    if self.opdict['syn_datalength']==None:	raise UserWarning('No synthetic datalength set')  
    if self.opdict['syn_samplefreq']==None:	raise UserWarning('No synthetic samplefreq set')  
    if self.opdict['syn_kwidth']==None:	raise UserWarning('No synthetic kwidth set')  
    if self.opdict['syn_otime']==None:	raise UserWarning('No synthetic otime set')  
    if self.opdict['syn_ix']==None:	raise UserWarning('No synthetic ix set')  
    if self.opdict['syn_iy']==None:	raise UserWarning('No synthetic iy set')  
    if self.opdict['syn_iz']==None:	raise UserWarning('No synthetic iz set')  
    if self.opdict['syn_filename']==None:	raise UserWarning('No filename set for synthetic grid')  



  def verify_plotting_options(self):

    self.verify_base_path()
    self._verify_lib_path()
    self._verify_datadir()
    self._verify_outdir()

    base_path=self.opdict['base_path']


    locfile=os.path.join(base_path,'out',self.opdict['outdir'],'loc','locations.dat')
    if not os.path.isfile(locfile): raise UserWarning('Locations file %s does not exist.'%locfile)

    self._verify_dataglob()
    self._verify_kurtglob()
    self._verify_gradglob()

    if not self.opdict.has_key('plot_tbefore') : raise UserWarning('Missing start time for plots (plot_tbefore)')
    if not self.opdict.has_key('plot_tafter') : raise UserWarning('Missing end time for plots (plot_tafter)')

    if not self.opdict.has_key('search_grid') or self.opdict['search_grid']==None:   raise UserWarning('Empty search grid filename') 
    search_grid=os.path.join(base_path,'lib',self.opdict['search_grid'])
    if not os.path.isfile(search_grid) : raise UserWarning('Cannot find %s'%search_grid)


    if not self.opdict.has_key('time_grid') or self.opdict['time_grid']==None:   raise UserWarning('Empty time grid base filename') 
    time_grid=os.path.join(base_path,'lib',self.opdict['time_grid'])
    tg_glob=time_grid+'*'
    tg_files=glob.glob(tg_glob)
    if len(tg_files) == 0 : raise UserWarning('No time grid files found %s'%tg_glob)


    if not self.opdict.has_key('stations') or self.opdict['stations']==None:   raise UserWarning('Empty stations coordinate file') 
    stations=os.path.join(base_path,'lib',self.opdict['stations'])
    if not os.path.isfile(stations) : raise UserWarning('Cannot find %s'%stations)
