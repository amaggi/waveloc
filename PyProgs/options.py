import os, glob, argparse

class WavelocOptions(object):

  def __init__(self):

    self.opdict={}

    base_path=os.getenv('WAVELOC_PATH')
    if not os.path.isdir(base_path): raise UserWarning('Environment variable WAVELOC_PATH not set correctly.')
    self.opdict['base_path']=base_path

    # check for existence of lib directory
    lib_path=os.path.join(base_path,'lib')
    if not os.path.isdir(lib_path): raise UserWarning('Directory %s does not exist.'%lib_path)
  
    self.p = argparse.ArgumentParser()

    self.p.add_argument('--time','-t',action='store_true',default=False,help='print timing information to stout')
    self.p.add_argument('--verbose','-v',action='store_true',default=False,help='print debugging information to stout')
  
    self.p.add_argument('--datadir',action='store',help="subdirectory of $WAVELOC_PATH/data")
    self.p.add_argument('--outdir', action='store', help='subdirectory of $WAVELOC_PATH/out for stocking output files')

    self.p.add_argument('--net_list', action='store',help="list of network codes (e.g. \"BE,G\") ")
    self.p.add_argument('--sta_list', action='store',help="list of station names (e.g. \"STA1,STA2\") ")
    self.p.add_argument('--comp_list',action='store',help="list of component names (e.g. \"HHZ,LHZ\") ")

    self.p.add_argument('--resample',action='store_true',default=False, help="resample data")
    self.p.add_argument('--fs',      action='store', type=float,  help="resample frequency")

    self.p.add_argument('--c1',action='store',type=float,  help="low frequency corner of band pass filter ")
    self.p.add_argument('--c2',action='store',type=float,  help="high frequency corner of band pass filter ")
    self.p.add_argument('--kwin',action='store',type=float, help="length of kurtosis window (seconds)")
    self.p.add_argument('--krec',action='store_true',default=False, help="use recursive kurtosis calculation (faster but less precise)")
    self.p.add_argument('--kderiv',action='store_true',default=False, help="use derivative of kurtosis")

    self.p.add_argument('--dataglob',action='store',help="data glob")
    self.p.add_argument('--kurtglob',action='store',help="kurtosis glob")
    self.p.add_argument('--gradglob',action='store',help="gradient glob")

    self.p.add_argument('--starttime',action='store',help="start time for data e.g. 2010-10-14T00:00:00.0Z")
    self.p.add_argument('--endtime',  action='store',help="end time for data e.g. 2010-10-14T10:00:00.0Z")
    self.p.add_argument('--data_length', action='store',type=float,help="length in seconds for data segments to analyse (e.g. 630)")
    self.p.add_argument('--data_overlap',action='store',type=float,help="length in seconds for overlapping data segments (e.g. 30)")

    self.p.add_argument('--stations',action='store',default='channels_HHZ.dat',help='station list (found in $WAVELOC_PATH/lib)')
    self.p.add_argument('--search_grid',action='store',help="search grid e.g. grid.500m.search.hdr (found in $WAVELOC_PATH/lib)")
    self.p.add_argument('--time_grid',  action='store',help="time grid basename e.g. belgium.P (found in $WAVELOC_PATH/lib)")
    self.p.add_argument('--load_ttimes_buf',action='store_true',default=True,help='load pre-calculated travel-times for the search grid from file')

    self.p.add_argument('--reloc', action='store_true', default=False, help='apply to relocated events')
    self.p.add_argument('--loclevel', action='store', default=50,   type=float,help='trigger stack level for locations (e.g. 50) ')
    self.p.add_argument('--snr_limit',action='store', default=10.0, type=float,help="signal_to_noise level for kurtosis acceptance")
    self.p.add_argument('--sn_time',action='store',   default=10.0, type=float,help="time over which to calculate the signal_to_noise ratio for kurtosis acceptance")
    self.p.add_argument('--n_kurt_min',action='store',default=4,    type=int,  help="min number of good kurtosis traces for a location")

   #self.p.add_argument('--2D',action='store_true',default=False,dest='twoD',help='use 2D time grids')


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
    self.opdict['loclevel']=args.loclevel
    self.opdict['snr_limit']=args.snr_limit
    self.opdict['sn_time']=args.sn_time
    self.opdict['n_kurt_min']=args.n_kurt_min

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
    self.opdict['loclevel']=50.0
    self.opdict['snr_limit']=10.0
    self.opdict['sn_time']=10.0
    self.opdict['n_kurt_min']=4

  def verify_SDS_processing_options(self):

    base_path=self.opdict['base_path']

    datadir=os.path.join(base_path,'data',self.opdict['datadir'])
    if self.opdict['datadir']==None:  raise UserWarning('Empty data directory name') 
    if not os.path.isdir(datadir):  raise UserWarning('Data directory %s does not exist'%datadir)

    if self.opdict['net_list']==None:  raise UserWarning('Empty network list') 
    if self.opdict['sta_list']==None:  raise UserWarning('Empty station list') 
    if self.opdict['comp_list']==None: raise UserWarning('Empty component list') 
   
    if self.opdict['starttime']==None: raise UserWarning('Missing start time') 
    if self.opdict['endtime']==None:   raise UserWarning('Missing end time') 
    
    if self.opdict['resample']:
      if self.opdict['fs'] ==None : raise UserWarning('Missing resampling frequency')

    if self.opdict['c1'] ==None :   raise UserWarning('Missing low frequency corner for filtering')
    if self.opdict['c2'] ==None :   raise UserWarning('Missing low frequency corner for filtering')
    if self.opdict['kwin'] ==None : raise UserWarning('Missing kurtosis window length')


  def verify_migration_options(self):

    base_path=self.opdict['base_path']

    if self.opdict['datadir']==None:  raise UserWarning('Empty data directory name') 
    datadir=os.path.join(base_path,'data',self.opdict['datadir'])
    if not os.path.isdir(datadir):  raise UserWarning('Data directory %s does not exist'%datadir)

    if self.opdict['outdir']==None:  raise UserWarning('Empty output directory name') 
    outdir=os.path.join(base_path,'out',self.opdict['outdir'])
    if not os.path.exists(outdir):  
      os.makedirs(outdir)
      os.makedirs(os.path.join(outdir,'stack'))

    if self.opdict['gradglob']==None:  raise UserWarning('Empty gradglob') 
    grad_names=glob.glob(os.path.join(datadir,self.opdict['gradglob']))
    if len(grad_names)==0: raise UserWarning('No kurtosis gradient files found : %s',grad_names)
    
    if self.opdict['starttime']==None: raise UserWarning('Empty start time') 
    if self.opdict['endtime']==None:   raise UserWarning('Empty end time') 
    if self.opdict['data_length']==None:   raise UserWarning('Empty data segment length') 
    if self.opdict['data_overlap']==None:   raise UserWarning('Empty data segment overlap') 

    if self.opdict['stations']==None:   raise UserWarning('Empty stations coordinate file') 
    stations=os.path.join(base_path,'lib',self.opdict['stations'])
    if not os.path.isfile(stations) : raise UserWarning('Cannot find %s'%stations)

    if self.opdict['search_grid']==None:   raise UserWarning('Empty search grid filename') 
    search_grid=os.path.join(base_path,'lib',self.opdict['search_grid'])
    if not os.path.isfile(search_grid) : raise UserWarning('Cannot find %s'%search_grid)

    if self.opdict['time_grid']==None:   raise UserWarning('Empty time grid base filename') 
    time_grid=os.path.join(base_path,'lib',self.opdict['time_grid'])
    tg_glob=time_grid+'*'
    tg_files=glob.glob(tg_glob)
    if len(tg_files) == 0 : raise UserWarning('No time grid files found %s'%tg_glob)

  def verify_location_options(self):

    base_path=self.opdict['base_path']

    if self.opdict['datadir']==None:  raise UserWarning('Empty data directory name') 
    datadir=os.path.join(base_path,'data',self.opdict['datadir'])
    if not os.path.isdir(datadir):  raise UserWarning('Data directory %s does not exist'%datadir)

    if self.opdict['kurtglob']==None:  raise UserWarning('Empty kurtglob') 
    kurt_names=glob.glob(os.path.join(datadir,self.opdict['kurtglob']))
    if len(kurt_names)==0: raise UserWarning('No kurtosis files found : %s',kurt_names)

    if self.opdict['gradglob']==None:  raise UserWarning('Empty gradglob') 
    grad_names=glob.glob(os.path.join(datadir,self.opdict['gradglob']))
    if len(grad_names)==0: raise UserWarning('No kurtosis gradient files found : %s',grad_names)

    out_path=os.path.join(base_path,'out',self.opdict['outdir'])
    if not os.path.isdir(out_path): raise UserWarning('Output directory %s does not exist.') 

    if self.opdict['outdir']==None:  raise UserWarning('Empty output directory name') 
    stackdir=os.path.join(base_path,'out',self.opdict['outdir'],'stack')
    if not os.path.isdir(stackdir): raise UserWarning('Stack directory %s does not exist.  Have you run migration correctly ?') 

    locdir=os.path.join(base_path,'out',self.opdict['outdir'],'loc')
    if not os.path.exists(locdir): os.makedirs(locdir)  

    relocdir=os.path.join(base_path,'out',self.opdict['outdir'],'reloc')
    if self.opdict['reloc'] and not os.path.exists(relocdir): os.makedirs(relocdir)  

    griddir=os.path.join(base_path,'out',self.opdict['outdir'],'grid')
    if not os.path.exists(griddir): os.makedirs(griddir)  

    figdir=os.path.join(base_path,'out',self.opdict['outdir'],'fig')
    if not os.path.exists(figdir): os.makedirs(figdir)  

    if self.opdict['loclevel']==None:   raise UserWarning('Empty location threshold') 
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


