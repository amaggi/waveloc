import os, glob, logging
import numpy as np

#@profile
def generateSyntheticDirac(opdict,time_grid=None):
    # Creates the synthetic dataset for us to work with

    from grids_paths import StationList, ChannelList, QDTimeGrid
    from hdf5_grids import migrate_4D_stack

    load_time_grids = False
    if time_grid==None : load_time_grids = True

    #define length and sampling frequency of synthetic data
    s_amplitude   = opdict['syn_amplitude']
    s_data_length = opdict['syn_datalength']
    s_sample_freq = opdict['syn_samplefreq']
    s_filename    = opdict['syn_filename']

    
    s_npts=int(s_data_length*s_sample_freq)
    s_delta=1/s_sample_freq
    s_kwidth=opdict['syn_kwidth']
    s_nkwidth=int(round(s_kwidth*s_sample_freq))

    # define origin time
    s_t0 = opdict['syn_otime']


    base_path=opdict['base_path']
    outdir=opdict['outdir']
    test_grid_file=os.path.join(base_path,'out',opdict['outdir'],'grid',s_filename)
    test_info_file=os.path.join(base_path,'out',opdict['outdir'],'grid','%s.info'%s_filename)

    fig_path = os.path.join(base_path,'out',outdir,'fig')

    # data (actual waveforms will not be read - will just use the metadata to set up synthetic problem)
    # if waveforms exist, then read them, else just use all the stations in the stations file
    use_data=True
    try:
      data_dir=os.path.join(base_path,'data',opdict['datadir'])
      data_glob=opdict['gradglob']
    except KeyError:
      logging.info('No data given, so use all stations in the coord_stations file')
      use_data=False
      

    out_dir=os.path.join(base_path,'out',opdict['outdir'])

    # get filenames for time-grids and search grids 
    grid_filename_base=os.path.join(base_path,'lib',opdict['time_grid'])
    search_grid_filename=os.path.join(base_path,'lib',opdict['search_grid'])
    stations_filename=os.path.join(base_path,'lib',opdict['stations'])

    # get parameters for noise etc
    syn_addnoise=opdict['syn_addnoise']

    #################################
    # start setting up synthetic data
    #################################

    if load_time_grids:
      sta=StationList()
      sta.read_from_file(stations_filename)

      cha=ChannelList()
      if use_data : 
        datafile_list=glob.glob(os.path.join(data_dir,data_glob))
        cha.populate_from_station_list_and_data_files(sta,datafile_list)
      else : 
        cha.populate_from_station_list(sta,comp_string=["HHZ"])

      time_grid=QDTimeGrid()
      time_grid.read_NLL_hdr_file(search_grid_filename)
      load_ttimes_buf=opdict['load_ttimes_buf']
      time_grid.populate_from_time_grids(grid_filename_base,cha,out_dir,load_ttimes_buf)

    
    #################################
    # create synthetic data
    #################################

    # choose hypocenter
    nx=time_grid.nx
    ny=time_grid.ny
    nz=time_grid.nz

    dx=time_grid.dx
    dy=time_grid.dy
    dz=time_grid.dz

    x_orig=time_grid.x_orig
    y_orig=time_grid.y_orig
    z_orig=time_grid.z_orig

    ix=opdict['syn_ix']
    iy=opdict['syn_iy']
    iz=opdict['syn_iz']
    it=int(round(s_t0/s_delta))

    # retrieve travel times for chosen hypocenter 
    ib= ix*ny*nz + iy*nz + iz
    ttimes=time_grid.buf[ib]
    logging.debug('ib for true hypocenter = %d'%ib)
    #print ttimes
    logging.debug('Travel-times for true hypocenter = %s'%ttimes)

    # construct data with these travel times
    integer_data={}
    for key,delay in ttimes.iteritems():
      if syn_addnoise:
        s_snr=opdict['syn_snr']       
        s=np.random.rand(s_npts)*s_amplitude/s_snr
      else:
        s=np.zeros(s_npts)
      atime=s_t0+delay
      i_atime=np.int(atime/s_delta)
      s[i_atime:i_atime+s_nkwidth]=s_amplitude-np.arange(s_nkwidth)*(s_amplitude/float(s_nkwidth))
      integer_data[key]=s
      

    # DO MIGRATION
    (n_buf, norm_stack_len, stack_shift_time, stack_grid) = migrate_4D_stack(integer_data,s_delta,search_grid_filename,time_grid)
    del integer_data
    del time_grid

    logging.info('Migration outputs : nx,ny,nz,norm_stack_len,stack_shift_time = %d %d %d %d %.3f'%(nx,ny,nz,norm_stack_len,stack_shift_time))
    stack_grid[:,:,:,0:norm_stack_len].tofile(test_grid_file)
    logging.info('Saved 4D grid to file %s'%test_grid_file)

    shifted_it=it+int(round(stack_shift_time/s_delta))

    # SETUP information to pass back
    test_info={}
    test_info['dat_file']=test_grid_file
    test_info['grid_shape']=stack_grid[:,:,:,0:norm_stack_len].shape
    test_info['grid_spacing']=dx,dy,dz,s_delta
    test_info['grid_orig']=x_orig,y_orig,z_orig
    test_info['true_indexes']=(ix,iy,iz,shifted_it)
    test_info['stack_shift_time']=stack_shift_time

    logging.info(test_info)
    f=open(test_info_file,'w')
    f.write(str(test_info))

    return test_info
    
 
if __name__ == '__main__':

  from options import WavelocOptions
  logging.basicConfig(level=logging.INFO, format='%(levelname)s : %(asctime)s : %(message)s')

  wo = WavelocOptions()
  args=wo.p.parse_args()

  wo.set_all_arguments(args)
  wo.verify_migration_options()
  wo.verify_synthetic_options()

  generateSynteticDirac(wo.opdict)


