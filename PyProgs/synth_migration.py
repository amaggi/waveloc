import os, glob, logging
import numpy as np

def generateSyntheticDirac(wo):
    # Creates the synthetic dataset for us to work with

    from grids_paths import StationList, ChannelList, QDTimeGrid, migrate_4D_stack

    opdict=wo.opdict

    wo.verify_migration_options()
    wo.verify_synthetic_options()

    #define length and sampling frequency of synthetic data
    s_amplitude   = 1.0 # amplitude on synthetic saveforms
    s_data_length = 20.0 # length of synthetic waveform in seconds
    s_sample_freq = 100.0 # frequency of synthetic waveform in Hz
    
    s_npts=s_data_length*s_sample_freq
    s_delta=1/s_sample_freq
    s_kwidth=0.1
    s_nkwidth=int(round(s_kwidth*s_sample_freq))

    # define origin time
    s_t0 = 6.0


    base_path=opdict['base_path']
    outdir=opdict['outdir']
    test_grid_file=os.path.join(base_path,'out',opdict['outdir'],'test_grid4D_hires.dat')

    fig_path = os.path.join(base_path,'out',outdir,'fig')

    # data (actual waveforms will not be read - will just use the metadata to set up synthetic problem)
    data_dir=os.path.join(base_path,'data',opdict['datadir'])
    out_dir=os.path.join(base_path,'out',opdict['outdir'])
    data_glob=opdict['gradglob']

    # get filenames for time-grids and search grids 
    grid_filename_base=os.path.join(base_path,'lib',opdict['time_grid'])
    search_grid_filename=os.path.join(base_path,'lib',opdict['search_grid'])
    stations_filename=os.path.join(base_path,'lib',opdict['stations'])

    # get parameters for noise etc
    syn_addnoise=opdict['syn_addnoise']

    #################################
    # start setting up synthetic data
    #################################

    sta=StationList()
    sta.read_from_file(stations_filename)

    datafile_list=glob.glob(os.path.join(data_dir,data_glob))

    cha=ChannelList()
    cha.populate_from_station_list_and_data_files(sta,datafile_list)

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

    ix=nx/2
    iy=ny/3
    iz=ny/4
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

    logging.info('Migration outputs : nx,ny,nz,norm_stack_len,stack_shift_time = %d %d %d %d %.3f'%(nx,ny,nz,norm_stack_len,stack_shift_time))
    stack_grid[:,:,:,0:norm_stack_len].tofile(test_grid_file)
    logging.info('Saved 4D grid to file %s'%test_grid_file)

    shifted_it=it+int(round(stack_shift_time/s_delta))

    # SETUP information to pass back
    test_info={}
    test_info['dat_file']=test_grid_file
    test_info['grid_shape']=stack_grid[:,:,:,0:norm_stack_len].shape
    test_info['grid_spacing']=dx,dy,dz,s_delta
    test_info['true_indexes']=(ix,iy,iz,shifted_it)
    test_info['stack_shift_time']=stack_shift_time

    logging.info(test_info)

    return test_info
    
 
