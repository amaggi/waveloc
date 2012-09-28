import os, glob, logging, h5py
import numpy as np

#@profile
def generateSyntheticDirac(opdict,time_grids=None):
    # Creates the synthetic dataset for us to work with

    from NllGridLib import read_stations_file, read_hdr_file
    from migration import migrate_4D_stack, extract_max_values
    from hdf5_grids import get_interpolated_time_grids

    load_time_grids = False
    if time_grids==None : load_time_grids = True

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
    test_stack_file=os.path.join(base_path,'out',opdict['outdir'],'stack','stack_'+s_filename)
    test_info_file=os.path.join(base_path,'out',opdict['outdir'],'grid','%s.info'%s_filename)

    fig_path = os.path.join(base_path,'out',outdir,'fig')

    # data (actual waveforms will not be read - will just use the metadata to set up synthetic problem)
    # if waveforms exist, then read them, else just use all the stations in the stations file
    use_data=True
    try:
      data_dir=os.path.join(base_path,'data',opdict['datadir'])
      data_glob=opdict['gradglob']
      datafile_list=glob.glob(os.path.join(data_dir,data_glob))
    except KeyError:
      logging.info('No data given, so use all stations in the coord_stations file')
      use_data=False
      

    out_dir=os.path.join(base_path,'out',opdict['outdir'])

    # get filenames for time-grids and search grids 
    grid_filename_base=os.path.join(base_path,'lib',opdict['time_grid'])
    search_grid_filename=os.path.join(base_path,'lib',opdict['search_grid'])
    stations_filename=os.path.join(base_path,'lib',opdict['stations'])
    stations=read_stations_file(stations_filename)

    # get parameters for noise etc
    syn_addnoise=opdict['syn_addnoise']

    #################################
    # start setting up synthetic data
    #################################

    grid_info=read_hdr_file(search_grid_filename)

    if load_time_grids:
      time_grids=get_interpolated_time_grids(opdict)

    #################################
    # create synthetic data
    #################################

    # choose hypocenter
    nx=grid_info['nx']
    ny=grid_info['ny']
    nz=grid_info['nz']

    dx=grid_info['dx']
    dy=grid_info['dy']
    dz=grid_info['dz']

    x_orig=grid_info['x_orig']
    y_orig=grid_info['y_orig']
    z_orig=grid_info['z_orig']

    ix=opdict['syn_ix']
    iy=opdict['syn_iy']
    iz=opdict['syn_iz']
    it=int(round(s_t0/s_delta))

    # retrieve travel times for chosen hypocenter 
    ib= ix*ny*nz + iy*nz + iz
    n_buf=nx*ny*nz
    logging.debug('ib for true hypocenter = %d'%ib)
    ttimes={}
    for key in time_grids.keys():
      ttimes[key]=time_grids[key].grid_data[ib]
    #print ttimes
    logging.debug('Travel-times for true hypocenter = %s'%ttimes)

    # construct data with these travel times
    data={}
    for key,delay in ttimes.iteritems():
      if syn_addnoise:
        s_snr=opdict['syn_snr']       
        s=np.random.rand(s_npts)*s_amplitude/s_snr
      else:
        s=np.zeros(s_npts)
      atime=s_t0+delay
      i_atime=np.int(atime/s_delta)
      s[i_atime:i_atime+s_nkwidth]=s_amplitude-np.arange(s_nkwidth)*(s_amplitude/float(s_nkwidth))
      data[key]=s
      

    # DO MIGRATION
    
    logging.info('Doing migration to %s'%test_grid_file)
    f=h5py.File(test_grid_file,'w')
    stack_grid=f.create_dataset('stack_grid',(n_buf,s_npts),'f',chunks=(1,s_npts))
    stack_shift_time = migrate_4D_stack(data,s_delta,time_grids,stack_grid)
    n_buf,nt=stack_grid.shape

    # add useful information to dataset
    for key,value in grid_info.iteritems():
      stack_grid.attrs[key]=value
    stack_grid.attrs['dt']=s_delta
    stack_grid.attrs['start_time']=-stack_shift_time

    # extract max-stack
    logging.info('Extracting max_val etc. to %s'%test_stack_file)
    f_stack = h5py.File(test_stack_file,'w')
    # extract maxima
    extract_max_values(stack_grid,grid_info,f_stack)
    for name in f_stack:
      dset=f_stack[name]
      logging.debug('After extract_max_values : %s %f %f'%(name,np.max(dset),np.sum(dset)))
      dset.attrs['start_time']=-stack_shift_time
      dset.attrs['dt']=s_delta

    # close the stack and grid files 
    f_stack.close()
    f.close()
    logging.info('Saved 4D grid to file %s'%test_grid_file)

    shifted_it=it+int(round(stack_shift_time/s_delta))

    # SETUP information to pass back
    test_info={}
    test_info['dat_file']=test_grid_file
    test_info['stack_file']=test_stack_file
    test_info['grid_shape']=nx,ny,nz,nt
    test_info['grid_spacing']=dx,dy,dz,s_delta
    test_info['grid_orig']=x_orig,y_orig,z_orig
    test_info['true_indexes']=(ix,iy,iz,shifted_it)
    test_info['start_time']=-stack_shift_time

    logging.debug(test_info)
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


