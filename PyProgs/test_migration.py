import os, glob, unittest
import numpy as np
from options import WavelocOptions
from OP_waveforms import Waveform
from migration import do_migration_setup_and_run
from test_processing import waveforms_to_signature
from integrate4D import * 

def generateSyntheticDirac(wo,add_noise=False):
    # Creates the synthetic dataset for us to work with

    from grids_paths import StationList, ChannelList, QDTimeGrid, migrate_4D_stack

    wo.verify_migration_options()

    #define length and sampling frequency of synthetic data
    s_amplitude   = 1.0 # amplitude on synthetic saveforms
    s_data_length = 20.0 # length of synthetic waveform in seconds
    s_sample_freq = 100.0 # frequency of synthetic waveform in Hz
    s_snr=3.0
    s_npts=s_data_length*s_sample_freq
    s_delta=1/s_sample_freq
    s_kwidth=0.1
    s_nkwidth=int(round(s_kwidth*s_sample_freq))

    # define origin time
    s_t0 = 6.0


    # DO NOT SET ANY MORE opdict options after this line
    opdict=wo.opdict

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
      if add_noise:
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
    
 

def suite():
  suite = unittest.TestSuite()
  suite.addTest(SyntheticMigrationTests('test_dirac_migration'))
  suite.addTest(MigrationTests('test_migration'))
  suite.addTest(MigrationTests('test_migration_fullRes'))
  return suite

class SyntheticMigrationTests(unittest.TestCase):

  def test_dirac_migration(self):
    from locations_trigger import trigger_locations_inner
    from plot_mpl import plot_probloc_mpl

    wo=WavelocOptions()
    wo.set_test_options()

    #wo.opdict['outdir'] = 'TEST_DiracNoisy'
    wo.opdict['outdir'] = 'TEST_Dirac'
    wo.opdict['search_grid']='grid.Taisne.search.hdr'
    wo.opdict['loclevel'] = 10
    wo.opdict['load_ttimes_buf'] = True # Optimized in time, but you must be usre you're reading the right grid for the test

    wo.verify_migration_options()
    wo.verify_location_options()

    # generate the test case and retrieve necessary information
    logging.info('Running synthetic test case generation...')
    #test_info=generateSyntheticDirac(wo,add_noise=True)
    test_info=generateSyntheticDirac(wo)

    # retrieve info
    dat_file=test_info['dat_file']
    nx,ny,nz,nt=test_info['grid_shape']
    dx,dy,dz,dt=test_info['grid_spacing']
    ix_true,iy_true,iz_true,it_true=test_info['true_indexes']
    stack_shift_time=test_info['stack_shift_time']

    # plot base filename
    base_path=wo.opdict['base_path']
    outdir=wo.opdict['outdir']
    plot_base_filename=os.path.join(base_path,'out',outdir,'fig','fig_synt_st_mpl')

    # loclevel for triggers
    loclevel=wo.opdict['loclevel']

    # set up x, y, z, t arrays
    x=np.arange(nx)*dx
    y=np.arange(ny)*dy
    z=np.arange(nz)*dz
    t=np.arange(nt)*dt-stack_shift_time

    # load grid
    stack_grid=np.fromfile(dat_file).reshape(nx,ny,nz,nt)

    # normalize grid for first probability density calculation
    stack_grid_int=compute_integral4D(stack_grid,x,y,z,t)
    stack_grid_norm=stack_grid / stack_grid_int

 
    # integrate normalized grid over all space dimensions to get marginal over time
    prob_t = si.trapz(si.trapz(si.trapz(stack_grid_norm,x=x,axis=0),x=y,axis=0),x=z,axis=0)
    exp_t = si.trapz(t*prob_t,x=t,axis=0)
    var_t = si.trapz((t-exp_t)*(t-exp_t)*prob_t,x=t,axis=0)
    logging.debug('var_t = %.3f'%var_t)
    sigma_t = np.sqrt(var_t)
    logging.debug('sigma_t = %.3f'%sigma_t)
    it_exp=int(round(exp_t/dt))
    nt_sigma=int(round(sigma_t/dt))
    it_left=it_exp-nt_sigma
    it_right=it_exp+nt_sigma
    t_slice=t[it_left:it_right]
    
    # simulate locations trigger
    max_val=stack_grid.max(0).max(0).max(0)
    max_x=stack_grid.max(2).max(1).argmax(0)*dx
    max_y=stack_grid.max(2).max(0).argmax(0)*dy
    max_z=stack_grid.max(1).max(0).argmax(0)*dz

    
    locs=trigger_locations_inner(max_val,max_x,max_y,max_z,loclevel,loclevel,dt)
    
    #print locs
    # This is a dirac test, so only have one element in locs
    trig_loc=locs[0]
    trig_max,trig_t,trig_sigma_t_left,trig_sigma_t_right,trig_x,trig_sigma_x,trig_y,trig_sigma_y,trig_z,trig_sigma_z = trig_loc

    logging.info("TRIGGER : Max = %.2f, Time %s s pm %.2fs, x=%.4f pm %.4f, y=%.4f pm %.4f, z=%.4f pm %.4f"%(trig_max,trig_t-stack_shift_time,max(trig_sigma_t_left,trig_sigma_t_right), trig_x, trig_sigma_x,trig_y,trig_sigma_y,trig_z,trig_sigma_z))
  
    # TODO - send to a plotter

    exp_x,exp_y,exp_z,exp_t,cov_matrix,prob_dict = compute_expected_coordinates4D(stack_grid[:,:,:,it_exp-nt_sigma:it_exp+nt_sigma],x,y,z,t_slice,return_2Dgrids=True)
    sigma_x=np.sqrt(cov_matrix[0,0])
    sigma_y=np.sqrt(cov_matrix[1,1])
    sigma_z=np.sqrt(cov_matrix[2,2])
    sigma_t=np.sqrt(cov_matrix[3,3])

    logging.info("PROB DENSITY : Time %s s pm %.2fs, x=%.4f pm %.4f, y=%.4f pm %.4f, z=%.4f pm %.4f"%(exp_t,sigma_t, exp_x, sigma_x,exp_y,sigma_y,exp_z,sigma_z))

    plot_probloc_mpl(prob_dict,[x,y,z,t_slice],plot_base_filename)


    self.assertTrue(True)
   
@unittest.skip('Not running real data migration tests')
class MigrationTests(unittest.TestCase):

  def setUp(self):

    self.wo=WavelocOptions()
    self.wo.set_test_options()
    self.wo.verify_migration_options()

 

#  @unittest.skip('Not running small test')
  def test_migration(self):

    self.wo.opdict['load_ttimes_buf'] = False
    self.wo.opdict['data_length'] = 300

    base_path=self.wo.opdict['base_path']
    test_datadir=self.wo.opdict['test_datadir']
    outdir=self.wo.opdict['outdir']

    expected_signature_filename = os.path.join(base_path,test_datadir,'TEST_stack_signature.dat')
    expected_signature_file = open(expected_signature_filename,'r') 
    expected_lines=expected_signature_file.readlines()

    do_migration_setup_and_run(self.wo.opdict)

    waveforms_to_signature(base_path,os.path.join('out',outdir,'stack'),'stack*mseed','stack_signature.dat')
    signature_filename=os.path.join(base_path,'out',outdir,'stack','stack_signature.dat')
    signature_file = open(signature_filename,'r') 
    lines=signature_file.readlines()

    self.assertSequenceEqual(lines,expected_lines)

#  @unittest.skip('Not running full resolution test')
  def test_migration_fullRes(self):

    self.wo.opdict['search_grid'] = 'grid.Taisne.search.hdr'
    self.wo.opdict['outdir'] = 'TEST_fullRes'
    self.wo.opdict['load_ttimes_buf'] = True
    self.wo.opdict['data_length'] = 100
    self.wo.verify_migration_options()

    base_path=self.wo.opdict['base_path']
    test_datadir=self.wo.opdict['test_datadir']
    outdir=self.wo.opdict['outdir']

    expected_signature_filename = os.path.join(base_path,test_datadir,'TEST_fullRes_stack_signature.dat')
    expected_signature_file = open(expected_signature_filename,'r') 
    expected_lines=expected_signature_file.readlines()

    do_migration_setup_and_run(self.wo.opdict)

    waveforms_to_signature(base_path,os.path.join('out',outdir,'stack'),'stack*mseed','stack_signature.dat')
    signature_filename=os.path.join(base_path,'out',outdir,'stack','stack_signature.dat')
    signature_file = open(signature_filename,'r') 
    lines=signature_file.readlines()

    self.assertSequenceEqual(lines,expected_lines)



if __name__ == '__main__':

  import logging
  logging.basicConfig(level=logging.INFO, format='%(levelname)s : %(asctime)s : %(message)s')
 
  unittest.TextTestRunner(verbosity=2).run(suite())
 
