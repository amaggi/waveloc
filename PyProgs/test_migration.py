import os, glob, unittest
import numpy as np
from options import WavelocOptions
from OP_waveforms import Waveform
from migration import do_migration_setup_and_run
from test_processing import waveforms_to_signature

def suite():
  suite = unittest.TestSuite()
  suite.addTest(MigrationTests('test_migration'))
  suite.addTest(MigrationTests('test_migration_fullRes'))
  suite.addTest(MigrationTests('test_migration_synthetic'))
  return suite

    
class MigrationTests(unittest.TestCase):

  def setUp(self):

    self.wo=WavelocOptions()
    self.wo.set_test_options()
    self.wo.verify_migration_options()

  def test_migration_synthetic(self):
    from grids_paths import StationList, ChannelList, QDTimeGrid, migrate_4D_stack
    from integrate4D import compute_expected_coordinates4D, compute_integral4D, compute_expected_coordinates1D 
    from plot_mpl import plot_locations_static_matplotlib, plot_test
    
    #define length and sampling frequency of synthetic data
    s_data_length = 20.0 # seconds
    s_sample_freq = 10.0 # Hz
    s_npts=s_data_length*s_sample_freq
    s_delta=1/s_sample_freq

    # define origin time
    s_t0 = 6.0

    # set some options for output - may not be needed
    self.wo.opdict['outdir'] = 'TEST_Dirac'
    #self.wo.opdict['load_ttimes_buf'] = False
    self.wo.opdict['load_ttimes_buf'] = True

    # verify consistency after changing default options
    self.wo.verify_migration_options()
    self.wo.verify_location_options()
   
    # DO NOT SET ANY MORE opdict options after this line
    opdict=self.wo.opdict

    base_path=opdict['base_path']
    outdir=opdict['outdir']

    fig_path = os.path.join(base_path,'out',outdir,'fig')

    # data (actual waveforms will not be read - will just use the metadata to set up synthetic problem)
    data_dir=os.path.join(base_path,'data',opdict['datadir'])
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
    time_grid.populate_from_time_grids(grid_filename_base,cha,load_ttimes_buf)

    
    #################################
    # create synthetic data
    #################################

    # choose hypocenter
    nx=time_grid.nx
    ny=time_grid.ny
    nz=time_grid.nz

    ix=nx/2
    iy=ny/3
    iz=ny/4
    it=s_t0/s_delta

    logging.debug('True ix, iy, iz, it = %d %d %d %s'%(ix,iy,iz,it))

    # retrieve travel times for chosen hypocenter 
    ib= ix*ny*nz + iy*nz + iz
    ttimes=time_grid.buf[ib]
    logging.debug('ib for true hypocenter = %d'%ib)
    #print ttimes
    logging.debug('Travel-times for true hypocenter = %s'%ttimes)

    # construct data with these travel times
    integer_data={}
    for key,delay in ttimes.iteritems():
      #s=np.random.randn(s_npts)
      s=np.zeros(s_npts)
      atime=s_t0+delay
      i_atime=np.int(atime/s_delta)
      #s[i_atime-2:i_atime+2]=5.0
      s[i_atime]=50.0
      integer_data[key]=s
      

    # DO MIGRATION
    (n_buf, norm_stack_len, stack_shift_time, stack_grid) = migrate_4D_stack(integer_data,s_delta,search_grid_filename,time_grid)

    stack_grid.tofile('broken_grid.dat')

    

    x0=np.arange(nx)*time_grid.dx
    x1=np.arange(ny)*time_grid.dy
    x2=np.arange(nz)*time_grid.dz
    xt=np.arange(norm_stack_len)*s_delta

    exp_x0,exp_x1,exp_x2,exp_xt,cov_matrix,prob_dict = compute_expected_coordinates4D(stack_grid[0:nx,0:ny,0:nz,0:norm_stack_len],x0,x1,x2,xt,return_2Dgrids=True)

    i_xt = int(round(exp_xt/s_delta)) 
    ib=np.argmax(stack_grid[:,:,:,i_xt])
    (ix,iy,iz)=np.unravel_index(ib,(nx,ny,nz))
    logging.debug('Absolute maximum of stack (%.3f) at ix=%d, iy=%d, iz=%d, it=%d, ib=%d'%(stack_grid[ix,iy,iz,i_xt],ix,iy,iz,i_xt,ib))
    stack_x=stack_grid[:,iy,iz,i_xt]
    stack_y=stack_grid[ix,:,iz,i_xt]
    stack_z=stack_grid[ix,iy,:,i_xt]

    exp_x, var_x = compute_expected_coordinates1D(stack_x,x0)
    exp_y, var_y = compute_expected_coordinates1D(stack_y,x1)
    exp_z, var_z = compute_expected_coordinates1D(stack_z,x2)

    sigma_x=np.sqrt(var_x)
    sigma_y=np.sqrt(var_y)
    sigma_z=np.sqrt(var_z)
    

    sigma_x0 = np.sqrt(cov_matrix[0,0])
    sigma_x1 = np.sqrt(cov_matrix[1,1])
    sigma_x2 = np.sqrt(cov_matrix[2,2])
    sigma_xt = np.sqrt(cov_matrix[3,3])

    logging.info("Time %s s pm %.2fs, x=%.4f pm %.4f, y=%.4f pm %.4f, z=%.4f pm %.4f"%(exp_xt,sigma_xt,exp_x,sigma_x, exp_y,sigma_y,exp_z,sigma_z))

    logging.info("Max = %.2f, Time %s s pm %.2fs, x=%.4f pm %.4f, y=%.4f pm %.4f, z=%.4f pm %.4f"%(stack_grid[0:nx,0:ny,0:nz,0:norm_stack_len].max(),exp_xt,sigma_xt, exp_x0,sigma_x0,exp_x1,sigma_x1,exp_x2,sigma_x2))

    fig_name=os.path.join(fig_path,'fig_synt_st_mpl')
    #plot_locations_static_matplotlib(prob_dict,[x0,x1,x2,xt],fig_name)
    plot_test((stack_x, stack_y, stack_z, prob_dict['prob_x3']),(x0,x1,x2,xt),fig_name)




  @unittest.skip('Not running small test')
  def test_migration(self):

    base_path=self.wo.opdict['base_path']
    test_datadir=self.wo.opdict['test_datadir']
    outdir=self.wo.opdict['outdir']

    expected_signature_filename = os.path.join(base_path,test_datadir,'test_stack_signature.dat')
    expected_signature_file = open(expected_signature_filename,'r') 
    expected_lines=expected_signature_file.readlines()

    do_migration_setup_and_run(self.wo.opdict)

    waveforms_to_signature(base_path,os.path.join('out',outdir,'stack'),'stack*mseed','stack_signature.dat')
    signature_filename=os.path.join(base_path,'out',outdir,'stack','stack_signature.dat')
    signature_file = open(signature_filename,'r') 
    lines=signature_file.readlines()

    self.assertSequenceEqual(lines,expected_lines)

  @unittest.skip('Not running full resolution test')
  def test_migration_fullRes(self):

    self.wo.opdict['search_grid'] = 'grid.Taisne.search.hdr'
    self.wo.opdict['outdir'] = 'TEST_fullRes'
    self.wo.opdict['load_ttimes_buf'] = False
    self.wo.opdict['data_length'] = 300

    base_path=self.wo.opdict['base_path']
    test_datadir=self.wo.opdict['test_datadir']
    outdir=self.wo.opdict['outdir']
    self.wo.verify_migration_options()

    expected_signature_filename = os.path.join(base_path,test_datadir,'test_stack_signature.dat')
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
  logging.basicConfig(level=logging.DEBUG, format='%(levelname)s : %(asctime)s : %(message)s')
 
  unittest.TextTestRunner(verbosity=2).run(suite())
 
