import os, glob, unittest
import numpy as np
from options import WavelocOptions
from OP_waveforms import Waveform
from migration import do_migration_setup_and_run
from test_processing import waveforms_to_signature
from integrate4D import * 
from synth_migration import generateSyntheticDirac



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

    wo.opdict['outdir'] = 'TEST_Dirac'
    wo.opdict['search_grid']='grid.Taisne.search.hdr'
    wo.opdict['loclevel'] = 10
    wo.opdict['load_ttimes_buf'] = True # Optimized in time, but you must be usre you're reading the right grid for the test
    wo.opdict['syn_addnoise']=False
    wo.opdict['syn_amplitude']=1.0
    wo.opdict['syn_datalength']=20.0
    wo.opdict['syn_samplefreq']=100.0
    wo.opdict['syn_kwidth']=0.1
    wo.opdict['syn_otime']=6.0
    wo.opdict['syn_ix']=16
    wo.opdict['syn_iy']=8
    wo.opdict['syn_iz']=6
    wo.opdict['syn_filename']='test_grid4D_hires.dat'

    wo.verify_migration_options()
    wo.verify_location_options()
    wo.verify_synthetic_options()

    ##########################
    # generate the test case and retrieve necessary information
    ##########################

    logging.info('Running synthetic test case generation...')
    test_info=generateSyntheticDirac(wo.opdict)
    logging.info(test_info)

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
   
#@unittest.skip('Not running real data migration tests')
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
 
