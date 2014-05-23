import os, logging, h5py
import numpy as np
from random import randint
from waveloc.NllGridLib import read_hdr_file
from waveloc.options import WavelocOptions
from waveloc.synth_migration import generateSyntheticDirac
from waveloc.locations_trigger import trigger_locations_inner

logging.basicConfig(level=logging.INFO, format='%(levelname)s : %(asctime)s : %(message)s')

def setUp() :

    logging.info('Setting up synthetic test case generation...')

    # set up default parameters
    wo = WavelocOptions()

    # set base path to $WAVELOC_PATH
    wo.verify_base_path()

    # set options for synthetic test
    wo.opdict['time']    = True
    wo.opdict['verbose'] = False

    wo.opdict['outdir'] = 'TEST_Dirac'
    wo.opdict['time_grid'] = 'Slow_len.100m.P'
    wo.opdict['load_ttimes_buf'] = True
    wo.opdict['search_grid'] = 'grid.Taisne.search.hdr'
    wo.opdict['stations'] = 'coord_stations_test'

    wo.opdict['syn_filename']   = 'test_grid4D_hires.hdf5'
    wo.opdict['syn_amplitude']  = 1.0
    wo.opdict['syn_datalength'] = 20.0
    wo.opdict['syn_samplefreq'] = 100.0
    wo.opdict['syn_kwidth']     = 0.1
    wo.opdict['syn_otime']      = 6.0

    # place default point at center of grid
    base_path = wo.opdict['base_path']
    search_grid_filename = os.path.join(base_path,'lib',wo.opdict['search_grid'])
    grid_info = read_hdr_file(search_grid_filename)
    wo.opdict['syn_ix'] = grid_info['nx']/2
    wo.opdict['syn_iy'] = grid_info['ny']/2
    wo.opdict['syn_iz'] = grid_info['nz']/2

    # sanity check
    wo.verify_synthetic_options()

    return (wo,grid_info)

def doPointTest(wo,loclevel=10.0) :
    
    logging.info('Doing synthetic test for point (%d,%d,%d)...'\
        %(wo.opdict['syn_ix'],wo.opdict['syn_iy'],wo.opdict['syn_iz']))

    # do the migration
    test_info = generateSyntheticDirac(wo.opdict)
    logging.info(test_info)

    # retrieve output info
    stack_filename = test_info['stack_file']
    nx,ny,nz,nt = test_info['grid_shape']
    dx,dy,dz,dt = test_info['grid_spacing']
    x_orig,y_orig,z_orig = test_info['grid_orig']
    ix_true,iy_true,iz_true,it_true = test_info['true_indexes']
    stack_start_time = test_info['start_time']

    # set up x, y, z, t arrays
    x = np.arange(nx)*dx
    y = np.arange(ny)*dy
    z = np.arange(nz)*dz
    t = np.arange(nt)*dt+stack_start_time

    # extract the max stacks
    f_stack=h5py.File(stack_filename,'r')
    max_val = f_stack['max_val']
    max_x = f_stack['max_x']
    max_y = f_stack['max_y']
    max_z = f_stack['max_z']

    # launch a location trivver
    locs=trigger_locations_inner(max_val,max_x,max_y,max_z,\
        loclevel,loclevel,stack_start_time,dt)
    
    f_stack.close()
    return test_info,locs

def analyseLocs(locs,wo,test_info) :

    dx,dy,dz,dt = test_info['grid_spacing']
    x_orig,y_orig,z_orig = test_info['grid_orig']
    ix_true,iy_true,iz_true,it_true = test_info['true_indexes']

    # This is a dirac test, but we may have more than one loc (secondary maxima)
    # so for safety, pull out best loc


    n_locs = len(locs)
    if n_locs > 0 :
        imax = np.argmax([loc['max_trig'] for loc in locs])
        trig_loc = locs[imax]
        print trig_loc['max_trig']

        loc_dist = np.sqrt( (ix_true*dx+x_orig - trig_loc['x_mean'])**2 +\
                            (iy_true*dy+y_orig - trig_loc['y_mean'])**2 +\
                            (iz_true*dy+z_orig - trig_loc['z_mean'])**2)
        loc_dt = trig_loc['o_time'] - wo.opdict['syn_otime']
        loc_dt_list = [aloc['o_time'] - wo.opdict['syn_otime'] for aloc in locs]
        print loc_dt_list
    else :
	loc_dist = None
        loc_dt = None
  

    return n_locs,loc_dist,loc_dt


if __name__ == '__main__' :


    wo,grid_info = setUp()

    # get information from search grid
    nx = grid_info['nx']
    ny = grid_info['ny']
    nz = grid_info['nz']

    # set up random test point
    ix = randint(0,nx)
    iy = randint(0,ny)
    iz = randint(0,nz)
    wo.opdict['syn_ix'] = ix
    wo.opdict['syn_iy'] = iy
    wo.opdict['syn_iz'] = iz

    # do synthetic test for this point
    test_info, locs = doPointTest(wo,loclevel=10.0)
    n_locs,best_dist,best_dt = analyseLocs(locs,wo,test_info)

    print ix,iy,iz,n_locs,best_dist,best_dt

