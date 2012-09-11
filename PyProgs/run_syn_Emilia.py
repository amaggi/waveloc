import os, logging
import numpy as np
from options import WavelocOptions
from synth_migration import generateSyntheticDirac
from grids_paths import QDGrid
from NllGridLib import *
from plot_mpl import plotDiracTest


logging.basicConfig(level=logging.INFO, format='%(levelname)s : %(asctime)s : %(message)s')


# set up default parameters
wo = WavelocOptions()

# cheat - the test set is from PdF
# (for tests with other geometries must replicate the set_test_options parameters here)

base_path = os.getenv('WAVELOC_PATH')
wo.opdict['outdir'] = 'TEST_DiracEmilia'
wo.opdict['time_grid'] = 'emilia.P'
wo.opdict['stations'] = 'coord_stations_temp_emilia'
wo.opdict['search_grid']='grid.emilia.search.hdr'
wo.opdict['loclevel'] = 10 
wo.opdict['load_ttimes_buf'] = True # Optimized in time, but you must be usre you're reading the right grid for the test
wo.opdict['syn_amplitude']=1.0
wo.opdict['syn_datalength']=50.0
wo.opdict['syn_samplefreq']=100.0
wo.opdict['syn_kwidth']=0.1
wo.opdict['syn_otime']=5.0
wo.opdict['syn_addnoise']=False


# set up basic grid information for test
grid_filename=os.path.join(base_path,'lib',wo.opdict['search_grid'])
dummy_grid=QDGrid()
dummy_grid.read_NLL_hdr_file(grid_filename)
# set up projection information for test
f=open(grid_filename)
lines=f.readlines()
f.close()
proj_line=lines[1]
proj_info={}
proj_info['orig_lon'] = np.float(proj_line.split()[3])
proj_info['orig_lat'] = np.float(proj_line.split()[5])
proj_info['map_rot'] = np.float(proj_line.split()[7])

print proj_info


event1ll=(44.9235 , 11.1418)
event2ll=(44.8833 , 11.1350)
event3ll=(44.7892 , 11.1705)
event4ll=(44.5860 , 10.9193)
event5ll=(44.8562 , 11.4182)

depths=[3.0, 5.0, 10.0]
events={}
events['ev1']=event1ll
events['ev2']=event2ll
events['ev3']=event3ll
events['ev4']=event4ll
events['ev5']=event5ll

for evname,ev in events.iteritems():
  for dep in depths:
    x,y=latlon2rect('TRANS_SIMPLE',ev[0],ev[1],proj_info)
    
    wo.opdict['syn_ix']=int(round((x-dummy_grid.x_orig)/dummy_grid.dx))
    wo.opdict['syn_iy']=int(round((y-dummy_grid.y_orig)/dummy_grid.dy))
    wo.opdict['syn_iz']=int(round((dep-dummy_grid.z_orig)/dummy_grid.dz))
    wo.opdict['syn_filename']='test_emilia_%s_%.1f.dat'%(evname,dep)
    #  No noise test
    wo.verify_synthetic_options()
    test_info=generateSyntheticDirac(wo.opdict)
    base_path=os.getenv('WAVELOC_PATH')
    figdir=os.path.join(base_dir,'out',wo.opdict['outdir'],'fig')
    plotDiracTest(test_info,figdir)
    exit()

    
