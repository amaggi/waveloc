import os, logging
from options import WavelocOptions
from synth_migration import generateSyntheticDirac
from plot_mpl import plotDiracTest


logging.basicConfig(level=logging.INFO, format='%(levelname)s : %(asctime)s : %(message)s')


# set up default parameters
wo = WavelocOptions()

# cheat - the test set is from PdF
# (for tests with other geometries must replicate the set_test_options parameters here)
wo.set_test_options()

base_path = os.getenv('WAVELOC_PATH')
wo.opdict['outdir'] = 'TEST_Dirac'
wo.opdict['search_grid']='grid.Taisne.search.hdr'
wo.opdict['loclevel'] = 10 
wo.opdict['load_ttimes_buf'] = True # Optimized in time, but you must be usre you're reading the right grid for the test
wo.opdict['syn_amplitude']=1.0
wo.opdict['syn_datalength']=20.0
wo.opdict['syn_samplefreq']=100.0
wo.opdict['syn_kwidth']=0.1
wo.opdict['syn_otime']=6.0
wo.opdict['syn_ix']=16
wo.opdict['syn_iy']=8
wo.opdict['syn_iz']=6
wo.opdict['syn_addnoise']=False
wo.opdict['syn_filename']='test_grid4D_hires.dat'

figdir=os.path.join(base_path,'out',wo.opdict['outdir'],'fig')

wo.verify_migration_options()
wo.verify_location_options()
wo.verify_synthetic_options()

#  No noise test
test_info=generateSyntheticDirac(wo.opdict)
test_info['grid_orig']=(0,0,-2.5)
plotDiracTest(test_info,figdir)

# Do noise tests

wo.opdict['outdir'] = 'TEST_DiracNoisy'
wo.opdict['syn_addnoise']=True
figdir=os.path.join(base_path,'out',wo.opdict['outdir'],'fig')

wo.opdict['syn_snr']=3.0
wo.opdict['syn_filename']='test_grid4D_hires_snr_3.0.dat'
test_info=generateSyntheticDirac(wo.opdict)
test_info['grid_orig']=(0,0,-2.5)
plotDiracTest(test_info,figdir)

