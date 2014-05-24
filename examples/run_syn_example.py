import os
import logging
from waveloc.options import WavelocOptions
from waveloc.synth_migration import generateSyntheticDirac
from waveloc.plot_mpl import plotDiracTest

logging.basicConfig(level=logging.INFO,
                    format='%(levelname)s : %(asctime)s : %(message)s')

# set up default parameters
wo = WavelocOptions()

# set base path to $WAVELOC_PATH
wo.verify_base_path()

# set options for synthetic test
wo.opdict['time'] = True
wo.opdict['verbose'] = False
wo.opdict['outdir'] = 'EXAMPLE_Dirac'
wo.opdict['time_grid'] = 'Slow_len.100m.P'
wo.opdict['load_ttimes_buf'] = True
wo.opdict['search_grid'] = 'grid.Taisne.search.hdr'
wo.opdict['stations'] = 'coord_stations_test'

#uncomment remove some stations to test response for fewer stations
#wo.opdict['sta_list'] = "FJS,FLR,FOR,HDL,RVL,SNE,UV01,UV02,UV03,UV04,UV05,\
#                         UV06,UV07,UV08,UV09,UV10,UV11,UV12,UV13,UV14,UV15"

wo.opdict['syn_amplitude'] = 1.0
wo.opdict['syn_datalength'] = 20.0
wo.opdict['syn_samplefreq'] = 100.0
wo.opdict['syn_kwidth'] = 0.1
wo.opdict['syn_otime'] = 6.0
wo.opdict['syn_ix'] = 16
wo.opdict['syn_iy'] = 8
wo.opdict['syn_iz'] = 6
wo.opdict['syn_filename'] = 'test_grid4D_hires.hdf5'

wo.opdict['plot_otime_window'] = 5.0

# sanity check on synthetic options
wo.verify_synthetic_options()

# Run the synthetic migration
test_info = generateSyntheticDirac(wo.opdict)

# Plot the synthetic migration
base_path = wo.opdict['base_path']
figdir = os.path.join(base_path, 'out', wo.opdict['outdir'], 'fig')
plotDiracTest(test_info, figdir, wo.opdict['plot_otime_window'])
