import os, logging, glob
import numpy as np
from options import *
from waveloc.SDS_processing import do_SDS_processing_setup_and_run
from waveloc.migration import do_migration_setup_and_run
from waveloc.locations_trigger import do_locations_trigger_setup_and_run 
from waveloc.plot_locations2 import do_plotting_setup_and_run

logging.basicConfig(level=logging.INFO, format='%(levelname)s : %(asctime)s : %(message)s')


# set up parameters parameters
wo = WavelocOptions()

base_path = os.getenv('WAVELOC_PATH')
wo.opdict['time']=True
wo.opdict['outdir'] = 'TEST_PdF'
wo.opdict['datadir'] = 'PdF'

wo.opdict['time_grid'] = 'Slow_len.100m.P'
wo.opdict['stations'] = 'coord_stations_piton'
wo.opdict['search_grid']='grid.Taisne.search.hdr'
wo.opdict['load_ttimes_buf'] = True # Optimized in time, but you must be usre you're reading the right grid for the test

wo.opdict['net_list']='YA'
wo.opdict['sta_list']="FJS,FLR,FOR,HDL,RVL,SNE,UV01,UV02,UV03,UV04,UV05,UV06,UV07,UV08,UV09,UV10,UV11,UV12,UV13,UV14,UV15"
wo.opdict['comp_list']="HHZ"

wo.opdict['starttime']="2010-10-14T00:00:00.0Z"
wo.opdict['endtime']="2010-10-14T16:00:00.0Z"

#wo.opdict['starttime']="2010-10-14T04:00:00.0Z"
#wo.opdict['endtime']="2010-10-14T05:00:00.0Z"

wo.opdict['resample']=False
wo.opdict['fs']=None

wo.opdict['c1']=4.0
wo.opdict['c2']=10.0

wo.opdict['kwin']=4
wo.opdict['krec']=True
wo.opdict['kderiv']=True

wo.opdict['data_length']=130
wo.opdict['data_overlap']=10

wo.opdict['dataglob']='*filt.mseed'
wo.opdict['kurtglob']='*kurt.mseed'
wo.opdict['gradglob']='*grad.mseed'

wo.opdict['auto_loclevel']=False
wo.opdict['snr_loclevel']=100
wo.opdict['loclevel']=50
wo.opdict['snr_limit']=10.0
wo.opdict['sn_time']=10.0
wo.opdict['n_kurt_min']=4.0

wo.opdict['plot_tbefore']=10.0
wo.opdict['plot_tafter']=20.0


# check processing options and run processing
wo.verify_SDS_processing_options()
#do_SDS_processing_setup_and_run(wo.opdict)

# check migration options and run migration
wo.verify_migration_options()
#do_migration_setup_and_run(wo.opdict)

# check location options and run location
wo.verify_location_options()
#locs=do_locations_trigger_setup_and_run(wo.opdict)

# check plotting options and run plotting
wo.verify_plotting_options()
do_plotting_setup_and_run(wo.opdict)

