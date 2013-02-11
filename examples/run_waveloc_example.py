import os, logging
from waveloc.options import WavelocOptions
from waveloc.SDS_processing import do_SDS_processing_setup_and_run
from waveloc.migration import do_migration_setup_and_run
from waveloc.locations_trigger import do_locations_trigger_setup_and_run
from waveloc.locations_prob import do_locations_prob_setup_and_run
from waveloc.plot_locations2 import do_plotting_setup_and_run, \
                                    do_probloc_plotting_setup_and_run


logging.basicConfig(level=logging.INFO, 
        format='%(levelname)s : %(asctime)s : %(message)s')


# set up default parameters
wo = WavelocOptions()

# set base path to $WAVELOC_PATH 
wo.verify_base_path()

##########################################
# set waveloc options
##########################################

wo.opdict['time']=True
wo.opdict['verbose']=False

wo.opdict['test_datadir']='test_data'
wo.opdict['datadir']='TEST'
wo.opdict['outdir']='TEST_fullRes'

wo.opdict['net_list']='YA'
wo.opdict['sta_list']="FJS,FLR,FOR,HDL,RVL,SNE,UV01,UV02,UV03,UV04,UV05,UV06,\
        UV07,UV08,UV09,UV10,UV11,UV12,UV13,UV14,UV15"
wo.opdict['comp_list']="HHZ"

wo.opdict['starttime']="2010-10-14T00:14:00.0Z"
wo.opdict['endtime']="2010-10-14T00:18:00.0Z"

wo.opdict['time_grid']='Slow_len.100m.P'
wo.opdict['search_grid'] = 'grid.Taisne.search.hdr'
wo.opdict['stations']='coord_stations_test'

wo.opdict['resample']=False
wo.opdict['fs']=None

wo.opdict['c1']=4.0
wo.opdict['c2']=10.0

wo.opdict['kwin']=4
wo.opdict['krec']=False
wo.opdict['kderiv']=True

wo.opdict['data_length']=300
wo.opdict['data_overlap']=20

wo.opdict['dataglob']='*filt.mseed'
wo.opdict['kurtglob']='*kurt.mseed'
wo.opdict['gradglob']='*grad.mseed'

wo.opdict['load_ttimes_buf']=True

wo.opdict['loclevel']=50.0
wo.opdict['snr_limit']=10.0
wo.opdict['sn_time']=10.0
wo.opdict['n_kurt_min']=4

wo.opdict['plot_tbefore']=4
wo.opdict['plot_tafter']=6
wo.opdict['plot_otime_window']=5

##########################################
# end of option setting - start processing
##########################################

wo.verify_SDS_processing_options()
do_SDS_processing_setup_and_run(wo.opdict)

wo.verify_migration_options()
do_migration_setup_and_run(wo.opdict)

# do trigger location 
wo.verify_location_options()
do_locations_trigger_setup_and_run(wo.opdict)

# This will do plotting of grids and stacks for locations
wo.verify_plotting_options()
do_plotting_setup_and_run(wo.opdict,plot_wfm=True,plot_grid=True)

