from waveloc.options import WavelocOptions
from waveloc.SDS_processing import do_SDS_processing_setup_and_run
from waveloc.migration import do_migration_setup_and_run
from waveloc.locations_trigger import do_locations_trigger_setup_and_run
from waveloc.plotting import do_plotting_setup_and_run

# set up default parameters
wo = WavelocOptions()

# set base path to $WAVELOC_PATH
wo.verify_base_path()

##########################################
# set waveloc options
##########################################

wo.opdict['time'] = True
wo.opdict['verbose'] = True
wo.opdict['ugrid_type'] = 'FULL'

wo.opdict['test_datadir'] = 'test_data'
wo.opdict['datadir'] = 'EXAMPLE'
wo.opdict['outdir'] = 'EXAMPLE_fullRes'

wo.opdict['net_list'] = 'YA'
wo.opdict['sta_list'] = "FJS,FLR,FOR,HDL,RVL,SNE,UV01,UV02,UV03,UV04,UV05,\
                         UV06,UV07,UV08,UV09,UV10,UV11,UV12,UV13,UV14,UV15"
wo.opdict['comp_list'] = "HHZ"

wo.opdict['starttime'] = "2010-10-14T00:14:00.0Z"
wo.opdict['endtime'] = "2010-10-14T00:18:00.0Z"

wo.opdict['time_grid'] = 'Slow_len.100m.P'
wo.opdict['search_grid'] = 'grid.Taisne.search.hdr'
wo.opdict['stations'] = 'coord_stations_test'

wo.opdict['resample'] = False
wo.opdict['fs'] = None

wo.opdict['c1'] = 4.0
wo.opdict['c2'] = 10.0

wo.opdict['kwin'] = 4
wo.opdict['krec'] = False
wo.opdict['kderiv'] = True

wo.opdict['data_length'] = 300
wo.opdict['data_overlap'] = 20

wo.opdict['dataglob'] = '*filt.mseed'
wo.opdict['kurtglob'] = '*kurt.mseed'
wo.opdict['gradglob'] = '*grad.mseed'

wo.opdict['load_ttimes_buf'] = True

wo.opdict['loclevel'] = 5000.0
wo.opdict['snr_limit'] = 10.0
wo.opdict['sn_time'] = 10.0
wo.opdict['n_kurt_min'] = 4

wo.opdict['plot_tbefore'] = 4
wo.opdict['plot_tafter'] = 6
wo.opdict['otime_window'] = 2

##########################################
# end of option setting - start processing
##########################################

# pre-process data from an SDS archive
wo.verify_SDS_processing_options()
do_SDS_processing_setup_and_run(wo.opdict)

# migrate
wo.verify_migration_options()
do_migration_setup_and_run(wo.opdict)

# detect and locate
wo.verify_location_options()
do_locations_trigger_setup_and_run(wo.opdict)

# plot results for located events
wo.verify_plotting_options()
do_plotting_setup_and_run(wo.opdict, plot_wfm=True, plot_grid=True)
