#!/bin/bash

# Script to run waveloc_pdf on the september 23 2010 crisis
data_dir="2010_266"
out_dir="2010_266"
start_time="2010-09-23T22:00:00.0Z"
#end_time="2010-09-23T22:10:00.0Z"
end_time="2010-09-24T00:00:00.0Z"
search_grid=grid.Taisne.search.hdr
search_grid_reloc=grid.Taisne.search.hdr

mkdir -p $WAVELOC_PATH_PDF/out/$out_dir
mkdir -p $WAVELOC_PATH_PDF/out/$out_dir/grid
mkdir -p $WAVELOC_PATH_PDF/out/$out_dir/stack
mkdir -p $WAVELOC_PATH_PDF/out/$out_dir/loc
mkdir -p $WAVELOC_PATH_PDF/out/$out_dir/reloc


#./run_PdF_processing_cleanup.py --datadir=$data_dir --dataglob=*.mseed --starttime=$start_time --endtime=$end_time --short_data_length=600 --c1=4 --c2=10 --kwin=3.0

./run_waveloc_threading.py -t -v --time_grid Slow_len.100m.P --search_grid $search_grid -s coord_stations_piton -o $out_dir --datadir=$data_dir --dataglob=*kurt.sac --starttime=$start_time --endtime=$end_time --data_length=200 --data_overlap=20

./locations_trigger.py --outdir=$out_dir --datadir=$data_dir --dataglob=*kurt.sac --loclevel=200


xvfb-run --server-args="-screen 0 1024x768x24" ./plot_locations.py --run_mayavi --max_stack=200 --outdir=$out_dir --datadir=$data_dir --search_grid=$search_grid --data_glob=*filt.sac --kurt_glob=*kurt.sac --time_grid Slow_len.100m.P --stations coord_stations_piton 
#./plot_locations.py --run_mayavi --max_stack=200 --outdir=$out_dir --datadir=$data_dir --search_grid=$search_grid --data_glob=*filt.sac --kurt_glob=*kurt.sac --time_grid Slow_len.100m.P --stations coord_stations_piton 



#./plot_only_dem_mayavi.py
