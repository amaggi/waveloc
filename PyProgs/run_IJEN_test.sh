#!/bin/bash

# Script to run waveloc_pdf on the september 23 2010 crisis
data_dir="NewDataIjen"
net_list="ID"
sta_list="IJEN,KWUI,MLLR,POS,POSI,PSG,RAUN,TRWI"
comp_list="EHZ,BHZ"
out_dir="IJEN_test_grid"
start_time="2011-10-06T00:00:00.0Z"
end_time="2011-10-07T00:00:00.0Z"

search_grid=grid.ijen.search.hdr

mkdir -p $WAVELOC_PATH_PDF/out/$out_dir
mkdir -p $WAVELOC_PATH_PDF/out/$out_dir/grid
mkdir -p $WAVELOC_PATH_PDF/out/$out_dir/stack
mkdir -p $WAVELOC_PATH_PDF/out/$out_dir/loc
mkdir -p $WAVELOC_PATH_PDF/out/$out_dir/reloc


#./run_SDS_processing_threading.py --datadir=$data_dir --net_list=$net_list --sta_list=$sta_list --comp_list=$comp_list  --starttime=$start_time --endtime=$end_time --c1=4 --c2=10 --kwin=3.0

./run_PdF_waveloc.py -t -v --time_grid ijen.P --search_grid $search_grid -s coord_stations_ijen -o $out_dir --datadir=$data_dir --dataglob=*kurt.mseed --starttime=$start_time --endtime=$end_time --data_length=3620 --data_overlap=20 --comp_list=$comp_list

./locations_trigger.py --outdir=$out_dir --loclevel=50 --datadir=$data_dir --dataglob=*kurt.mseed --n_kurt_min=3 


xvfb-run --server-args="-screen 0 1024x768x24" ./plot_locations.py --run_mayavi --max_stack=100 --outdir=$out_dir --datadir=$data_dir --search_grid=$search_grid --data_glob=*filt.mseed --kurt_glob=*kurt.mseed --time_grid ijen.P --stations coord_stations_ijen   --comp_list=$comp_list


#./plot_only_dem_mayavi.py
