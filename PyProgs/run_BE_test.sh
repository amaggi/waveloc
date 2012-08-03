#!/bin/bash

# Script to run waveloc_pdf on the september 23 2010 crisis
data_dir="BW"
net_list="BE"
sta_list="CLA,DOU,GES,GRZ,OT1,OT2,OT3,OT5,OTT,RQR,SKQ,UCC"
comp_list="HHZ"
out_dir="BW_test_3D"
start_time="2009-03-03T03:00:00.0Z"
end_time="2009-03-03T04:00:00.0Z"

time_grid=belgium.P
search_grid=grid.belgium.search.hdr

mkdir -p $WAVELOC_PATH_PDF/out/$out_dir
mkdir -p $WAVELOC_PATH_PDF/out/$out_dir/grid
mkdir -p $WAVELOC_PATH_PDF/out/$out_dir/stack
mkdir -p $WAVELOC_PATH_PDF/out/$out_dir/loc
mkdir -p $WAVELOC_PATH_PDF/out/$out_dir/reloc


#./run_SDS_processing_threading.py --datadir=$data_dir --net_list=$net_list --sta_list=$sta_list --comp_list=$comp_list  --starttime=$start_time --endtime=$end_time --c1=5 --c2=25 --kwin=1.0 --krec --resample --fs=100 --kderiv

#./run_waveloc_threading.py -t -v --time_grid $time_grid --search_grid $search_grid -s coord_stations_belgium -o $out_dir --datadir=$data_dir --dataglob=*kurt.mseed --starttime=$start_time --endtime=$end_time --data_length=2000 --data_overlap=20 

#./locations_trigger.py --outdir=$out_dir --loclevel=100 --datadir=$data_dir --dataglob=*kurt.mseed --n_kurt_min=4
./locations_prob.py --outdir=$out_dir --loclevel=100 --datadir=$data_dir --dataglob=*kurt.mseed --n_kurt_min=4


#xvfb-run -e err_file.txt --server-args="-screen 0 1024x768x24" ./plot_locations.py --run_mayavi --max_stack=2000 --outdir=$out_dir --datadir=$data_dir --search_grid=$search_grid --data_glob=*filt.mseed --kurt_glob=*kurt.mseed --time_grid $time_grid --stations coord_stations_belgium


#./plot_only_dem_mayavi.py
