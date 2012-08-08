#!/bin/bash

# Script to run waveloc on the test dataset for Piton de la Fournaise
data_dir="PdF"
net_list="YA"
sta_list="FJS,FLR,FOR,HDL,RVL,SNE,UV01,UV02,UV03,UV04,UV05,UV06,UV07,UV08,UV09,UV10,UV11,UV12,UV13,UV14,UV15"
comp_list="HHZ"
out_dir="PdF_test"
start_time="2010-10-14T00:00:00.0Z"
end_time="2010-10-14T01:00:00.0Z"

snr_limit=10.0
n_kurt_min=4

time_grid=Slow_len.100m.P
search_grid=grid.Taisne.search.hdr
coord_stations=coord_stations_piton

mkdir -p $WAVELOC_PATH/out/$out_dir
mkdir -p $WAVELOC_PATH/out/$out_dir/grid
mkdir -p $WAVELOC_PATH/out/$out_dir/stack
mkdir -p $WAVELOC_PATH/out/$out_dir/loc
mkdir -p $WAVELOC_PATH/out/$out_dir/reloc

./run_SDS_processing_threading.py -n 10 --datadir=$data_dir --net_list=$net_list --sta_list=$sta_list --comp_list=$comp_list  --starttime=$start_time --endtime=$end_time --c1=4 --c2=10 --kwin=5.0 --krec --resample --fs=50 --kderiv

#./run_waveloc_threading.py -t -v -n 1 --time_grid $time_grid --search_grid $search_grid -s $coord_stations -o $out_dir --datadir=$data_dir --dataglob=*kurt_grad.mseed --starttime=$start_time --endtime=$end_time --data_length=600 --data_overlap=20
#./run_waveloc_threading.py -t -v -n 4 --time_grid $time_grid --search_grid $search_grid -s $coord_stations -o $out_dir --datadir=$data_dir --dataglob=*kurt_grad.mseed --starttime=$start_time --endtime=$end_time --data_length=600 --data_overlap=20 --load_ttimes_buf 

#./locations_trigger.py --outdir=$out_dir --loclevel=50 --datadir=$data_dir --dataglob=*kurt.mseed --n_kurt_min=$n_kurt_min --snr_limit=$snr_limit
#./locations_prob.py --outdir=$out_dir --loclevel=50 --datadir=$data_dir --dataglob=*kurt.mseed --n_kurt_min=$n_kurt_min --snr_limit=$snr_limit

#xvfb-run -e err_file.txt --server-args="-screen 0 1024x768x24" ./plot_locations.py --run_mayavi --max_stack=120 --outdir=$out_dir --datadir=$data_dir --search_grid=$search_grid --data_glob=*filt.mseed --kurt_glob=*kurt.mseed --grad_glob=*kurt_grad.mseed --time_grid $time_grid --stations $coord_stations --snr_limit=$snr_limit



#./plot_only_dem_mayavi.py
