#!/bin/bash

# Script to run waveloc on the tutorial data
data_dir="Tutorial"
net_list="DJ,DO"
sta_list="ARO,AROG,ATA,ATD,FIK,GBR,KSD,LDL,MAO,MCA,OBO,SGH,TDD,ASE,BLS,DAF2,DAY,FIEA,MOY2,PCKD,PTRN,RAND,SBLV,WADA,WAY"
comp_list="HHZ,BHZ"
out_dir="Tutorial"
start_time="2010-04-20T00:00:00.0Z"
end_time="2010-04-21T00:00:00.0Z"

time_grid=tutorial.P
search_grid=grid.tutorial.search.hdr
coord_stations=coord_stations_tutorial

mkdir -p $WAVELOC_PATH_PDF/out/$out_dir
mkdir -p $WAVELOC_PATH_PDF/out/$out_dir/grid
mkdir -p $WAVELOC_PATH_PDF/out/$out_dir/stack
mkdir -p $WAVELOC_PATH_PDF/out/$out_dir/loc
mkdir -p $WAVELOC_PATH_PDF/out/$out_dir/reloc


#./run_SDS_processing_threading.py -n 10 --datadir=$data_dir --net_list=$net_list --sta_list=$sta_list --comp_list=$comp_list  --starttime=$start_time --endtime=$end_time --c1=6 --c2=10 --kwin=5.0 --krec --resample --fs=40 --kderiv

#./run_waveloc_threading.py -t -v -n 1 --time_grid $time_grid --search_grid $search_grid -s $coord_stations -o $out_dir --datadir=$data_dir --dataglob=*kurt_grad.mseed --starttime=$start_time --endtime=$end_time --data_length=600 --data_overlap=20
#./run_waveloc_threading.py -t -v -n 1 --time_grid $time_grid --search_grid $search_grid -s $coord_stations -o $out_dir --datadir=$data_dir --dataglob=*kurt_grad.mseed --starttime=$start_time --endtime=$end_time --data_length=600 --data_overlap=20 --load_ttimes_buf 

./locations_trigger.py --outdir=$out_dir --loclevel=50 --datadir=$data_dir --dataglob=*kurt.mseed --n_kurt_min=4
#./locations_prob.py --outdir=$out_dir --loclevel=100 --datadir=$data_dir --dataglob=*kurt.mseed --n_kurt_min=4


xvfb-run -e err_file.txt --server-args="-screen 0 1024x768x24" ./plot_locations.py --run_mayavi --max_stack=120 --outdir=$out_dir --datadir=$data_dir --search_grid=$search_grid --data_glob=*filt.mseed --kurt_glob=*kurt.mseed --grad_glob=*kurt_grad.mseed --time_grid $time_grid --stations $coord_stations


