#!/bin/bash

# Script to run waveloc_pdf on the september 23 2010 crisis
data_dir="LAquila"
net_list="IV,MN,TV"
sta_list="CAFR,CAMP,CERA,CERT,CESI,CESX,FAGN,FIAM,GUAR,GUMA,INTR,LNSS,LPEL,MA9,MNS,MTCE,NRCA,OFFI,POFI,RM01,RM02,RM03,RM04,RM05,RM06,RM07,RM08,RM09,RM10,RM11,RM13,RM14,RM18,RM20,RM21,RM22,RM23,RM24,RMP,RNI2,ROM9,TERO,TRTR,AQU,T0101,T0102,T0103,T0105,T0106,T0108"
comp_list="HHZ,EHZ,HHZ,BNZ"
out_dir="LAquila"
start_time="2009-04-06T00:00:00.0Z"
end_time="2009-04-06T10:00:00.0Z"

search_grid=grid.ijen.search.hdr

mkdir -p $WAVELOC_PATH_PDF/out/$out_dir
mkdir -p $WAVELOC_PATH_PDF/out/$out_dir/grid
mkdir -p $WAVELOC_PATH_PDF/out/$out_dir/stack
mkdir -p $WAVELOC_PATH_PDF/out/$out_dir/loc
mkdir -p $WAVELOC_PATH_PDF/out/$out_dir/reloc


./run_SDS_processing_threading.py --datadir=$data_dir --net_list=$net_list --sta_list=$sta_list --comp_list=$comp_list  --starttime=$start_time --endtime=$end_time --c1=1 --c2=4 --kwin=15.0

#./run_PdF_waveloc.py -t -v --time_grid ijen.P --search_grid $search_grid -s coord_stations_ijen -o $out_dir --datadir=$data_dir --dataglob=*kurt.mseed --starttime=$start_time --endtime=$end_time --data_length=920 --data_overlap=20 

#./locations_trigger.py --outdir=$out_dir --loclevel=50 --datadir=$data_dir --dataglob=*kurt.mseed


#xvfb-run --server-args="-screen 0 1024x768x24" ./plot_locations.py --run_mayavi --max_stack=200 --outdir=$out_dir --datadir=$data_dir --search_grid=$search_grid --data_glob=*filt.mseed --kurt_glob=*kurt.mseed --time_grid ijen.P --stations coord_stations_ijen


#./plot_only_dem_mayavi.py
