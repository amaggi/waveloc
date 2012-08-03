#!/bin/bash

export WAVELOC_PATH_PDF="$HOME/am_test"

data_dir=14_oct
out_dir=testing

./extract_located_events.py -o $out_dir -d $data_dir --data_glob *.filt.sac --left_pad 10 --right_pad 20

