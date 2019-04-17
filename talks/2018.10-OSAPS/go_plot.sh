#!/bin/bash

option=$1
#option=0 means "short"
#option=1 means "all 2013"


if [ $option -eq 0 ]
then
	echo "Running single run"
	./plot_stuff 2 /data/user/pfendner/A23/outputs/Processed_final/processed_station_2_run_1798_joined_bins_6_19.root
elif [ $option -eq 2 ]
then
	echo "Running medium 2013"
	./plot_stuff 2 /data/user/pfendner/A23/outputs/Processed_final/processed_station_2_run_14*_joined_bins_6_19.root
fi