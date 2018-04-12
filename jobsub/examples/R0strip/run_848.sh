#!/bin/sh

first="000848"
last="000848"
CONFIG="./config/config_848.cfg"
RUNLIST="./runlist/runlist.csv"

alias js="jobsub -c $CONFIG -csv $RUNLIST"

for i in `seq $first $last`; do
         js noisypixel $i 
	 js clustering $i
	 js hitmaker $i 
	 js patternRecognition $i
         #js GBLAlign $i
         #js patternRecognition $i 	
	 js GBLTrackFit $i

done
