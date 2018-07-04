#!/bin/sh

first="001304"
last="001304"
CONFIG="./config/config.cfg"
RUNLIST="./runlist/runlist.csv"

alias js="jobsub -c $CONFIG -csv $RUNLIST"

for i in `seq $first $last`; do
         js noisypixel $i 
	 js clustering $i
	 js hitmaker $i 
	 js patternRecognition $i
         #js GBLAlign $i
	 js GBLTrackFit $i

done
