#!/bin/sh

first="125127"
last="125127"
CONFIG="./config/config.cfg"
RUNLIST="./runlist/runlist.csv"

alias js="jobsub -c $CONFIG -csv $RUNLIST"

for i in `seq $first $last`; do
         js noisypixel $i 
         #js clustering $i
         #js hitmaker $i 
         #js patternRecognition $i
         #js GBLAlign $i
         #js patternRecognition $i      
         #js GBLTrackFit $i

done

