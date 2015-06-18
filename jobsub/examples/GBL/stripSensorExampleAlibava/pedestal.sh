#!/bin/sh

################################
#                              #
# ALiBaVa Analysis - Pedestal  #
#                              #
################################

# usage: sh pedestal.sh runnumber

#To calculate Pedestals
jobsub -c config/config.cfg -csv runlist/runlist.csv alibava-convert-ped $1
jobsub -c config/config.cfg -csv runlist/runlist.csv alibava-pedestal $1
jobsub -c config/config.cfg -csv runlist/runlist.csv alibava-commonmode $1
jobsub -c config/config.cfg -csv runlist/runlist.csv alibava-pedestal2 $1

