#!/bin/bash

#   Runs the SSL2 EISMINT-Greenland experiment.  Requires green20km_Tsteady.nc,
# the output of bootstrap.sh.  Saves state every 10000 model years.  Use
# series.py to see volume time series to evaluate EISMINT-Greenland
# steady state criterion of 0.01% change in volume in 10000 model years.
#   See the PISM User's Manual for the modeling choices in "pgrn -ssl2".

#   A recommended way to run with 2 processors is "./ssl2.sh 2 >& ssl2.txt &"
# which saves a transcript in ssl2.txt.

NN=2  # default number of processors
if [ $# -gt 0 ] ; then  # if user says "ssl2.sh 8" then NN = 8
  NN="$1"
fi
set -e  # exit on error

mpiexec -n $NN pgrn -ssl2 -if green20km_Tsteady.nc \
       -ys 0 -y 10000 -o green_SSL2_10k

for ((kyear=20; kyear <= 110 ; kyear+=10))
do
  (( oldkyear = kyear - 10 ))
  mpiexec -n $NN pgrn -ssl2 -if green_SSL2_${oldkyear}k.nc \
       -y 10000 -o green_SSL2_${kyear}k
done


