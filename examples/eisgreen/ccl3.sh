#!/bin/bash

#   Runs the CCL3 EISMINT-Greenland experiment.  Requires GRIP ice core 
# temperature data in grip_dT.nc and SPECMAP sea bed core (stack) sea level
# data in specmap_dSL.nc (from preprocess.sh) and final result from ssl.sh.
# See the PISM User's Manual.

#   A recommended way to run with 2 processors is "./ccl3.sh 2 >& out.ccl3 &"
# which saves a transcript in out.ccl3

NN=2  # default number of processors
if [ $# -gt 0 ] ; then  # if user says "ccl3.sh 8" then NN = 8
  NN="$1"
fi

set -e  # exit on error

SHOWONLY=0
if [ $# -gt 1 ] ; then  # if user says "ccl3.sh 8 D" then *only* shows; debug mode
  SHOWONLY=1
fi

DTNAME=grip_dT.nc
DSLNAME=specmap_dSL.nc

# function to run "pgrn -ccl3 -dTforcing ... -dSLforcing ..." on NN processors
mpgrn()
{
    # change this if "mpirun" or "bin/pgrn", etc.:
    cmd="mpiexec -n $NN pgrn -ccl3 -dTforcing $DTNAME -dSLforcing $DSLNAME $1"
    
    echo 
    echo "date = '`date`' on host '`uname -n`':"
    echo "trying '$cmd'"
    if [ $SHOWONLY = 0 ] ; then
      $cmd
      echo
    fi
}

# the EISMINT-Greenland CCL3 experiment:

mpgrn "-i green_SSL2_110k.nc -ys -249900 -y 9900 -o green_CCL3_y-240k.nc"

for ((kyear=230; kyear >= 0 ; kyear-=10)); do
  (( oldkyear = kyear + 10 ))
  mpgrn "-i green_CCL3_y-${oldkyear}k.nc -y 10000 -o green_CCL3_y-${kyear}k.nc"
done

