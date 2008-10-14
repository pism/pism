#!/bin/bash

#   Runs the SSL2 EISMINT-Greenland experiment from the output of bootstrap.sh.
# Saves state every 10000 model years.  See the PISM User's Manual.

#   Recommended way to run with 2 processors is "./ssl2.sh 2 >& out.ssl2 &"
# which saves a transcript in out.ssl2

NN=2  # default number of processors
if [ $# -gt 0 ] ; then  # if user says "ssl2.sh 8" then NN = 8
  NN="$1"
fi

set -e  # exit on error

SHOWONLY=0
if [ $# -gt 1 ] ; then  # if user says "ssl2.sh 8 D" then NN = 8 and *only* shows
                        #   what will happen; NO RUN (debug mode)
  SHOWONLY=1
fi

# function to run "pgrn -ssl2" on NN processors
mpgrn()
{
    cmd="mpiexec -n $NN pgrn -ssl2 $1"  # change if "mpirun" or "bin/pgrn", etc.
    
    echo 
    echo "date = '`date`' on host '`uname -n`':"
    echo "trying '$cmd'"
    if [ $SHOWONLY = 0 ] ; then
      $cmd
      echo
    fi
}

# the EISMINT-Greenland SSL2 experiment:

mpgrn "-if green20km_Tsteady.nc -ys 0 -y 10000 -o green_SSL2_10k"

for ((kyear=20; kyear <= 110 ; kyear+=10)); do
  (( oldkyear = kyear - 10 ))
  mpgrn "-if green_SSL2_${oldkyear}k.nc -y 10000 -o green_SSL2_${kyear}k"
done

