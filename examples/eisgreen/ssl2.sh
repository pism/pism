#!/bin/bash

#   Runs the SSL2 EISMINT-Greenland experiment from the output of bootstrap.sh.
# Saves state every 10000 model years.  See the PISM User's Manual.

#   Recommended way to run with 8 processors is "./ssl2.sh 8 >& out.ssl2 &"
# which saves a transcript in out.ssl2

NN=2  # default number of processors
if [ $# -gt 0 ] ; then  # if user says "ssl2.sh 8" then NN = 8
  NN="$1"
fi

set -e  # exit on error

SHOWONLY=0
if [ $# -gt 1 ] ; then  # if user says "./ssl2.sh 8 D" then NN = 8 and *only* shows
                        #   what will happen; NO RUN (debug mode)
  SHOWONLY=1
fi

# function to run "pgrn -ssl2 -skip 10" on NN processors
mpgrn()
{
    cmd="mpiexec -n $NN pgrn -ssl2 -skip 10 $1"  # change if "mpirun" or "bin/pgrn", etc.
    
    echo 
    echo "date = '`date`' on host '`uname -n`':"
    echo "trying '$cmd'"
    if [ $SHOWONLY = 0 ] ; then
      $cmd
      echo
    fi
}

# the EISMINT-Greenland SSL2 experiment:

INAME=green20km_Tsteady.nc
RESET="-ys 0"
for ((kyear=10; kyear <= 250 ; kyear+=10)); do
  (( oldkyear = kyear - 10 ))
  ENDNAME=green_SSL2_${kyear}k.nc
  mpgrn "-i ${INAME} ${RESET} -y 10000 -ts_file ssl2vol_${kyear}k.nc -ts_times ${oldkyear}000:1000:${kyear}000 -ts_vars ivol -o ${ENDNAME}"
  INAME=${ENDNAME}
  RESET=""
done

# use ncrcat (an NCO = NetCDF Operator), if available, to build single time series for ice volume:
# ncrcat -o ssl2vol_all.nc ssl2vol_?0k.nc ssl2vol_??0k.nc

