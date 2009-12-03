#!/bin/bash

# Runs the SSL2 EISMINT-Greenland experiment, starting from the output file
# produced by bootstrap.sh.  Saves state every 10000 model years, and saves a timeseries of
# volume every 100 years.  See the PISM User's Manual.
#
# Recommended way to run with 8 processors is "./ssl2.sh 8 >& out.ssl2 &",
# which saves a transcript in out.ssl2.

NN=2  # default number of processors
if [ $# -gt 0 ] ; then  # if user says "./ssl2.sh 8" then NN = 8
  NN="$1"
fi

SHOW=
if [ $# -gt 1 ] ; then  # if user says "./ssl2.sh 8 D" then NN = 8 and *only* shows
                        #   what will happen; NO RUN (debug mode)
  SHOW="echo "
fi

set -e  # exit on error

MPIDO="mpiexec -n"

INFILE=green20km_Tsteady.nc

snaps="-save_file snaps_ssl2.nc -save_times 10000:10000:100000"

ts="-ts_file vol_ssl2.nc -ts_vars ivol -ts_times 100:100:110000"

$SHOW $MPIDO $NN pgrn -ssl2 -skip 10 -i $INFILE -ys 0 -ye 110000 ${snaps} ${ts} -o green_ssl2_110ka.nc

