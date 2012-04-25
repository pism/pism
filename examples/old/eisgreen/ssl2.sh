#!/bin/bash

# Runs the SSL2 EISMINT-Greenland experiment, starting from the output file
# produced by bootstrap.sh.  Saves map-plane diagnostics every 1000 model years,
# and saves a timeseries of volume every year.  See the PISM User's Manual.
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

exs="-extra_file ex_ssl2.nc -extra_times 1000:1000:100000 -extra_vars diffusivity,temppabase,csurf,hardav,bmelt,tempicethk_basal,mask,dHdt,thk,topg,usurf"

ts="-ts_file vol_ssl2.nc -ts_times 0:yearly:110000"

PGRN="pismr -e 3 -ocean_kill -atmosphere eismint_greenland -surface pdd"

$SHOW $MPIDO $NN $PGRN -skip 20 -i $INFILE -ys 0 -ye 110000 ${exs} ${ts} -o green_ssl2_110ka.nc

