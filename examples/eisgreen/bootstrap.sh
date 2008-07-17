#!/bin/bash

#   Creates a credible initial state for the EISMINT-Greenland experiments.
# Starts by applying scripts eis_green.py and eis_core.py to data downloaded from
#       http://homepages.vub.ac.be/~phuybrec/eismint/
# Requires
#       grid20-EISMINT, suaq20-EISMINT, specmap.017, sum89-92-ss09-50yr.stp
# from that site.  (Uncomment the "wget" lines below to download.)

#   Recommended way to run with 2 processors is "./bootstrap.sh 2 >& bootstrap.txt &"
# which gives a transcript in bootstrap.txt.

NN=2  # default number of processors
if [ $# -gt 0 ] ; then  # if user says "ccl3.sh 8" then NN = 8
  NN="$1"
fi
set -e  # exit on error

# Uncomment these lines to download the files using wget:
#for fname in "grid20-EISMINT" "suaq20-EISMINT" "specmap.017" "sum89-92-ss09-50yr.stp"
#do
#  wget http://homepages.vub.ac.be/~phuybrec/eismint/$fname
#done

./eis_green.py  # creates  eis_green20.nc

fill_missing.py -f eis_green20.nc -v topg -o eis_green_smoothed.nc

./eis_core.py  # creates  grip_dT.nc  and  specmap_dSL.nc

pgrn -bif eis_green_smoothed.nc -Mx 141 -My 83 -Lz 4000 -Mz 51 -quadZ \
       -ocean_kill -tempskip 0 -y 1 -o green20km_y1

mpiexec -n $NN pgrn -if green20km_y1.nc -no_mass -y 10000 -o green20km_Tsteady

# uncomment these lines if 100k model year temperature relaxation is desired before SSL2:
#mv green20km_Tsteady.nc green20km_10k.nc
#mpiexec -n $NN pgrn -if green20km_10k.nc -no_mass -y 90000 -o green20km_Tsteady

