#!/bin/bash

#   Runs the CCL3 EISMINT-Greenland experiment.  Requires GRIP ice core 
# temperature data in grip_dT.nc and SPECMAP sea bed core (stack) sea level
# data in specmap_dSL.nc.  Saves state every 50000 model years.
#   Also requires green_SSL2_110k.nc, the final state of the EISMINT-Greenland 
# SSL2 run.  See examples/eisgreen/ssl2.sh.
#   See the PISM User's Manual for the modeling choices in "pgrn -ccl3".

#   A recommended way to run with 2 processors is "./ccl3.sh 2 >& ccl3.txt &"
# which saves a transcript in ccl3.txt.

NN=2  # default number of processors
if [ $# -gt 0 ] ; then  # if user says "ccl3.sh 8" then NN = 8
  NN="$1"
fi
set -e  # exit on error

mpiexec -n $NN pgrn -ccl3 -if green_SSL2_110k.nc \
   -dTforcing grip_dT.nc -dSLforcing specmap_dSL.nc \
   -ys -249900 -y 49900 -o green_CCL3_y-200000
   
mpiexec -n $NN pgrn -ccl3 -if green_CCL3_y-200000.nc \
   -dTforcing grip_dT.nc -dSLforcing specmap_dSL.nc \
   -y 50000 -o green_CCL3_y-150000
   
mpiexec -n $NN pgrn -ccl3 -if green_CCL3_y-150000.nc \
   -dTforcing grip_dT.nc -dSLforcing specmap_dSL.nc \
   -y 50000 -o green_CCL3_y-100000
   
mpiexec -n $NN pgrn -ccl3 -if green_CCL3_y-100000.nc \
   -dTforcing grip_dT.nc -dSLforcing specmap_dSL.nc \
   -y 50000 -o green_CCL3_y-50000

mpiexec -n $NN pgrn -ccl3 -if green_CCL3_y-50000.nc \
   -dTforcing grip_dT.nc -dSLforcing specmap_dSL.nc \
   -ye 0 -o green_CCL3_y0

