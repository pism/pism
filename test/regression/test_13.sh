#!/bin/bash

PISM_PATH=$1
MPIEXEC=$2

# Test name:
test="Test #13: temperature continuity at ice-bed interface (cold case)."
# The list of files to delete when done.
files="foo.nc temp.nc litho_temp.nc"

rm -f $files

# generate sample state with cold base, but which has seen some strain-heating
$MPIEXEC -n 2 $PISM_PATH/pisms -eisII A -Mx 31 -My 31 -Mz 31 -Lbz 1000 -Mbz 11 \
    -y 6e3 -temp_pa -o foo.nc

# extract only what is needed for comparison
ncks -v temp -d z,0m foo.nc temp.nc
ncks -v litho_temp -d zb,-0.000001m foo.nc litho_temp.nc
# neither of these seems to work:
#run ncks -v litho_temp -d zb,0m foo.nc litho_temp.nc
#run ncks -v litho_temp -d zb,-0m foo.nc litho_temp.nc

# compare
ncrename -O -v litho_temp,temp litho_temp.nc
$PISM_PATH/nccmp.py -v temp temp.nc litho_temp.nc

if [ $? != 0 ];
then
    exit 1
fi

rm -f $files; exit 0

