#!/bin/bash

PISM_PATH=$1
MPIEXEC=$2

# Test name:
test="Test #15: verif test E regression: isothermal SIA with sliding."
# The list of files to delete when done.
files="test_15-E-out.txt verify.nc verify.nc~"

rm -f $files

# run test E
OPTS="-test E -y 100 -o_size small -verbose 1 -Mbz 1"
$PISM_PATH/pismv -Mx 21 -My 21 -Mz 3 $OPTS  > test_15-E-out.txt
$PISM_PATH/pismv -Mx 41 -My 41 -Mz 3 $OPTS >> test_15-E-out.txt

# compare results
diff test_15-E-out.txt -  <<END-OF-OUTPUT
NUMERICAL ERRORS evaluated at final time (relative to exact solution):
geometry  :    prcntVOL        maxH         avH   relmaxETA
               0.827276  368.801098   32.698040    0.040785
base vels :  maxvector   avvector  prcntavvec     maxub     maxvb
                4.2915    0.22484     0.44288    4.1961    2.2109
NUM ERRORS DONE
NUMERICAL ERRORS evaluated at final time (relative to exact solution):
geometry  :    prcntVOL        maxH         avH   relmaxETA
               0.329772  525.609347   26.226880    0.040169
base vels :  maxvector   avvector  prcntavvec     maxub     maxvb
                1.3480    0.10223     0.20137    1.3375    1.0205
NUM ERRORS DONE
END-OF-OUTPUT

if [ $? != 0 ];
then
    exit 1
fi

rm -f $files; exit 0

