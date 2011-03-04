#!/bin/bash

PISM_PATH=$1
MPIEXEC=$2

# Test name:
test="Test #15: verif test E regression: isothermal SIA with sliding."
# The list of files to delete when done.
files="test_15-E-out.txt verify.nc verify.nc~"

rm -f $files

# run test E
OPTS="-test E -y 1000 -o_size small -verbose 1 -Mbz 1"
$PISM_PATH/pismv -Mx 21 -My 21 -Mz 21 $OPTS  > test_15-E-out.txt
$PISM_PATH/pismv -Mx 41 -My 41 -Mz 41 $OPTS >> test_15-E-out.txt

# compare results
diff test_15-E-out.txt - > /dev/null <<END-OF-OUTPUT
NUMERICAL ERRORS evaluated at final time (relative to exact solution):
geometry  :    prcntVOL        maxH         avH   relmaxETA
               6.571588  747.089432   89.407623    0.109291
base vels :  maxvector   avvector  prcntavvec     maxub     maxvb
                5.9483    0.29471     0.58051    5.9475    2.7752
NUM ERRORS DONE
NUMERICAL ERRORS evaluated at final time (relative to exact solution):
geometry  :    prcntVOL        maxH         avH   relmaxETA
               3.609672  720.379860   51.405684    0.058540
base vels :  maxvector   avvector  prcntavvec     maxub     maxvb
                1.8179    0.14171     0.27914    1.6160    1.2159
NUM ERRORS DONE
END-OF-OUTPUT

if [ $? != 0 ];
then
    exit 1
fi

rm -f $files; exit 0

