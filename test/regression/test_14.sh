#!/bin/bash

PISM_PATH=$1
MPIEXEC=$2

set -x

# Test name:
echo "Test #14: verif test E regression: isothermal SIA with sliding."
# The list of files to delete when done.
files="test_14-E-out.txt verify.nc verify.nc~"

rm -f $files

# run test E
OPTS="-test E -y 100 -o_size small -verbose 1 -Mbz 1"
$PISM_PATH/pismv -Mx 21 -My 21 -Mz 3 $OPTS  > test_14-E-out.txt
$PISM_PATH/pismv -Mx 41 -My 41 -Mz 3 $OPTS >> test_14-E-out.txt

# compare results
diff test_14-E-out.txt -  <<END-OF-OUTPUT
NUMERICAL ERRORS evaluated at final time (relative to exact solution):
geometry  :    prcntVOL        maxH         avH   relmaxETA
               0.848651  368.057068   32.711326    0.041252
base vels :  maxvector   avvector  prcntavvec     maxub     maxvb
                4.1538    0.22008     0.43350    4.0511    2.0747
NUM ERRORS DONE
NUMERICAL ERRORS evaluated at final time (relative to exact solution):
geometry  :    prcntVOL        maxH         avH   relmaxETA
               0.352172  526.115752   26.312808    0.040403
base vels :  maxvector   avvector  prcntavvec     maxub     maxvb
                1.3782    0.10231     0.20151    1.3573    1.0278
NUM ERRORS DONE
END-OF-OUTPUT

if [ $? != 0 ];
then
    exit 1
fi

rm -f $files; exit 0

