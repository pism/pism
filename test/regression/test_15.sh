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
diff test_15-E-out.txt -  <<END-OF-OUTPUT
NUMERICAL ERRORS evaluated at final time (relative to exact solution):
geometry  :    prcntVOL        maxH         avH   relmaxETA
               6.589296  747.331102   89.665903    0.109347
base vels :  maxvector   avvector  prcntavvec     maxub     maxvb
                6.0136    0.30921     0.60906    6.0130    2.8310
NUM ERRORS DONE
NUMERICAL ERRORS evaluated at final time (relative to exact solution):
geometry  :    prcntVOL        maxH         avH   relmaxETA
               3.627740  720.374602   51.662995    0.058539
base vels :  maxvector   avvector  prcntavvec     maxub     maxvb
                1.9234    0.15022     0.29589    1.6931    1.1956
NUM ERRORS DONE
END-OF-OUTPUT

if [ $? != 0 ];
then
    exit 1
fi

rm -f $files; exit 0

