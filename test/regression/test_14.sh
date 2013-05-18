#!/bin/bash

PISM_PATH=$1
MPIEXEC=$2

set -x

# Test name:
echo "Test #14: verif test E regression: isothermal SIA with sliding."
# The list of files to delete when done.
files="test_14-E-out.txt"

rm -f $files

# run test E
OPTS="-test E -y 100 -o_size none -verbose 1 -Mbz 1"
$PISM_PATH/pismv -Mx 21 -My 21 -Mz 3 $OPTS  > test_14-E-out.txt
$PISM_PATH/pismv -Mx 41 -My 41 -Mz 3 $OPTS >> test_14-E-out.txt

# compare results
diff test_14-E-out.txt -  <<END-OF-OUTPUT
NUMERICAL ERRORS evaluated at final time (relative to exact solution):
geometry  :    prcntVOL        maxH         avH   relmaxETA
               0.848937  368.203505   32.724467    0.041279
base vels :  maxvector   avvector  prcntavvec     maxub     maxvb
                4.1525    0.22000     0.43334    4.0495    2.0718
NUM ERRORS DONE
NUMERICAL ERRORS evaluated at final time (relative to exact solution):
geometry  :    prcntVOL        maxH         avH   relmaxETA
               0.352302  526.145853   26.314556    0.040406
base vels :  maxvector   avvector  prcntavvec     maxub     maxvb
                1.3801    0.10244     0.20178    1.3599    1.0328
NUM ERRORS DONE
END-OF-OUTPUT

if [ $? != 0 ];
then
    exit 1
fi

rm -f $files; exit 0

