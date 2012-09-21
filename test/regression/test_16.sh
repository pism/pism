#!/bin/bash

PISM_PATH=$1
MPIEXEC=$2

# Test name:
echo "Test #16: verif test L regression: isothermal SIA with non-flat bed."
# The list of files to delete when done.
files="test_16-L-out.txt"

rm -f $files

# run test L
OPTS="-test L -Mbz 1 -Mz 31 -y 1000 -o_size none -verbose 1"
$PISM_PATH/pismv -Mx 31 -My 31 $OPTS   > test_16-L-out.txt
$PISM_PATH/pismv -Mx 41 -My 41 $OPTS  >> test_16-L-out.txt

# compare results
diff test_16-L-out.txt -  <<END-OF-OUTPUT
NUMERICAL ERRORS evaluated at final time (relative to exact solution):
geometry  :    prcntVOL        maxH         avH   relmaxETA
               0.029956  143.583244    3.835807    0.002941
NUM ERRORS DONE
NUMERICAL ERRORS evaluated at final time (relative to exact solution):
geometry  :    prcntVOL        maxH         avH   relmaxETA
               0.080578  154.548155    3.830053    0.002503
NUM ERRORS DONE
END-OF-OUTPUT

if [ $? != 0 ];
then
    exit 1
fi

rm -f $files; exit 0

