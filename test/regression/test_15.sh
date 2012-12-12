#!/bin/bash

PISM_PATH=$1
MPIEXEC=$2

# Test name:
echo "Test #15: verif test C regression: isothermal SIA w. time-dependent SMB."
# The list of files to delete when done.
files="test_15-C-out.txt"

rm -f $files

# run test C
OPTS="-test C -Mbz 1 -Mz 31 -y 5000 -o_size none -verbose 1"
$PISM_PATH/pismv -Mx 31 -My 31 $OPTS  > test_15-C-out.txt
$PISM_PATH/pismv -Mx 41 -My 41 $OPTS >> test_15-C-out.txt

# compare results
diff test_15-C-out.txt -  <<END-OF-OUTPUT
NUMERICAL ERRORS evaluated at final time (relative to exact solution):
geometry  :    prcntVOL        maxH         avH   relmaxETA
              80.824124  503.131175    3.114691    0.820828
NUM ERRORS DONE
NUMERICAL ERRORS evaluated at final time (relative to exact solution):
geometry  :    prcntVOL        maxH         avH   relmaxETA
              32.783573  193.022555    1.330304    0.405692
NUM ERRORS DONE
END-OF-OUTPUT

if [ $? != 0 ];
then
    exit 1
fi

rm -f $files; exit 0

