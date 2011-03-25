#!/bin/bash

PISM_PATH=$1
MPIEXEC=$2

# Test name:
echo "Test #21: verif test K regression: cold ice method, bedrock thermal layer."
# The list of files to delete when done.
files="test-K-out.txt verify.nc verify.nc~"

rm -f $files

# run test K
OPTS="-test K -Mx 4 -My 4 -y 13000.0 -Lbz 1000 -verbose 1 -o_size small"
$PISM_PATH/pismv -Mz 41 -Mbz 11 $OPTS  > test-K-out.txt
$PISM_PATH/pismv -Mz 81 -Mbz 21 $OPTS >> test-K-out.txt

# compare results
diff test-K-out.txt -  <<END-OF-OUTPUT
NUMERICAL ERRORS evaluated at final time (relative to exact solution):
temp      :        maxT         avT       maxTb        avTb
               0.436647    0.165911    0.878458    0.633755
NUM ERRORS DONE
NUMERICAL ERRORS evaluated at final time (relative to exact solution):
temp      :        maxT         avT       maxTb        avTb
               0.128707    0.037709    0.290305    0.196936
NUM ERRORS DONE
END-OF-OUTPUT

if [ $? != 0 ];
then
    exit 1
fi

rm -f $files; exit 0

