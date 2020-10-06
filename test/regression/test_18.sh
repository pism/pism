#!/bin/bash

PISM_PATH=$1
MPIEXEC=$2

# Test name:
echo "Test #18: verif test K regression: cold ice method, bedrock thermal layer."
# The list of files to delete when done.
files="test-K-out.txt"

rm -f $files

# run test K
OPTS="-test K -Mx 4 -My 4 -y 13000.0 -Lbz 1000 -z_spacing equal -verbose 1 -o_size none"
$PISM_PATH/pismv -Mz 41 -Mbz 11 -max_dt 60.0 $OPTS  > test-K-out.txt
$PISM_PATH/pismv -Mz 81 -Mbz 21 -max_dt 30.0 $OPTS >> test-K-out.txt

# compare results
diff test-K-out.txt -  <<END-OF-OUTPUT
NUMERICAL ERRORS evaluated at final time (relative to exact solution):
temp      :        maxT         avT       maxTb        avTb
               0.005298    0.002065    0.001899    0.001297
NUM ERRORS DONE
NUMERICAL ERRORS evaluated at final time (relative to exact solution):
temp      :        maxT         avT       maxTb        avTb
               0.001537    0.000618    0.001074    0.000664
NUM ERRORS DONE
END-OF-OUTPUT

if [ $? != 0 ];
then
    exit 1
fi

rm -f $files; exit 0

