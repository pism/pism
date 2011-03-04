#!/bin/bash

PISM_PATH=$1
MPIEXEC=$2

# Test name:
test="Test #19: verif test I regression: nonlinear SSA with plastic bed."
# The list of files to delete when done.
files="test-I-out.txt verify.nc verify.nc~"

rm -f $files

# run test I
OPTS="-test I -Mx 5 -Mz 31 -Mbz 1 -ssa_rtol 5e-07 -ksp_rtol 1e-12 -verbose 1 -o_size small"
$PISM_PATH/pismv -My 49  $OPTS  > test-I-out.txt
$PISM_PATH/pismv -My 193 $OPTS >> test-I-out.txt

# compare results
diff test-I-out.txt - > /dev/null <<END-OF-OUTPUT
NUMERICAL ERRORS in velocity relative to exact solution:
velocity  :  maxvector   prcntavvec      maxu       avu
               23.2923      0.98759   23.2923    7.6788
NUM ERRORS DONE
NUMERICAL ERRORS in velocity relative to exact solution:
velocity  :  maxvector   prcntavvec      maxu       avu
                1.3241      0.05675    1.3241    0.4413
NUM ERRORS DONE
END-OF-OUTPUT

if [ $? != 0 ];
then
    exit 1
fi

rm -f $files; exit 0

