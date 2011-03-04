#!/bin/bash

PISM_PATH=$1
MPIEXEC=$2

# Test name:
test="Test #20: verif test J regression: linearized SSA, floating."
# The list of files to delete when done.
files="test-J-out.txt verify.nc verify.nc~"

rm -f $files

# run test J
OPTS="-test J -Mz 11 -Mbz 1 -pc_type asm -sub_pc_type lu -ksp_rtol 1e-12 -verbose 1 -o_size small"
$PISM_PATH/pismv -Mx 49 -My 49 $OPTS  > test-J-out.txt
$PISM_PATH/pismv -Mx 98 -My 98 $OPTS >> test-J-out.txt

# compare results
diff test-J-out.txt - > /dev/null <<END-OF-OUTPUT
NUMERICAL ERRORS in velocity relative to exact solution:
velocity  :  maxvector   prcntavvec      maxu      maxv       avu       avv
                0.2526      0.08339    0.2421    0.0725    0.1432    0.0385
NUM ERRORS DONE
NUMERICAL ERRORS in velocity relative to exact solution:
velocity  :  maxvector   prcntavvec      maxu      maxv       avu       avv
                0.0604      0.02081    0.0578    0.0174    0.0357    0.0096
NUM ERRORS DONE
END-OF-OUTPUT

if [ $? != 0 ];
then
    exit 1
fi

rm -f $files; exit 0

