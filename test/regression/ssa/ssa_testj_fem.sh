#!/bin/bash

# SSAFEM verification test J regression test

PISM_PATH=$1
MPIEXEC=$2
MPIEXEC_COMMAND="$MPIEXEC -n 2"
PISM_SOURCE_DIR=$3
EXT=""
if [ $# -ge 4 ] && [ "$4" == "-python" ]
then
  PYTHONEXEC=$5
  MPIEXEC_COMMAND="$MPIEXEC_COMMAND $PYTHONEXEC"
  PYTHONPATH=${PISM_PATH}:${PYTHONPATH}
  PISM_PATH=${PISM_SOURCE_DIR}/examples/python/ssa_tests
  EXT=".py"
fi

# List of files to remove when done:
files="foo-fem-j.nc foo-fem-j.nc~ test-J-out-fem.txt"

rm -f $files

set -e

OPTS="-verbose 1 -ssa_method fem -o foo-fem-j.nc -ksp_type cg"

# do stuff
$MPIEXEC_COMMAND $PISM_PATH/ssa_testj${EXT} -Mx 61 -My 61 $OPTS > test-J-out-fem.txt
$MPIEXEC_COMMAND $PISM_PATH/ssa_testj${EXT} -Mx 121 -My 121 $OPTS >> test-J-out-fem.txt

set +e

# Check results:
diff test-J-out-fem.txt -  <<END-OF-OUTPUT
NUMERICAL ERRORS in velocity relative to exact solution:
velocity  :  maxvector   prcntavvec      maxu      maxv       avu       avv
                0.1709      0.05526    0.1679    0.0540    0.0934    0.0285
NUM ERRORS DONE
NUMERICAL ERRORS in velocity relative to exact solution:
velocity  :  maxvector   prcntavvec      maxu      maxv       avu       avv
                0.0422      0.01401    0.0415    0.0134    0.0237    0.0072
NUM ERRORS DONE
END-OF-OUTPUT

if [ $? != 0 ];
then
    exit 1
fi

rm -f $files; exit 0
