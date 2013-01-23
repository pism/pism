#!/bin/bash

# SSAFEM linear flow regression test

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
files="foo-fem-linear.nc foo-fem-linear.nc~ test-out-linear.txt"

rm -f $files

set -e

OPTS="-verbose 1 -ssa_method fem -o foo-fem-linear.nc -ksp_type cg"

# do stuff
$MPIEXEC_COMMAND $PISM_PATH/ssa_test_linear${EXT} -Mx 61 -My 61 $OPTS > test-out-linear.txt
$MPIEXEC_COMMAND $PISM_PATH/ssa_test_linear${EXT} -Mx 121 -My 121 $OPTS >> test-out-linear.txt

set +e

# Check results:
diff test-out-linear.txt -  <<END-OF-OUTPUT
NUMERICAL ERRORS in velocity relative to exact solution:
velocity  :  maxvector   prcntavvec      maxu      maxv       avu       avv
                0.4207      0.00319    0.4207    0.0444    0.1884    0.0093
NUM ERRORS DONE
NUMERICAL ERRORS in velocity relative to exact solution:
velocity  :  maxvector   prcntavvec      maxu      maxv       avu       avv
                0.1051      0.00081    0.1051    0.0111    0.0479    0.0024
NUM ERRORS DONE
END-OF-OUTPUT

if [ $? != 0 ];
then
    exit 1
fi

rm -f $files; exit 0
