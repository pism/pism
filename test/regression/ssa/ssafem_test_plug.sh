#!/bin/bash

# SSAFEM plug flow regression test

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
files="foo-fem-plug.nc foo-fem-plug.nc~ test-out-plug.txt"

rm -f $files

set -e
set -x

OPTS="-verbose 1 -ssa_method fem -o foo-fem-plug.nc -ksp_type cg"

# do stuff
$MPIEXEC_COMMAND $PISM_PATH/ssa_test_plug${EXT} -Mx 22 -My 31 $OPTS > test-out-plug.txt
$MPIEXEC_COMMAND $PISM_PATH/ssa_test_plug${EXT} -Mx 61 -My 61 $OPTS >> test-out-plug.txt

set +e

# Check results:
diff test-out-plug.txt -  <<END-OF-OUTPUT
NUMERICAL ERRORS in velocity relative to exact solution:
velocity  :  maxvector   prcntavvec      maxu      maxv       avu       avv
                0.2024      0.00559    0.2024    0.0325    0.0765    0.0069
NUM ERRORS DONE
NUMERICAL ERRORS in velocity relative to exact solution:
velocity  :  maxvector   prcntavvec      maxu      maxv       avu       avv
                0.0510      0.00147    0.0510    0.0082    0.0201    0.0019
NUM ERRORS DONE
END-OF-OUTPUT

if [ $? != 0 ];
then
    exit 1
fi

rm -f $files; exit 0
