#!/bin/bash

# SSAFD verification test V (van der Veen) regression test

PISM_PATH=$1
MPIEXEC=$2
MPIEXEC_COMMAND="$MPIEXEC -n 1"
PISM_SOURCE_DIR=$3
EXT=""
if [ $# -ge 4 ] && [ "$4" == "-python" ]
then
  PYTHONEXEC=$5
  MPIEXEC_COMMAND="$MPIEXEC_COMMAND $PYTHONEXEC"
  PYTHONPATH=${PISM_PATH}/site-packages:${PYTHONPATH}
  PISM_PATH=${PISM_SOURCE_DIR}/examples/python/ssa_tests
  EXT=".py"
fi

# List of files to remove when done:
files="foo-V-fem.nc foo-V-fem.nc~ test-v-out-fem.txt"

rm -f $files

set -e
set -x

OPTS="-verbose 1 -o foo-V-fem.nc -My 3 -ssa_method fem"

# do stuff
$MPIEXEC_COMMAND $PISM_PATH/ssa_test_cfbc${EXT} -Mx 201 $OPTS > test-v-out-fem.txt
$MPIEXEC_COMMAND $PISM_PATH/ssa_test_cfbc${EXT} -Mx 401 $OPTS >> test-v-out-fem.txt

set +e

# Check results:
diff test-v-out-fem.txt -  <<END-OF-OUTPUT
NUMERICAL ERRORS in velocity relative to exact solution:
velocity  :  maxvector   prcntavvec      maxu      maxv       avu       avv
                1.1792      0.11288    1.1792    0.0000    1.0998    0.0000
NUM ERRORS DONE
NUMERICAL ERRORS in velocity relative to exact solution:
velocity  :  maxvector   prcntavvec      maxu      maxv       avu       avv
                0.4323      0.03780    0.4323    0.0000    0.3687    0.0000
NUM ERRORS DONE
END-OF-OUTPUT

if [ $? != 0 ];
then
    exit 1
fi

rm -f $files; exit 0
