#!/bin/bash

# SSAFD verification test V (van der Veen) regression test

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
files="foo-V.nc foo-V.nc~ test-V-out.txt"

rm -f $files

set -e
set -x

OPTS="-verbose 1 -o foo-V.nc -My 5"

# do stuff
$MPIEXEC_COMMAND $PISM_PATH/ssa_test_cfbc${EXT} -Mx 201 $OPTS > test-V-out.txt
$MPIEXEC_COMMAND $PISM_PATH/ssa_test_cfbc${EXT} -Mx 401 $OPTS >> test-V-out.txt

set +e

# Check results:
diff test-V-out.txt -  <<END-OF-OUTPUT
NUMERICAL ERRORS in velocity relative to exact solution:
velocity  :  maxvector   prcntavvec      maxu      maxv       avu       avv
                1.3840      0.08724    1.3223    0.5238    0.7898    0.3128
NUM ERRORS DONE
NUMERICAL ERRORS in velocity relative to exact solution:
velocity  :  maxvector   prcntavvec      maxu      maxv       avu       avv
                0.5034      0.02569    0.4890    0.1395    0.2315    0.0829
NUM ERRORS DONE
END-OF-OUTPUT

if [ $? != 0 ];
then
    exit 1
fi

rm -f $files; exit 0
