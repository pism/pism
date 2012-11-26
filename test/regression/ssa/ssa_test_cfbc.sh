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
  PYTHONPATH=${PISM_PATH}
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
                3.9224      0.06646    1.6010    3.9218    0.4326    0.3798
NUM ERRORS DONE
NUMERICAL ERRORS in velocity relative to exact solution:
velocity  :  maxvector   prcntavvec      maxu      maxv       avu       avv
                2.2977      0.05555    1.2724    2.2907    0.4731    0.1565
NUM ERRORS DONE
END-OF-OUTPUT

if [ $? != 0 ];
then
    exit 1
fi

rm -f $files; exit 0
