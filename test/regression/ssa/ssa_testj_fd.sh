#!/bin/bash

# SSAFD verification test J regression test

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
files="foo-fd-j.nc foo-fd-j.nc~ test-J-out-fd.txt"

rm -f $files

set -e

OPTS="-verbose 1 -ssa_method fd -o foo-fd-j.nc"

# do stuff
$MPIEXEC_COMMAND $PISM_PATH/ssa_testj${EXT} -Mx 61 -My 61 $OPTS > test-J-out-fd.txt
$MPIEXEC_COMMAND $PISM_PATH/ssa_testj${EXT} -Mx 121 -My 121 $OPTS >> test-J-out-fd.txt

set +e

# Check results:
diff test-J-out-fd.txt -  <<END-OF-OUTPUT
NUMERICAL ERRORS in velocity relative to exact solution:
velocity  :  maxvector   prcntavvec      maxu      maxv       avu       avv
                0.2699      0.10496    0.1903    0.2128    0.0949    0.1578
NUM ERRORS DONE
NUMERICAL ERRORS in velocity relative to exact solution:
velocity  :  maxvector   prcntavvec      maxu      maxv       avu       avv
                0.1410      0.05888    0.1067    0.1011    0.0631    0.0832
NUM ERRORS DONE
END-OF-OUTPUT

if [ $? != 0 ];
then
    exit 1
fi

rm -f $files; exit 0
