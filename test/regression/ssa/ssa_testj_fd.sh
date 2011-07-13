#!/bin/bash

# SSAFD verification test J regression test

PISM_PATH=$1
MPIEXEC=$2
PISM_SOURCE_DIR=$3
EXT=""
if [ $# -ge 4 ] && [ "$4" == "-python" ]
then
  PISM_PATH=${PISM_SOURCE_DIR}/examples/python/ssa_tests
  EXT=".py"
fi

# List of files to remove when done:
files="foo.nc foo.nc~ test-J-out.txt"

rm -f $files

set -e

OPTS="-verbose 1 -ssa_method fd -o foo.nc"

# do stuff
$PISM_PATH/ssa_testj${EXT} -Mx 61 -My 61 $OPTS > test-J-out.txt
$PISM_PATH/ssa_testj${EXT} -Mx 121 -My 121 $OPTS >> test-J-out.txt

set +e

# Check results:
diff test-J-out.txt -  <<END-OF-OUTPUT
NUMERICAL ERRORS in velocity relative to exact solution:
velocity  :  maxvector   prcntavvec      maxu      maxv       avu       avv
                0.1825      0.05477    0.1750    0.0615    0.0931    0.0265
NUM ERRORS DONE
NUMERICAL ERRORS in velocity relative to exact solution:
velocity  :  maxvector   prcntavvec      maxu      maxv       avu       avv
                0.1331      0.04919    0.1314    0.0397    0.0850    0.0245
NUM ERRORS DONE
END-OF-OUTPUT

if [ $? != 0 ];
then
    exit 1
fi

rm -f $files; exit 0
