#!/bin/bash

# SSAFEM verification test J regression test

PISM_PATH=$1
MPIEXEC=$2
PISM_SOURCE_DIR=$3

# List of files to remove when done:
files="foo.nc foo.nc~ test-J-out.txt"

rm -f $files

set -e

OPTS="-verbose 1 -ssa_method fem -o foo.nc"

# do stuff
$PISM_PATH/ssa_testj -Mx 61 -My 61 $OPTS > test-J-out.txt
$PISM_PATH/ssa_testj -Mx 121 -My 121 $OPTS >> test-J-out.txt

set +e

# Check results:
diff test-J-out.txt - > /dev/null <<END-OF-OUTPUT
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
