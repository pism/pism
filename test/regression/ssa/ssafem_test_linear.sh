#!/bin/bash

# SSAFEM linear flow regression test

PISM_PATH=$1
MPIEXEC=$2
PISM_SOURCE_DIR=$3

# List of files to remove when done:
files="foo.nc foo.nc~ test-out.txt"

rm -f $files

set -e

OPTS="-verbose 1 -ssa_method fem -o foo.nc"

# do stuff
$PISM_PATH/ssa_test_linear -Mx 61 -My 61 $OPTS > test-out.txt
$PISM_PATH/ssa_test_linear -Mx 121 -My 121 $OPTS >> test-out.txt

set +e

# Check results:
diff test-out.txt -  <<END-OF-OUTPUT
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
