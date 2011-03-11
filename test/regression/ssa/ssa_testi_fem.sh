#!/bin/bash

# SSAFEM verification test I regression test

PISM_PATH=$1
MPIEXEC=$2
PISM_SOURCE_DIR=$3

# List of files to remove when done:
files="foo.nc foo.nc~ test-I-out.txt"

rm -f $files

set -e

OPTS="-verbose 1 -ssa_method fem -o foo.nc -Mx 5"

# do stuff
$PISM_PATH/ssa_testi -My 61 $OPTS > test-I-out.txt
$PISM_PATH/ssa_testi -My 121 $OPTS >> test-I-out.txt

set +e

# Check results:
diff test-I-out.txt - > /dev/null <<END-OF-OUTPUT
NUMERICAL ERRORS in velocity relative to exact solution:
velocity  :  maxvector   prcntavvec      maxu      maxv       avu       avv
               16.2024      0.14888   16.2024    0.7544    1.1522    0.0513
NUM ERRORS DONE
NUMERICAL ERRORS in velocity relative to exact solution:
velocity  :  maxvector   prcntavvec      maxu      maxv       avu       avv
                4.2045      0.03669    4.2045    0.1967    0.2838    0.0134
NUM ERRORS DONE
END-OF-OUTPUT

if [ $? != 0 ];
then
    exit 1
fi

rm -f $files; exit 0
