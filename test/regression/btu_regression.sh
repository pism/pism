#!/bin/bash

# PISMBedThermalUnit regression test using z=0 heat flux output of verification test K

PISM_PATH=$1
MPIEXEC=$2
PISM_SOURCE_DIR=$3

# List of files to remove when done:
files="foo-btu.nc foo-btu.nc~ btu_test_out.txt"

rm -f $files

set -e -x

OPTS="-verbose 1 -ys 0.0 -ye 1.0 -dt 0.1 -o foo-btu.nc"

# do stuff: test with a litho_temp state variable, and without
$PISM_PATH/btutest -Mbz 11 -Lbz 1000 $OPTS > btu_test_out.txt
$PISM_PATH/btutest -Mbz 1 $OPTS >> btu_test_out.txt

set +e

# Check results:
diff btu_test_out.txt -  <<END-OF-OUTPUT
NUMERICAL ERRORS in upward heat flux at z=0 relative to exact solution:
bheatflx0  :       max    prcntmax          av
             0.0034644   11.2927879    0.0034644
NUM ERRORS DONE
NUMERICAL ERRORS in upward heat flux at z=0 relative to exact solution:
bheatflx0  :       max    prcntmax          av
             0.0113218   36.9052106    0.0113218
NUM ERRORS DONE
END-OF-OUTPUT

if [ $? != 0 ];
then
    exit 1
fi

rm -f $files; exit 0

