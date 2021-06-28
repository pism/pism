#!/bin/bash

# PISMBedThermalUnit regression test using z=0 heat flux output of verification test K

PISM_PATH=$1
MPIEXEC=$2
PISM_SOURCE_DIR=$3

output=`mktemp pism-btutest.XXXX` || exit 1

set -e -x

OPTS="-verbose 1 -ys 0.0 -ye 1.0 -dt 0.1 -o_size none"

# do stuff: test with a litho_temp state variable, and without
$PISM_PATH/btutest -Mbz 11 -Lbz 1000 $OPTS > ${output}
$PISM_PATH/btutest -Mbz 1 $OPTS >> ${output}

set +e

# Check results:
diff ${output} -  <<END-OF-OUTPUT
NUMERICAL ERRORS in upward heat flux at z=0 relative to exact solution:
bheatflx0  :       max    prcntmax          av
             0.0034644   11.2928460    0.0034644
NUM ERRORS DONE
NUMERICAL ERRORS in upward heat flux at z=0 relative to exact solution:
bheatflx0  :       max    prcntmax          av
             0.0113218   36.9052047    0.0113218
NUM ERRORS DONE
END-OF-OUTPUT

if [ $? != 0 ];
then
  cat ${output}
  exit 1
fi

exit 0
