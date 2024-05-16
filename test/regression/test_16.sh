#!/bin/bash

PISM_PATH=$1
MPIEXEC=$2

# Test name:
echo "Test #16: verif test L regression: isothermal SIA with non-flat bed."

output=`mktemp pism-test-L.XXXX` || exit 1

set -x
set -e

# run test L
OPTS="-test L -Mbz 1 -Mz 31 -y 1000years -max_dt 60years -o_size none -verbose 1"
$MPIEXEC -n 2 $PISM_PATH/pismr -Mx 21 -My 21 $OPTS   > ${output}
$MPIEXEC -n 2 $PISM_PATH/pismr -Mx 31 -My 31 $OPTS  >> ${output}

# compare results
diff ${output} -  <<END-OF-OUTPUT
NUMERICAL ERRORS evaluated at final time (relative to exact solution):
geometry  :    prcntVOL        maxH         avH   relmaxETA
               0.283121  168.847973    8.958281    0.004678
NUM ERRORS DONE
NUMERICAL ERRORS evaluated at final time (relative to exact solution):
geometry  :    prcntVOL        maxH         avH   relmaxETA
               0.030774  143.588110    3.831173    0.002931
NUM ERRORS DONE
END-OF-OUTPUT

if [ $? != 0 ];
then
  cat ${output}
  exit 1
fi

exit 0
