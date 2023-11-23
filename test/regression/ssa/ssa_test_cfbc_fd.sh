#!/bin/bash

# SSAFD verification test V (van der Veen) regression test

PISM_PATH=$1
MPIEXEC=$2
MPIEXEC_COMMAND="$MPIEXEC -n 1"
PISM_SOURCE_DIR=$3
EXT=""
if [ $# -ge 4 ] && [ "$4" == "-python" ]
then
  PYTHONEXEC=$5
  MPIEXEC_COMMAND="$MPIEXEC_COMMAND $PYTHONEXEC"
  PYTHONPATH=${PISM_PATH}/site-packages:${PYTHONPATH}
  PISM_PATH=${PISM_SOURCE_DIR}/examples/python/ssa_tests
  EXT=".py"
fi

output=`mktemp pism-testv-XXXX` || exit 1

set -e
set -x

OPTS="
-verbose 1
-o_size none
-My 3
-ssafd_ksp_type richardson
-ssafd_pc_type lu
-stress_balance.ssa.epsilon 0"

# do stuff
$MPIEXEC_COMMAND $PISM_PATH/ssa_test_cfbc${EXT} -Mx 201 $OPTS > ${output}
$MPIEXEC_COMMAND $PISM_PATH/ssa_test_cfbc${EXT} -Mx 401 $OPTS >> ${output}

set +e

# Check results:
diff ${output} -  <<END-OF-OUTPUT
NUMERICAL ERRORS in velocity relative to exact solution:
velocity  :  maxvector   prcntavvec      maxu      maxv       avu       avv
                0.1347      0.00684    0.1347    0.0000    0.0667    0.0000
NUM ERRORS DONE
NUMERICAL ERRORS in velocity relative to exact solution:
velocity  :  maxvector   prcntavvec      maxu      maxv       avu       avv
                0.0436      0.00359    0.0436    0.0000    0.0350    0.0000
NUM ERRORS DONE
END-OF-OUTPUT

if [ $? != 0 ];
then
  cat ${output}
  exit 1
fi

exit 0
