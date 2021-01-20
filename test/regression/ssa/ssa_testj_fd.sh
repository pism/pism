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
  PYTHONPATH=${PISM_PATH}/site-packages:${PYTHONPATH}
  PISM_PATH=${PISM_SOURCE_DIR}/examples/python/ssa_tests
  EXT=".py"
fi

output=`mktemp pism-ssa-test-j.XXXX` || exit 1

set -e

OPTS="-verbose 1 -ssa_method fd -o_size none"

# do stuff
$MPIEXEC_COMMAND $PISM_PATH/ssa_testj${EXT} -Mx 61 -My 61 $OPTS > ${output}
$MPIEXEC_COMMAND $PISM_PATH/ssa_testj${EXT} -Mx 121 -My 121 $OPTS >> ${output}

set +e

# Check results:
diff ${output} -  <<END-OF-OUTPUT
NUMERICAL ERRORS in velocity relative to exact solution:
velocity  :  maxvector   prcntavvec      maxu      maxv       avu       avv
                0.2888      0.08782    0.2584    0.1457    0.1158    0.0933
NUM ERRORS DONE
NUMERICAL ERRORS in velocity relative to exact solution:
velocity  :  maxvector   prcntavvec      maxu      maxv       avu       avv
                0.1689      0.07538    0.1181    0.1339    0.0719    0.1139
NUM ERRORS DONE
END-OF-OUTPUT

if [ $? != 0 ];
then
  cat ${output}
  exit 1
fi

exit 0
