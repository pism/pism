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

OPTS="-verbose 1 -ssa_method fd -o_size none -ssafd_pc_type gamg"

# do stuff
$MPIEXEC_COMMAND $PISM_PATH/ssa_testj${EXT} -Mx 61 -My 61 $OPTS > ${output}
$MPIEXEC_COMMAND $PISM_PATH/ssa_testj${EXT} -Mx 121 -My 121 $OPTS >> ${output}

set +e

# Check results:
diff ${output} -  <<END-OF-OUTPUT
NUMERICAL ERRORS in velocity relative to exact solution:
velocity  :  maxvector   prcntavvec      maxu      maxv       avu       avv
                0.1565      0.05377    0.1502    0.0455    0.0923    0.0248
NUM ERRORS DONE
NUMERICAL ERRORS in velocity relative to exact solution:
velocity  :  maxvector   prcntavvec      maxu      maxv       avu       avv
                0.0426      0.01377    0.0407    0.0132    0.0236    0.0064
NUM ERRORS DONE
END-OF-OUTPUT

if [ $? != 0 ];
then
  cat ${output}
  exit 1
fi

exit 0
