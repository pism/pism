#!/bin/bash

# SSAFEM plug flow regression test

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

output=`mktemp pism-ssa-test-plug.XXXX` || exit 1

set -e
set -x

OPTS="-verbose 1 -ssa_method fem -o_size none -ksp_type cg"

# do stuff
$MPIEXEC_COMMAND $PISM_PATH/ssa_test_plug${EXT} -Mx 22 -My 31 $OPTS > ${output}
$MPIEXEC_COMMAND $PISM_PATH/ssa_test_plug${EXT} -Mx 61 -My 61 $OPTS >> ${output}

set +e

# Check results:
diff ${output} -  <<END-OF-OUTPUT
NUMERICAL ERRORS in velocity relative to exact solution:
velocity  :  maxvector   prcntavvec      maxu      maxv       avu       avv
                0.2024      0.00559    0.2024    0.0325    0.0765    0.0069
NUM ERRORS DONE
NUMERICAL ERRORS in velocity relative to exact solution:
velocity  :  maxvector   prcntavvec      maxu      maxv       avu       avv
                0.0510      0.00147    0.0510    0.0082    0.0201    0.0019
NUM ERRORS DONE
END-OF-OUTPUT

if [ $? != 0 ];
then
  cat ${output}
  exit 1
fi

exit 0
