#!/bin/bash

# SSAFEM linear flow regression test

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

output=`mktemp pism-ssa-test-linear.XXXX` || exit 1

set -e

OPTS="-verbose 1 -ssa_method fem -o_size none -ksp_type cg"

# do stuff
$MPIEXEC_COMMAND $PISM_PATH/pism_ssa_test_linear${EXT} -Mx 61 -My 61 $OPTS > ${output}
$MPIEXEC_COMMAND $PISM_PATH/pism_ssa_test_linear${EXT} -Mx 121 -My 121 $OPTS >> ${output}

set +e

# Check results:
diff ${output} -  <<END-OF-OUTPUT
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
  cat ${output}
  exit 1
fi

exit 0
