#!/bin/bash

# SSAFD verification test I regression test

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

output=`mktemp pism-test-i.XXXX` || exit 1

set -e
set -x

OPTS="-verbose 1 -ssa_method fd -o_size none -ssafd_picard_rtol 5e-07 -ssafd_ksp_rtol 1e-12 -Mx 5"

# do stuff
$MPIEXEC_COMMAND $PISM_PATH/ssa_testi${EXT} -My 61 $OPTS > ${output}
$MPIEXEC_COMMAND $PISM_PATH/ssa_testi${EXT} -My 121 $OPTS >> ${output}

set +e

# Check results:
diff ${output} -  <<END-OF-OUTPUT
NUMERICAL ERRORS in velocity relative to exact solution:
velocity  :  maxvector   prcntavvec      maxu      maxv       avu       avv
                4.7417      0.05219    4.7417    0.1976    0.4041    0.0087
NUM ERRORS DONE
NUMERICAL ERRORS in velocity relative to exact solution:
velocity  :  maxvector   prcntavvec      maxu      maxv       avu       avv
                1.3907      0.01351    1.3907    0.0385    0.1050    0.0018
NUM ERRORS DONE
END-OF-OUTPUT

if [ $? != 0 ];
then
  cat ${output}
  exit 1
fi

exit 0
