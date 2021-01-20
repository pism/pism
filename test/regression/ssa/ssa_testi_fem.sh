#!/bin/bash

# SSAFEM verification test I regression test

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

output=`mktemp pism-ssa-test-i.XXXX` || exit 1

set -e

OPTS="-verbose 1 -ssa_method fem -o_size none -Mx 5 -ksp_type cg"

# do stuff
$MPIEXEC_COMMAND $PISM_PATH/ssa_testi${EXT} -My 61 $OPTS > ${output}
$MPIEXEC_COMMAND $PISM_PATH/ssa_testi${EXT} -My 121 $OPTS >> ${output}

set +e

# Check results:
diff ${output} -  <<END-OF-OUTPUT
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
  cat ${output}
  exit 1
fi

exit 0
