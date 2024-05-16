#!/bin/bash

PISM_PATH=$1
MPIEXEC=$2

# Test name:
echo "Test #15: verif test C regression: isothermal SIA w. time-dependent SMB."

output=`mktemp pism-test-C.XXXX` || exit 1

# run test C
OPTS="-test C -Mbz 1 -Mz 31 -y 5000years -o_size none -verbose 1 -max_dt 60years"
$MPIEXEC -n 4 $PISM_PATH/pismr -Mx 31 -My 31 $OPTS  > ${output}
$MPIEXEC -n 4 $PISM_PATH/pismr -Mx 41 -My 41 $OPTS >> ${output}

# compare results
diff ${output} -  <<END-OF-OUTPUT
NUMERICAL ERRORS evaluated at final time (relative to exact solution):
geometry  :    prcntVOL        maxH         avH   relmaxETA
              80.824124  503.131175    3.114691    0.820828
NUM ERRORS DONE
NUMERICAL ERRORS evaluated at final time (relative to exact solution):
geometry  :    prcntVOL        maxH         avH   relmaxETA
              32.783573  193.022555    1.330304    0.405692
NUM ERRORS DONE
END-OF-OUTPUT

if [ $? != 0 ];
then
  cat ${output}
  exit 1
fi

exit 0
