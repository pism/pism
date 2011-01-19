#!/bin/bash

source ../functions.sh

# Test name:
test="Test #19: verif test I regression: nonlinear SSA with plastic bed."
# The list of files to delete when done.
files="test-I-out.txt verify.nc"
dir=`pwd`

run_test ()
{
    cleanup

    # run test I
    OPTS="-test I -Mx 5 -Mz 31 -Mbz 1 -ssa_rtol 5e-07 -ksp_rtol 1e-12 -verbose 1 -o_size small"
    pismv -My 49  $OPTS  > test-I-out.txt
    pismv -My 193 $OPTS >> test-I-out.txt

    # compare results
    diff test-I-out.txt - > /dev/null <<END-OF-OUTPUT
NUMERICAL ERRORS in velocity relative to exact solution:
velocity  :  maxvector   prcntavvec      maxu       avu
               23.2923      0.98759   23.2923    7.6788
NUM ERRORS DONE
NUMERICAL ERRORS in velocity relative to exact solution:
velocity  :  maxvector   prcntavvec      maxu       avu
                1.3241      0.05675    1.3241    0.4413
NUM ERRORS DONE
END-OF-OUTPUT

    if [ $? != 0 ];
    then
	fail "numerical error report does not match the one stored"
	# the return statement *is* needed here, because 'fail' does not
	# terminate the test execution
	return 1
    fi

    pass
    return 0
}

run_test
