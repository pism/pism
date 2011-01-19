#!/bin/bash

source ../functions.sh

# Test name:
test="Test #20: verif test J regression: linearized SSA, floating."
# The list of files to delete when done.
files="test-J-out.txt verify.nc"
dir=`pwd`

run_test ()
{
    cleanup

    # run test J
    OPTS="-test J -Mz 11 -Mbz 1 -pc_type asm -sub_pc_type lu -ksp_rtol 1e-12 -verbose 1 -o_size small"
    pismv -Mx 49 -My 49 $OPTS  > test-J-out.txt
    pismv -Mx 98 -My 98 $OPTS >> test-J-out.txt

    # compare results
    diff test-J-out.txt - > /dev/null <<END-OF-OUTPUT
NUMERICAL ERRORS in velocity relative to exact solution:
velocity  :  maxvector   prcntavvec      maxu      maxv       avu       avv
                0.2412      0.08327    0.2311    0.0696    0.1430    0.0384
NUM ERRORS DONE
NUMERICAL ERRORS in velocity relative to exact solution:
velocity  :  maxvector   prcntavvec      maxu      maxv       avu       avv
                0.0604      0.02081    0.0578    0.0174    0.0357    0.0096
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
