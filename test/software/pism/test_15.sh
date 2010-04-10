#!/bin/bash

source ../functions.sh

# Test name:
test="Test #15: Verification test B regression."
# The list of files to delete when done.
files="test-B-out.txt verify.nc"
dir=`pwd`

run_test ()
{
    cleanup

    # run test B
    OPTS="-test B -Mbz 1 -ys 422.45 -y 2500 -o_size small -verbose 1"
    pismv -Mx 31 -My 31 -Mz 31 $OPTS  > test-B-out.txt
    pismv -Mx 41 -My 41 -Mz 41 $OPTS >> test-B-out.txt

    # compare results
    diff test-B-out.txt - > /dev/null <<END-OF-OUTPUT
NUMERICAL ERRORS evaluated at final time (relative to exact solution):
geometry  :    prcntVOL        maxH         avH   relmaxETA
               0.075974  276.557369   11.445854    0.030467
NUM ERRORS DONE
NUMERICAL ERRORS evaluated at final time (relative to exact solution):
geometry  :    prcntVOL        maxH         avH   relmaxETA
               0.384904  211.821582    5.809095    0.011945
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
