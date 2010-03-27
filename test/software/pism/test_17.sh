#!/bin/bash

source ../functions.sh

# Test name:
test="Test #17: Verification test L regression."
# The list of files to delete when done.
files="test-L-out.txt verify.nc"
dir=`pwd`

run_test ()
{
    cleanup

    # run test L
    OPTS="-test L -Mbz 1 -Mz 31 -y 1000 -o_size small -verbose 1"
    pismv -Mx 31 -My 31 $OPTS   > test-L-out.txt
    pismv -Mx 41 -My 41 $OPTS  >> test-L-out.txt

    # compare results
    diff test-L-out.txt - > /dev/null <<END-OF-OUTPUT
NUMERICAL ERRORS evaluated at final time (relative to exact solution):
geometry  :    prcntVOL        maxH         avH   relmaxETA
               0.030759  143.592042    3.831632    0.002932
NUM ERRORS DONE
NUMERICAL ERRORS evaluated at final time (relative to exact solution):
geometry  :    prcntVOL        maxH         avH   relmaxETA
               0.080920  154.561756    3.828041    0.002495
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
