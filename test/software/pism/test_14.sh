#!/bin/bash

source ../functions.sh

# Test name:
test="Test #14: Verification test A regression."
# The list of files to delete when done.
files="test-A-out.txt verify.nc"
dir=`pwd`

run_test ()
{
    cleanup

    # run test A
    OPTS="-test A -y 1000 -o_size small -verbose 1 -Mbz 1"
    pismv -Mx 21 -My 21 -Mz 21 $OPTS  > test-A-out.txt
    pismv -Mx 41 -My 41 -Mz 41 $OPTS >> test-A-out.txt

    # compare results
    diff test-A-out.txt - > /dev/null <<END-OF-OUTPUT
NUMERICAL ERRORS evaluated at final time (relative to exact solution):
geometry  :    prcntVOL        maxH         avH   relmaxETA
               6.761003  746.932951   91.953306    0.109256
NUM ERRORS DONE
NUMERICAL ERRORS evaluated at final time (relative to exact solution):
geometry  :    prcntVOL        maxH         avH   relmaxETA
               3.703364  720.526587   52.739966    0.058561
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
