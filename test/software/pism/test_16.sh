#!/bin/bash

source ../functions.sh

# Test name:
test="Test #16: Verification test C regression."
# The list of files to delete when done.
files="test_16-C-out.txt verify.nc"
dir=`pwd`

run_test ()
{
    cleanup

    # run test C
    OPTS="-test C -Mbz 1 -Mz 31 -y 5000 -o_size small -verbose 1"
    pismv -Mx 31 -My 31 $OPTS  > test_16-C-out.txt
    pismv -Mx 41 -My 41 $OPTS >> test_16-C-out.txt

    # compare results
    diff test_16-C-out.txt - > /dev/null <<END-OF-OUTPUT
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
	fail "numerical error report does not match the one stored"
	# the return statement *is* needed here, because 'fail' does not
	# terminate the test execution
	return 1
    fi

    pass
    return 0
}

run_test
