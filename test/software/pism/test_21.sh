#!/bin/bash

source ../functions.sh

# Test name:
test="Test #21: Verification test K regression."
# The list of files to delete when done.
files="test-K-out.txt verify.nc"
dir=`pwd`

run_test ()
{
    cleanup

    # run test K
    OPTS="-test K -Mx 4 -My 4 -y 13000.0 -Lbz 1000 -verbose 1 -o_size small"
    pismv -Mz 41 -Mbz 11 $OPTS  > test-K-out.txt
    pismv -Mz 81 -Mbz 21 $OPTS >> test-K-out.txt

    # compare results
    diff test-K-out.txt - > /dev/null <<END-OF-OUTPUT
NUMERICAL ERRORS evaluated at final time (relative to exact solution):
temp      :        maxT         avT       maxTb        avTb
               0.376380    0.069534    0.037720    0.032864
NUM ERRORS DONE
NUMERICAL ERRORS evaluated at final time (relative to exact solution):
temp      :        maxT         avT       maxTb        avTb
               0.052752    0.008728    0.007214    0.006248
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
