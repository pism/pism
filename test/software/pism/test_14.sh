#!/bin/bash

source ../functions.sh

# Test name:
test="Test #14: Verification test E regression."
# The list of files to delete when done.
files="test-E-out.txt verify.nc"
dir=`pwd`

run_test ()
{
    cleanup

    # run test E
    OPTS="-test E -y 1000 -o_size small -verbose 1 -Mbz 1"
    pismv -Mx 21 -My 21 -Mz 21 $OPTS  > test-E-out.txt
    pismv -Mx 41 -My 41 -Mz 41 $OPTS >> test-E-out.txt

    # compare results
    diff test-E-out.txt - > /dev/null <<END-OF-OUTPUT
NUMERICAL ERRORS evaluated at final time (relative to exact solution):
geometry  :    prcntVOL        maxH         avH   relmaxETA
               6.573233  747.157559   89.430387    0.109307
base vels :  maxvector   avvector  prcntavvec     maxub     maxvb
                5.9575    0.29524     0.58154    5.9567    2.7744
NUM ERRORS DONE
NUMERICAL ERRORS evaluated at final time (relative to exact solution):
geometry  :    prcntVOL        maxH         avH   relmaxETA
               3.610547  720.403968   51.418143    0.058543
base vels :  maxvector   avvector  prcntavvec     maxub     maxvb
                1.8189    0.14184     0.27939    1.6172    1.2165
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
