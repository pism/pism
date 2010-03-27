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
               6.585156  747.610002   89.594089    0.109410
base vels :  maxvector   avvector  prcntavvec     maxub     maxvb
                5.9669    0.29640     0.58383    5.9661    2.7798
NUM ERRORS DONE
NUMERICAL ERRORS evaluated at final time (relative to exact solution):
geometry  :    prcntVOL        maxH         avH   relmaxETA
               3.611327  720.424340   51.429248    0.058546
base vels :  maxvector   avvector  prcntavvec     maxub     maxvb
                1.8194    0.14192     0.27953    1.6178    1.2167
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
