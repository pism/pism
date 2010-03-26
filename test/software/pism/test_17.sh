#!/bin/bash

source ../functions.sh

# Test name:
test="Test #17: Verification test F regression."
# The list of files to delete when done.
files="test-F-out.txt verify.nc"
dir=`pwd`

run_test ()
{
    cleanup

    # run test F
    OPTS="-test F -Mbz 1 -Mz 31 -y 1000 -o_size small -verbose 1"
    pismv -Mx 31 -My 31 $OPTS   > test-F-out.txt
    pismv -Mx 41 -My 41 $OPTS  >> test-F-out.txt

    # compare results
    diff test-F-out.txt - > /dev/null <<END-OF-OUTPUT
NUMERICAL ERRORS evaluated at final time (relative to exact solution):
geometry  :    prcntVOL        maxH         avH   relmaxETA
               0.000792    7.969510    0.170437    0.007100
temp      :        maxT         avT    basemaxT     baseavT
               0.774360    0.054981    0.204768    0.014161
Sigma     :      maxSig       avSig
               2.756683    0.096017
surf vels :     maxUvec      avUvec        maxW         avW
               0.150476    0.016189    0.004862    0.000260
NUM ERRORS DONE
NUMERICAL ERRORS evaluated at final time (relative to exact solution):
geometry  :    prcntVOL        maxH         avH   relmaxETA
               0.001215    5.628235    0.112949    0.005011
temp      :        maxT         avT    basemaxT     baseavT
               0.789185    0.054923    0.414673    0.014309
Sigma     :      maxSig       avSig
               2.458317    0.056954
surf vels :     maxUvec      avUvec        maxW         avW
               0.117184    0.009892    0.004205    0.000170
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
