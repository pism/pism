#!/bin/bash

source ../functions.sh

# Test name:
test="Test #18: verification test G regression: thermo SIA w. time-dependent mass balance."
# The list of files to delete when done.
files="test_18-G-out.txt verify.nc"
dir=`pwd`

run_test ()
{
    cleanup

    # run test G
    OPTS="-test G -Mbz 1 -Mz 31 -y 1000 -o_size small -verbose 1"
    pismv -Mx 31 -My 31 $OPTS   > test_18-G-out.txt
    pismv -Mx 41 -My 41 $OPTS  >> test_18-G-out.txt

    # compare results
    diff test_18-G-out.txt - > /dev/null <<END-OF-OUTPUT
NUMERICAL ERRORS evaluated at final time (relative to exact solution):
geometry  :    prcntVOL        maxH         avH   relmaxETA
               0.780389   32.430770    7.148950    0.016288
temp      :        maxT         avT    basemaxT     baseavT
               0.835718    0.250497    0.759924    0.153262
Sigma     :      maxSig       avSig
               8.570017    0.996601
surf vels :     maxUvec      avUvec        maxW         avW
               0.944836    0.200079    0.028356    0.004029
NUM ERRORS DONE
NUMERICAL ERRORS evaluated at final time (relative to exact solution):
geometry  :    prcntVOL        maxH         avH   relmaxETA
               0.782888   30.763934    7.242048    0.015469
temp      :        maxT         avT    basemaxT     baseavT
               0.895269    0.249036    0.754972    0.157030
Sigma     :      maxSig       avSig
               8.211503    0.945999
surf vels :     maxUvec      avUvec        maxW         avW
               0.885743    0.194988    0.027397    0.004220
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
