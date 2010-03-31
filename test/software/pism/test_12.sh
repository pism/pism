#!/bin/bash

source ../functions.sh

# Test name:
test="Test #12: mass conservation within rounding error for nonsliding SIA moving margin."
# The list of files to delete when done.
files="verify.nc ivol.nc"
dir=`pwd`

run_test ()
{
    cleanup
    GRID="-Mx 31 -My 31 -Mz 31 -y 5000 -ys 1000"
    OPTS="-test B -max_dt 25 -o_size small"
    TS_OPTS="-ts_file ivol.nc -ts_vars ivol -ts_times 1000:25:1e4"
    # run test B 
    run -n 2 pismv -test B $GRID $OPTS $TS_OPTS

    python <<EOF
from numpy import diff, log10, floor
from sys import exit
try:
    from netCDF3 import Dataset
except:
    from netCDF4 import Dataset

nc = Dataset("ivol.nc", 'r')
ivol = nc.variables['ivol'][:]
ivol_max = ivol.max()
threshold = 10**(floor(log10(ivol_max)) - 14) # 14 digits of accuracy
diff_max = diff(ivol).max()

if diff_max < threshold:
    exit(0)
else:
    print "diff(ivol).max() = %f > %f" % (diff_max, threshold)
    exit(1)
EOF

    if [ $? != 0 ];
    then
	fail "ivol time series in ivol.nc is not constant"
	# the return statement *is* needed here, because 'fail' does not
	# terminate the test execution
	return 1
    fi

    pass
    return 0
}

run_test
