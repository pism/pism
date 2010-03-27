#!/bin/bash

source ../functions.sh

# Test name:
test="Test #12: mass conservation within rounding error for nonsliding SIA moving margin."
# The list of files to delete when done.
files="verify.nc ivol.nc dump.txt"
dir=`pwd`

run_test ()
{
    cleanup
    GRID="-Mx 31 -My 31 -Mz 31 -y 5000 -ys 1000"
    OPTS="-test B -max_dt 25 -o_size small"
    TS_OPTS="-ts_file ivol.nc -ts_vars ivol -ts_times 1000:25:1e4"
    # run test B 
    run -n 2 pismv -test B $GRID $OPTS $TS_OPTS

    ncdump -p 1,14 -v ivol -l 30 ivol.nc |sed -e "s/,//g; s/;//g" | uniq > dump.txt

    diff dump.txt - > /dev/null <<EOF
netcdf ivol {
dimensions:
	t = UNLIMITED  // (201 currently)
variables:
	double t(t) 
		t:units = "years" 
	double ivol(t) 
		ivol:units = "m3" 
		ivol:valid_min = 0.f 
		ivol:long_name = "total ice volume" 
data:

 ivol = 
    3.9740987142713e+15 
}
EOF

    if [ $? != 0 ];
    then
	fail "temperature fields in bar.nc and baz.nc are different"
	# the return statement *is* needed here, because 'fail' does not
	# terminate the test execution
	return 1
    fi

    pass
    return 0
}

run_test
