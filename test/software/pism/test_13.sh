#!/bin/bash

source ../functions.sh

# Test name:
test="Test #13: temperature continuity at ice-bed interface (cold case)."
# The list of files to delete when done.
files="foo.nc temp.nc litho_temp.nc"
dir=`pwd`

run_test ()
{
    cleanup

    # generate sample state with cold base, but which has seen some strain-heating
    run -n 2 pisms -eisII A -Mx 31 -My 31 -Mz 31 -Lbz 1000 -Mbz 11 \
               -y 6e3 -temp_pa -o foo.nc

    # extract only what is needed for comparison
    run ncks -v temp -d z,0m foo.nc temp.nc
    run ncks -v litho_temp -d zb,-0.000001m foo.nc litho_temp.nc
    # neither of these seems to work:
    #run ncks -v litho_temp -d zb,0m foo.nc litho_temp.nc
    #run ncks -v litho_temp -d zb,-0m foo.nc litho_temp.nc

    # compare
    run ncrename -O -v litho_temp,temp litho_temp.nc
    run nccmp.py -v temp temp.nc litho_temp.nc

    if [ $? != 0 ];
    then
	fail "Basal ice temperature and bedrock temperature are different"
	# the return statement *is* needed here, because 'fail' does not
	# terminate the test execution
	return 1
    fi

    pass
    return 0
}

run_test
