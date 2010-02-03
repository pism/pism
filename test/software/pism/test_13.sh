#!/bin/bash

source ../functions.sh

# Test name:
test="Test #13: consistency of the temperature field at the ice-bedrock interface."
# The list of files to delete when done.
files="foo.nc temp.nc litho_temp.nc"
dir=`pwd`

run_test ()
{
    cleanup

    # do stuff
    
    run -n 2 pisms -y 1e3 -o foo.nc

    # do more stuff

    run ncks -v temp -d z,0m foo.nc temp.nc
    run ncks -v litho_temp -d zb,0m foo.nc litho_temp.nc
    run ncrename -O -v litho_temp,temp litho_temp.nc

    # compare results

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
