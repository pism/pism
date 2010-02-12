#!/bin/bash

source ../functions.sh

test="Test #4: regridding during bootstrapping: coarse -> fine -> coarse."
files="foo.nc bar.nc baz.nc"
dir=`pwd`

OPTS="-surface constant -Lz 4000 -y 0"
COARSE="-Mx 11 -My 11 -Mz 11"
FINE=" -Mx 21 -My 21 -Mz 21"
test_04 ()
{
    cleanup

    set -e

    # Create a file to bootstrap from:
    run -n 1 pismv -test G $COARSE -y 0 -o foo.nc

    # Coarse -> fine:
    run -n 1 pismr -boot_from foo.nc $FINE   $OPTS -o bar.nc

    # Fine -> coarse:
    run -n 1 pismr -boot_from bar.nc $COARSE $OPTS -o baz.nc

    set +e

    # Compare:
    run nccmp.py -v topg foo.nc baz.nc
    if [ $? != 0 ];
    then
	fail "foo.nc and baz.nc are different."
	return 1
    fi

    pass
    return 0
}

test_04
