#!/bin/bash

source ../functions.sh

test="Test #10: regridding with different number of processors."
dir=`pwd`
files="foo0.nc foo1.nc foo2.nc foo3.nc foo4.nc foo5.nc"

run_test ()
{
    cleanup

    set -e

    NRANGE="1 2 3 4 5"

    # Create a file to bootstrap from:
    run -n 1 pisms -eisII I -Mx 101 -My 201 -y 0 -o foo0.nc

    # Bootstrap:
    for NN in $NRANGE;
    do 
	run -n $NN pismr -boot_file foo0.nc -surface constant -Mx 101 -My 201 -Mz 11 -Lz 5000 -y 0 -o foo$NN.nc -o_size small
    done

    set +e

    # Compare:
    for i in $NRANGE;
    do
	for j in $NRANGE;
	do
	    if [ $i -le $j ]; then continue; fi
	    
	    run nccmp.py foo$i.nc foo$j.nc
	    if [ $? != 0 ];
	    then
		fail "Output files foo$i.nc and foo$j.nc are different."
		return 1
	    fi
	done
    done

    pass
    return 0
}

run_test
