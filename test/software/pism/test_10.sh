#!/bin/bash

source ../functions.sh

test="Test #10: regridding with different number of processors."
dir=`pwd`
files="foo0.nc foo1.nc foo2.nc foo3.nc foo8.nc foo10.nc"

run_test ()
{
    cleanup

    set -e

    # Create a file to bootstrap from:
    run -n 1 pisms -eisII I -Mx 101 -My 201 -y 0 -o foo0.nc

    # Bootstrap:
    for NN in 1 2 3 4 5;
    do 
	run -n $NN pismr -boot_from foo0.nc -surface constant -Mx 101 -My 201 -Mz 11 -Lz 5000 -y 0 -o foo$NN.nc -o_size small
    done

    set +e

    # Compare:
    for i in 1 2 3 4 5;
    do
	for j in 1 2 3 4 5;
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