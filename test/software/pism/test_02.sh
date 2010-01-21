#!/bin/bash

source ../functions.sh

test="Test #2: pismv exact processor independence (SIA only; test G)."
dir=`pwd`
files="foo1.nc foo2.nc foo3.nc foo4.nc foo6.nc"

test_02 ()
{
    cleanup

    set -e

    # Create the files:
    for NN in 1 2 3 4 6;
    do 
	run -n $NN pismv -test G -Mx 61 -My 61 -Mz 61 -y 1 -verbose 1 -o foo$NN.nc
    done

    set +e

    # Compare:
    for i in 1 2 3 4 6;
    do
	for j in 1 2 3 4 6;
	do
	    if [ $i -le $j ]; then continue; fi
	    
	    run nccmp.py -x -v rank foo$i.nc foo$j.nc
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

test_02
