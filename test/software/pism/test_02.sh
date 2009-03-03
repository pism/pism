#!/bin/bash

source ../functions.sh

test="Test #2: pismv exact processor independence for SIA runs."
dir=`pwd`
files="foo1.nc foo2.nc foo3.nc foo8.nc foo10.nc"

test_02 ()
{
    cleanup

    # Create the files:
    for NN in 1 2 3 8 10;
    do 
	run mpiexec -n $NN pismv -test G -Mx 61 -My 61 -Mz 61 -y 1 -verbose 1 -o foo$NN.nc
    done

    # Compare:
    for i in 1 2 3 8 10;
    do
	for j in 1 2 3 8 10;
	do
	    if [ $i == $j ]; then continue; fi
	    
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

test_02