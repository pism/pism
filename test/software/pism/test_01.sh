#!/bin/bash

source ../functions.sh

test="Test #1: pismv exact restartability."
files="verify.nc foo.nc bar.nc foo.txt bar.txt"
dir=`pwd`

test_01 ()
{
    cleanup

    pismv -test G -Mx 61 -My 61 -Mz 61 -y 10 -max_dt 1 -verbose 1 -o foo.nc > foo.txt
    run pismv -test G -Mx 61 -My 61 -Mz 61 -y 5 -max_dt 1 -verbose 1
    # Fix the following line and uncomment the block below:
    pismv -test G -Mx 61 -My 61 -Mz 61 -if verify.nc -y 5 -max_dt 1 -verbose 1 -o bar.nc > bar.txt

    # Compare output files:
    nccmp.py -t 1e-6 foo.nc bar.nc > /dev/null
    if [ ! $? ];
    then
	fail "Output files are different."
	return 1
    fi

    # Compare numerical error reports:
#     if [ -n "$(diff foo.txt bar.txt)" ];
#     then
# 	fail "Numerical error reports are different."
# 	return 1
#     fi

    pass
    return 0
}

test_01