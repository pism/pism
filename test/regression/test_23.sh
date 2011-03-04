#!/bin/bash

PISM_PATH=$1
MPIEXEC=$2
PISM_SOURCE_DIR=$3

# Test name:
test="Test #23: PST regression P1: SIA+SSA hybrid, flat bed, cold ice."
# The list of files to delete when done.
files="simp_exper.nc ser_pst_unnamedpst.nc P0A_small.nc P1_small50.nc"

rm -f $files

set +e

# prepare
cp $PISM_SOURCE_DIR/examples/pst/test/P0A_small.nc.gz .
cp $PISM_SOURCE_DIR/examples/pst/test/P1_small50.nc.gz .
gunzip P0A_small.nc.gz
gunzip P1_small50.nc.gz

# run pisms for PST exper:
$MPIEXEC -n 2 $PISM_PATH/pisms -pst -P1 -i P0A_small.nc -ys 0 -y 50 -skip 5

# compare; goal is 7 digit agreement; comment at right gives scale
$PISM_PATH/nccmp.py -v bwat -t 1e-7 simp_exper.nc P1_small50.nc # 1 m
if [ $? != 0 ]; then
    exit 1;
fi

$PISM_PATH/nccmp.py -v tauc -t 1e-0 simp_exper.nc P1_small50.nc # 10 bar = 1e6 Pa
if [ $? != 0 ]; then
    exit 2;
fi

$PISM_PATH/nccmp.py -v thk -t 1e-3 simp_exper.nc P1_small50.nc # 1000 m
if [ $? != 0 ]; then
    exit 3;
fi

$PISM_PATH/nccmp.py -v csurf -t 1e-3 simp_exper.nc P1_small50.nc # 100 m/a
if [ $? != 0 ]; then
    exit 4;
fi

$PISM_PATH/nccmp.py -v cbase -t 1e-3 simp_exper.nc P1_small50.nc # 100 m/a
if [ $? != 0 ]; then
    exit 5;
fi

rm -f $files; exit 0

