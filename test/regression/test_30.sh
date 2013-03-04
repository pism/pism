#!/bin/bash

# Test PISM's PDD model. Compare output from pclimate and pypdd.py.

PISM_PATH=$1
MPIEXEC=$2
PISM_SOURCE_DIR=$3

# List of files to remove when done:
files="i.nc atm.nc pypdd_smb.nc pism_smb.nc"

rm -f $files

set -e

# Run pypdd
$PISM_PATH/pypdd.py -b -o pypdd_smb.nc

# Prepare an input file for pclimate
$PISM_PATH/pisms -Lx 1 -Ly 1 -Mx 21 -My 21 -y 0 -o i.nc

# Run pclimate
$PISM_PATH/pclimate -i i.nc -o pism_smb.nc -times 0:1:1 -atmosphere given -atmosphere_given_file atm.nc -atmosphere_given_period 1 -surface pdd -pdd_annualize -pdd_std_dev 0.01

set +e

# Check results:
$PISM_PATH/nccmp.py -v climatic_mass_balance pypdd_smb.nc pism_smb.nc
if [ $? != 0 ];
then
    exit 1
fi

rm -f $files; exit 0
