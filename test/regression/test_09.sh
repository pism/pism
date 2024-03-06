#!/bin/bash

PISM_PATH=$1
MPIEXEC=$2

echo "Test # 9: 2D and 3D regridding from files with different variable orders."

# create a temporary directory and set up automatic cleanup
temp_dir=$(mktemp -d --tmpdir pism-test-09-XXXX)
trap 'rm -rf "$temp_dir"' EXIT
cd $temp_dir

set -e
set -x

# Create a file to bootstrap from (with a non-trivial bed topography):
# Note: Mx, My, and Mz should all be different.
$MPIEXEC -n 1 $PISM_PATH/pismr -eisII J -Mx 51 -My 60 -Mz 21 -Mbz 21 -Lbz 1000 -y 0 -o input-y,x,z.nc

# There are 6 possible variable orders.
#
# The y,x,z is the default that does not require transposing data. Note that this order is
# handled by the first iteration, creating the output file to compare all the other
# outputs to.
for order in y,x,z y,z,x x,y,z x,z,y z,x,y z,y,x;
do
  echo ${order}
  # Note: this is a no-op the first time through the loop.
  ncpdq -4 -L1 -O -a ${order} input-y,x,z.nc input-${order}.nc

  # We don't need to regrid all the possible variables: one 2D and one 3D variable is
  # enough.
  variables="topg,enthalpy"

  # Bootstrap from this file and run for 0 years:
  $MPIEXEC -n 2 $PISM_PATH/pismr \
           -Lz 4000 \
           -Mx 61 \
           -My 51 \
           -Mz 41 \
           -bootstrap \
           -i input-${order}.nc \
           -o o-${order}.nc \
           -regrid_file input-${order}.nc \
           -regrid_vars ${variables} \
           -y 0 \
           ""

  # Compare results
  $PISM_PATH/nccmp.py -v ${variables} o-y,x,z.nc o-${order}.nc
done
