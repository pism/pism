#!/bin/bash

# Test YAC-based interpolation from a grid too small to distribute across 2 MPI processes.

PISM_PATH=$1
MPIEXEC=$2
PISM_SOURCE_DIR=$3

# create a temporary directory and set up automatic cleanup
temp_dir=$(mktemp -d --tmpdir pism-interp-small-grid-XXXX)
trap 'rm -rf "$temp_dir"' EXIT
cd $temp_dir

set -e

# Note: a 3x3 grid cannot be distributed across 2 or more MPI processes (in PISM) because
# each sub-domain has to be at least stencil_width wide and PISM uses stencil widths of 2
# (most runs) and 3 (runs using eigen_calving).
#
# This restriction comes from PETSc: if the width of a sub-domain is smaller than the
# stencil width communicating ghosts would require exchanging messages with other ranks
# *in addition to* direct neighbors.

cat > bed.cdl <<EOF
netcdf bed {
dimensions:
 y = 3 ;
 x = 3 ;
variables:
 short topg(y, x) ;
  topg:grid_mapping = "mapping" ;
  topg:standard_name = "bedrock_altitude" ;
  topg:units = "m" ;
 double mapping ;
  mapping:ellipsoid = "WGS84" ;
  mapping:grid_mapping_name = "polar_stereographic" ;
  mapping:false_easting = 0. ;
  mapping:false_northing = 0. ;
  mapping:latitude_of_projection_origin = -90. ;
  mapping:standard_parallel = -71. ;
  mapping:straight_vertical_longitude_from_pole = 0. ;
 int x(x) ;
  x:axis = "X" ;
  x:standard_name = "projection_x_coordinate" ;
  x:units = "m" ;
 int y(y) ;
  y:axis = "Y" ;
  y:standard_name = "projection_y_coordinate" ;
  y:units = "m" ;

data:
   x = -3022600, 0, 3022600;
   y = -3022600, 0, 3022600;
   topg = 0, 1, 2,
          1, 2, 3,
          2, 3, 4;
}

EOF

# generate a regridding file:
ncgen -o bed.nc bed.cdl

input_file=${PISM_SOURCE_DIR}/test/regression/pico_split/bedmap2_schmidtko14_50km.nc

options="
   -Lz 5000
   -atmosphere uniform
   -bootstrap
   -energy none
   -i ${input_file}
   -regrid_file bed.nc
   -regrid_vars topg
   -stress_balance none
   -surface simple
   -y 1s
"

# Note: we don't check interpolation results here, just ensure that PISM finishes
# without an error message.
${MPIEXEC} -n 2 ${PISM_PATH}/pism ${options} -o foo.nc
