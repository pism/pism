#!/bin/bash

# Testing temporal averaging of the load used by bed deformation models.

PISM_PATH=$1
MPIEXEC=$2
PISM_SOURCE_DIR=$3

if [ $# -ge 4 ] && [ "$4" == "-python" ]
then
  PYTHONEXEC=$5
  export PYTHONPATH=${PISM_PATH}/site-packages:${PYTHONPATH}
else
  exit 1
fi

pism="$PISM_PATH/pism"

# create a temporary directory and set up automatic cleanup
temp_dir=$(mktemp -d --tmpdir pism-test-XXXX)
trap 'rm -rf "$temp_dir"' EXIT
cd $temp_dir

# Make sure PISM can find the configuration file
echo "
-config ${PISM_PATH}/pism_config.nc
" > .petscrc

set -e
set -u
set -x

${PYTHONEXEC} <<EOF
import PISM
import PISM.testing as pt
import numpy as np

ctx = PISM.Context().ctx

ice_density = ctx.config().get_number("constants.ice.density")

Lx = 500e3
M = 5
grid = PISM.Grid.Shallow(ctx, Lx, Lx, 0, 0, M, M, PISM.CELL_CENTER, PISM.NOT_PERIODIC)

def forcing(values, filename):
    days_per_year = 360
    times = np.array([2, 6, 10, 14]) / 16 * days_per_year
    time_bounds = np.array([0, 4, 8, 12, 16]) / 16 * days_per_year

    pt.create_forcing(grid,
                      filename,
                      "delta_P",
                      "kg / (m^2 360 day)",
                      values,
                      time_units="days since 1-1-1",
                      times=times,
                      time_bounds=time_bounds)


def input_data(thickness, filename):
    bed = PISM.Scalar(grid, "bed")
    # well above sea level to avoid "flooding" a part of it and losing SMB (it is not applied
    # to ice free ocean)
    bed.set(1000.0)
    bed.metadata(0).standard_name("bedrock_altitude").units("m")

    thk = PISM.Scalar(grid, "thk")
    thk.set(thickness)
    thk.metadata(0).standard_name("land_ice_thickness").units("m")

    try:
        f = PISM.util.prepare_output(filename, calendar="360_day")
        thk.write(f)
        bed.write(f)
    finally:
        f.close()

forcing(np.array([10.0, -10, -10.0, 10.0]) * ice_density, "dP_PMMP.nc")
forcing(np.array([-10, -10, 10, 10]) * ice_density, "dP_MMPP.nc")
forcing(np.array([0, 0, 0, 0]) * ice_density, "dP_0.nc")

input_data(10, "H-10.nc")
input_data(12.5, "H-12.5.nc")
EOF


common_options="
-atmosphere uniform,delta_P
-atmosphere.delta_P.periodic
-atmosphere.uniform.precipitation 0
-bed_def lc
-bed_deformation.update_interval 360days
-bootstrap
-calendar 360_day
-energy none
-max_dt 90days
-stress_balance none
-surface simple
-ys 1-1-1
"

extra="
-extra_times 720days
-extra_vars topg
"

# 1. Bootstrap from the file with H=10m and run for N years with SMB=[10, -10, -10, 10]
#
# 2.1. Bootstrap from the file with H=10m. Run for 0s to get viscous and elastic bed
#      displacements corresponding to equilibrium with the ice thickness H=10m. Save to
#      H-10-full.nc
#
#      This step is needed so that in 2.2 below the bed deformation model is in
#      equilibrium with H=10m *at the beginning of the run*. Otherwise it would be in
#      equilibrium with H=12.5m and the load of 10m would result in some uplift.
#
# 2.2. Bootstrap from H-12.5.nc and regrid elastic and viscous displacement from H-10-full.nc
#      and run for N years with SMB=[-10, 10, 10, -10]
#
# 3. Bootstrap from the file with H=10m. Run for N years with SMB=0.
#
# In the first two cases the SMB over 1 year is zero and ice thickness oscillates between
# 7.5 m and 12.5 m with the average of 10m. The average load is the same, but the load at
# the end of the year is different.

# Run for 2 years. Here we update bed elevation every year, so this is enough.
N=2
run_length=$(( N * 360 ))

# 1. Bootstrap from the file with H=10m and run for N years with SMB=[10, -10, -10, 10]
#
${pism} \
   ${common_options} \
   -atmosphere.delta_P.file dP_PMMP.nc \
   ${extra} -extra_file ex_PMMP.nc \
   -i H-10.nc \
   -time.run_length ${run_length}days \
   -o_size none \
   ""

# 2.1. Bootstrap from the file with H=10m. Run for 0s to get viscous and elastic bed
#      displacements corresponding to equilibrium with the ice thickness H=10m. Save to
#      H-10-full.nc
#
${pism} \
   ${common_options} \
   -atmosphere.delta_P.file dP_PMMP.nc \
   -i H-10.nc \
   -o H-10-full.nc \
   -time.run_length 0 \
   ""

# 2.2. Bootstrap from H-12.5.nc, then regrid elastic and viscous displacement from
#      H-10-full.nc and run for N years with SMB=[-10, -10, 10, 10]
#
${pism} \
   ${common_options} \
   -atmosphere.delta_P.file dP_MMPP.nc \
   ${extra} -extra_file ex_MMPP.nc \
   -i H-12.5.nc \
   -regrid_file H-10-full.nc \
   -regrid_vars viscous_bed_displacement,elastic_bed_displacement \
   -time.run_length ${run_length}days \
   -o_size none \
   ""

# 3. Bootstrap from the file with H=10m. Run for N years with SMB=0.
#
${pism} \
   ${common_options} \
   -atmosphere.delta_P.file dP_0.nc \
   ${extra} -extra_file ex_0.nc \
   -i H-10.nc \
   -time.run_length ${run_length}days \
   -o_size none \
   ""

# Compare results
$PISM_PATH/pism_nccmp -v topg ex_0.nc ex_PMMP.nc && $PISM_PATH/pism_nccmp -v topg ex_0.nc ex_MMPP.nc
