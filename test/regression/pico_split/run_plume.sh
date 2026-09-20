#!/bin/bash

# Copyright (C) 2026 PISM authors

set -e
set -u
set -x

###############################################################################
# Smoke test of the plume ocean model (PICOP's plume without the PICO box model):
# run it with theta_ocean as potential temperature and as thermal forcing and check
# that both produce sub-shelf melt.
###############################################################################

PISM_PATH=$1
MPIEXEC=$2
PISM_SOURCE_DIR=$3

input_file=${PISM_SOURCE_DIR}/test/regression/pico_split/bedmap2_schmidtko14_50km.nc

temp_dir=$(mktemp -d -t pism_plume.XXXX) || exit 1

pushd ${temp_dir}

# A thermal forcing version of the input: theta_ocean is in degrees Celsius and lies
# roughly 1.9 degrees above the freezing point.
ncap2 -O -s "theta_ocean=theta_ocean+1.9f" ${input_file} input_tf.nc

grid="-bootstrap -Mx 120 -My 120 -Lz 6000 -Lbz 2000 -Mz 81 -Mbz 21 -grid.recompute_longitude_and_latitude false"

stressbalance="-pik -stress_balance ssa+sia -ssa_method fd"
surface="-atmosphere uniform -surface simple"

max_melt() {
# Prints the maximum of plume_basal_melt_rate in a spatial output file.
ncwa -O -y max -v plume_basal_melt_rate $1 max.nc
ncks -H -s "%g" -v plume_basal_melt_rate max.nc
}

check() {
# Fails the test unless the melt rate in a spatial output file is positive somewhere.
m=$(max_melt $1)
if awk -v m="$m" 'BEGIN { exit !(m > 0) }'
then
  echo "$1: max melt rate $m m/year, OK"
else
  echo "$1: max melt rate $m m/year, FAILED"
  exit 1
fi
}

${PISM_PATH}/pism -verbose 2 -i ${input_file} \
            -config ${PISM_PATH}/pism_config.nc \
            $grid $stressbalance $surface \
            -ocean plume -ocean.plume.file ${input_file} \
            -y 0.001 -spatial_file ex_theta.nc -spatial_times 0.001 \
            -spatial_vars plume_basal_melt_rate,plume_temperature \
            -o o_theta.nc

${PISM_PATH}/pism -verbose 2 -i ${input_file} \
            -config ${PISM_PATH}/pism_config.nc \
            $grid $stressbalance $surface \
            -ocean plume -ocean.plume.file input_tf.nc \
            -ocean.plume.temperature_as_thermal_forcing \
            -y 0.001 -spatial_file ex_tf.nc -spatial_times 0.001 \
            -spatial_vars plume_basal_melt_rate,plume_temperature \
            -o o_tf.nc

check ex_theta.nc
check ex_tf.nc

rm -rf ${temp_dir}
