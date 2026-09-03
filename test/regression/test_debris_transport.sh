#!/bin/bash

# Test the debris transport model on the idealized valley glacier of examples/debris:
# - a run can be stopped and re-started without changing the result,
# - the debris mass budget closes.

set -u

PISM_PATH=$1
MPIEXEC=$2
PISM_SOURCE_DIR=$3
PYTHON=${5:-python3}

# create a temporary directory and set up automatic cleanup
temp_dir=$(mktemp -d --tmpdir pism-test-XXXX)
trap 'rm -rf "$temp_dir"' EXIT
cd $temp_dir

set -e

# a coarse version of the example
${PYTHON} ${PISM_SOURCE_DIR}/examples/debris/create_input.py --dx 100 --debris-file debris_input.nc input.nc

bootstrap="
        -input.bootstrap
        -input.file input.nc
        -grid.Mz 21
        -grid.Lz 300
        -time.start 0
"

common_options="
        -stress_balance.model sia
        -stress_balance.sia.flow_law isothermal_glen
        -stress_balance.sia.bed_smoother.range 0
        -flow_law.isothermal_Glen.ice_softness 1.0e-24
        -energy.model none
        -surface.models given
        -geometry.front_retreat.prescribed.file input.nc
        -debris.models transport
        -debris.transport.input.file debris_input.nc
        -debris.transport.marginal_length_scale 100
        -debris.ice_melt_enhancement.model verhaegen
        -time_stepping.maximum_time_step 1
        -output.sizes.medium debris_thickness,englacial_debris_concentration,ice_melt_enhancement
"

scalars=total_debris_mass,englacial_debris_mass,supraglacial_debris_mass,debris_input_mass_flux,debris_mass_conservation_error

# number of MPI processes:
N=4

end_time=20
stop_time=10

set -x

# from 0 to A:
${MPIEXEC} -n ${N} ${PISM_PATH}/pism \
           ${bootstrap} ${common_options} \
           -output.file o_part1.nc \
           -time.end ${stop_time}
# and from A to B:
${MPIEXEC} -n ${N} ${PISM_PATH}/pism \
           ${common_options} \
           -input.file o_part1.nc \
           -time.end ${end_time} \
           -output.file o_interrupted.nc

# Now run from 0 to B and compare:
${MPIEXEC} -n ${N} ${PISM_PATH}/pism \
           ${bootstrap} ${common_options} \
           -output.file o_uninterrupted.nc \
           -output.scalar.file scalar.nc \
           -output.scalar.times 1 \
           -output.scalar.variables ${scalars} \
           -time.end ${end_time}

set +e

ignored_vars=wall_clock_time,step_counter,pism_config,model_years_per_processor_hour

# Compare results:
$PISM_PATH/pism_nccmp -x -v ${ignored_vars} o_uninterrupted.nc o_interrupted.nc
if [ $? != 0 ]; then
  exit 1
fi

# Check the mass budget: the source is switched on at year 5, so debris has to be present
# and the budget has to close.
${PYTHON} - <<'EOF'
import netCDF4, numpy as np, sys
f = netCDF4.Dataset("scalar.nc")
total = f.variables["total_debris_mass"][:]
error = f.variables["debris_mass_conservation_error"][:]
flux  = f.variables["debris_input_mass_flux"]
# the flux is written in "glaciological" units (per year)
per_year = flux[:] * (1.0 if "year" in flux.units else 365.0 * 86400.0)
print("total debris mass (kg):", np.asarray(total))
print("conservation error (kg):", np.asarray(error))
print("input mass flux (kg/year):", np.asarray(per_year))
if not total[-1] > 1e6:
    print("FAILED: too little debris at the end of the run")
    sys.exit(1)
if not abs(per_year[-1] - 2e6) < 1e-6 * 2e6:
    print("FAILED: the input mass flux should be 2e6 kg/year, got", per_year[-1])
    sys.exit(1)
if not np.all(np.abs(error) < 1e-6 * max(total[-1], 1.0)):
    print("FAILED: the debris mass budget does not close")
    sys.exit(1)
print("debris mass budget OK")
EOF
