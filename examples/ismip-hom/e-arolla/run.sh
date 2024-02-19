#!/bin/bash

# Runs experiment E from the ISMIP-HOM benchmark using PISM's Blatter-Pattyn solver.

# The ISMIP-HOM paper uses 1e-16
A=$(echo "scale=50; 10^(-16) / (365 * 86400.0)" | bc -l)

input=$1
output=$2
# number of grid points in the X direction
Mx=$3
# number of MG levels
M=$4

# coarsening factor
C=4
# number of vertical levels
bp_Mz=$(echo "$C^($M - 1) + 1" | bc)

echo "Using ${bp_Mz} vertical levels..."

mpiexec -n 8 pismr -i ${input} -bootstrap \
      -Mx ${Mx} \
      -Mz 216 \
      -Lz 215 \
      -z_spacing equal \
      -grid.registration corner \
      -grid.periodicity y \
      -stress_balance.model blatter \
      -stress_balance.blatter.flow_law isothermal_glen \
      -flow_law.isothermal_Glen.ice_softness ${A} \
      -stress_balance.blatter.coarsening_factor ${C} \
      -blatter_Mz ${bp_Mz} \
      -bp_snes_monitor_ratio \
      -bp_ksp_type gmres \
      -bp_pc_type mg \
      -bp_pc_mg_levels ${M} \
      -bp_mg_levels_ksp_type richardson \
      -bp_mg_levels_pc_type sor \
      -bp_mg_coarse_ksp_type cg \
      -bp_mg_coarse_pc_type gamg \
      -basal_resistance.pseudo_plastic.enabled \
      -basal_resistance.pseudo_plastic.q 1.0 \
      -basal_resistance.pseudo_plastic.u_threshold 1m.s-1 \
      -basal_yield_stress.model constant \
      -energy none \
      -geometry.update.enabled false \
      -atmosphere uniform \
      -atmosphere.uniform.precipitation 0 \
      -surface simple \
      -y 1e-16 \
      -o_size big \
      -o ${output}
