#!/bin/bash

# Runs experiment E from the ISMIP-HOM benchmark using PISM's Blatter-Pattyn solver.

# The ISMIP-HOM paper uses 1e-16
A=$(echo "scale=50; 10^(-16) / (365 * 86400.0)" | bc -l)

input=$1
output=$2

# This run is *serial* to make sure it works with most PETSc installations. The equivalent
# parallel run would require PETSc built with a parallel direct solver (e.g. MUMPS).

pismr -i ${input} -bootstrap \
      -Mx 401 \
      -grid.registration corner \
      -grid.periodicity y \
      -stress_balance.model blatter \
      -stress_balance.blatter.flow_law isothermal_glen \
      -flow_law.isothermal_Glen.ice_softness ${A} \
      -stress_balance.blatter.coarsening_factor 7 \
      -blatter_Mz 50 \
      -bp_snes_monitor_ratio \
      -bp_ksp_type gmres \
      -bp_pc_type mg \
      -bp_pc_mg_levels 3 \
      -bp_mg_levels_ksp_type richardson \
      -bp_mg_levels_pc_type sor \
      -bp_mg_coarse_ksp_type preonly \
      -bp_mg_coarse_pc_type lu \
      -basal_resistance.pseudo_plastic.enabled \
      -basal_resistance.pseudo_plastic.q 1.0 \
      -basal_resistance.pseudo_plastic.u_threshold 3.1556926e7 \
      -basal_yield_stress.model constant \
      -energy none \
      -geometry.update.enabled false \
      -atmosphere uniform \
      -atmosphere.uniform.precipitation 0 \
      -surface simple \
      -y 1e-16 \
      -o ${output}
