#!/bin/bash
# Tests a simple inversion setup.
# Requires PISM's Python bindings and siple.
PYTHONEXEC=$5
PISM_BUILD_DIR=$1

# make sure that Python imports the right modules
export PYTHONPATH=${PISM_BUILD_DIR}/site-packages:$PYTHONPATH

# check if siple is installed
$PYTHONEXEC -c 'import siple'
if [ $? != 0 ];
then
    echo "Please make sure that siple is installed!"
    exit 1
fi

set -x
set -e

# Create input files
$PYTHONEXEC build_tiny.py -Mx 9 -My 9 -o tiny_nlcg.nc

$PYTHONEXEC make_synth_ssa.py -i tiny_nlcg.nc -o inv_data.nc \
              -pseudo_plastic -pseudo_plastic_q 0.25 -regional \
              -ssa_dirichlet_bc -generate_ssa_observed -ssa_method fem \
              -design_prior_const 70000 -inv_ssa tauc

# Run the inversion code
$PYTHONEXEC pismi.py \
              -i tiny_nlcg.nc -pseudo_plastic -pseudo_plastic_q 0.25 -inv_data inv_data.nc \
              -o tiny_inv.nc -regional -ssa_dirichlet_bc -inv_use_tauc_prior \
              -inv_design_param trunc -inv_design_cL2 1 -inv_design_cH1 0 \
              -inv_method nlcg -inv_target_misfit 100 -report_coverage

# Check if we succeeded
$PYTHONEXEC verify_ssa_inv.py tiny_inv.nc --desired_misfit 100 --iter_max 40 --morozov
