#!/bin/bash
# Tests a simple inversion setup.
# Requires PISM's Python bindings
PYTHONEXEC=$5
PISM_BUILD_DIR=$1

# make sure that Python imports the right modules
export PYTHONPATH=${PISM_BUILD_DIR}/site-packages:$PYTHONPATH

set -x
set -e

# Create input files
$PYTHONEXEC build_tiny.py -Mx 9 -My 9

$PYTHONEXEC make_synth_ssa.py -i tiny.nc -o inv_data.nc \
              -pseudo_plastic -pseudo_plastic_q 0.25 -regional \
              -ssa_dirichlet_bc -generate_ssa_observed -ssa_method fem \
              -design_prior_const 70000 -inv_ssa tauc

# Run the inversion code
$PYTHONEXEC pismi.py \
              -i tiny.nc -pseudo_plastic -pseudo_plastic_q 0.25 -inv_data inv_data.nc \
              -o tiny_tikhonov_lmvm.nc -regional -ssa_dirichlet_bc -inv_use_tauc_prior \
              -inv_design_param trunc -inv_design_cL2 1 -inv_design_cH1 0 \
              -inv_method tikhonov_lmvm -tikhonov_penalty 6e-2 -report_coverage

# Check if we succeeded
$PYTHONEXEC verify_ssa_inv.py tiny_tikhonov_lmvm.nc --desired_misfit 10 --misfit_tolerance .5 --iter_max 100
