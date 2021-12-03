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
tiny=`mktemp -u tiny-XXXX.nc` || exit 1
$PYTHONEXEC build_tiny.py -Mx 9 -My 9 -o ${tiny}

inv_data=`mktemp -u inv-data-XXXX.nc` || exit 1
$PYTHONEXEC make_synth_ssa.py -i ${tiny} -o ${inv_data} \
              -pseudo_plastic -pseudo_plastic_q 0.25 -regional \
              -ssa_dirichlet_bc -generate_ssa_observed -ssa_method fem \
              -design_prior_const 70000 -inv_ssa tauc

# Run the inversion code
output=`mktemp -u tiny-tikhonov-lmvm-XXXX.nc` || exit 1
$PYTHONEXEC pismi.py \
              -i ${tiny} -pseudo_plastic -pseudo_plastic_q 0.25 -inv_data ${inv_data} \
              -o ${output} -regional -ssa_dirichlet_bc -inv_use_tauc_prior \
              -inv_design_param trunc -inv_design_cL2 1 -inv_design_cH1 0 \
              -inv_method tikhonov_lmvm -tikhonov_penalty 6e-2 -report_coverage

# Check if we succeeded
$PYTHONEXEC verify_ssa_inv.py ${output} --desired_misfit 10 --misfit_tolerance .5 --iter_max 100

# Clean up
rm -f ${tiny} ${inv_data} ${output}
