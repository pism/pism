#!/bin/bash

# set computer and user specific variables,
# this is sourced from the pism run script, so that
# the pism run script can be shared without effort.

export pismcode_dir=/home/pism/version
export pism_exec=$pismcode_dir/bin/pismr
export working_dir=/home/pism/results
export input_data_dir=/home/pism/inputdata
export pism_mpi_do="srun -n"

