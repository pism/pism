#!/bin/bash

# set computer and user specific variables,
# this is sourced from the pism run script, so that
# the pism run script can be shared without effort.

export pismcode_dir=$PISM_DIR
export pism_exec=$PISM_DIR/bin/pismr
export working_dir=/Users/flo/Recherche/MODELS/WORK_PISM/ANT_paleo_spinup_30km
export input_data_dir=/Users/flo/Recherche/MODELS/PISM_DATA/pism-ais
export pism_mpi_do="mpiexec -n"

