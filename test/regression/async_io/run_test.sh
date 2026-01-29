#!/bin/bash

pism_dir=${pism_dir:-$HOME/github/pism/pism}

run="python3"
# run="coverage run"

# Note: async_io_test.py should use "-n X" for some X > 1. This is needed to make sure
# that asynchronous_output_server.py uses the code assembling parts of the grid received
# from different ranks (on the PISM side).
mpirun -n 3 python3 ./async_io_test.py :\
       -n 1 ${run} ${pism_dir}/src/util/io/asynchronous_output_server.py -d

# coverage html
