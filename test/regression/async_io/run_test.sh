#!/bin/bash

pism_dir=${pism_dir:-$HOME/github/pism/pism}

run="python3"
run="coverage run"

mpirun -n 3 python3 ./async_io_test.py :\
       -n 1 ${run} ${pism_dir}/src/util/io/asynchronous_output_server.py

coverage html
