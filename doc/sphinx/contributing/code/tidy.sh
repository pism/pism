#!/bin/bash

set -e
set -u

build_dir=${1:?ERROR: PISM build directory not specified.}
input_file=${2:?ERROR: Input file not specified.}

# Extract compiler options specifying MPI's include directories
mpi=$(mpicc -show | grep -E -o -e "-I[^ ]+")

# Run clang-tidy
clang-tidy --extra-arg=${mpi} -p ${build_dir} ${input_file}
