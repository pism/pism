#!/bin/bash

set -e
set -u

build_dir=${1:?ERROR: PISM build directory not specified.}
input_file=${2:?ERROR: Input file not specified.}

# Extract compiler options specifying MPI's include directories and
# pass each one as a separate --extra-arg to clang-tidy
extra_args=$(mpicc -show | grep -E -o -e "-I[^ ]+" | sed 's/^/--extra-arg=/')

# Run clang-tidy
clang-tidy ${extra_args} -p ${build_dir} ${input_file}
