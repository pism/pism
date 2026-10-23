#!/bin/bash

set -e
set -u

build_dir=${1:?ERROR: PISM build directory not specified.}
input_file=${2:?ERROR: Input file not specified.}

# Run clang-tidy
clang-tidy -p ${build_dir} ${input_file}
