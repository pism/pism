#!/bin/bash

set +x
# Install xarray into the user-site for the *system* Python BEFORE
# activating the venv. The Sphinx documentation extension at
# doc/sphinx/pism_config.py is loaded by the apt-installed
# /usr/bin/sphinx-build, whose shebang hardcodes /usr/bin/python3 and
# wouldn't see a venv-only install. The venv was created with
# --system-site-packages, so it picks up the user-site copy too. Drop
# this once ckhrulev/pism-ubuntu is rebuilt to include python3-xarray.
/usr/bin/python3 -m pip install --user --break-system-packages --quiet xarray

# Activate the environment containing mpi4py that will be needed for testing:
source $HOME/local/pism/bin/activate
set -x

set -e
set -u

# Uses environment variables CC and CXX to select the compiler

# The Make target to build
target=${target:-all}

# Compile using this many jobs
N=${N:-4}

# CMake installation prefix
CMAKE_PREFIX=${CMAKE_PREFIX:-/usr}

# PISM's source directory in the container
source_dir=${source_dir:-$HOME/project}

# The build directory to use
build_dir=${build_dir:-/tmp/pism-build}

# Installation prefix for all prerequisites built from source
lib_prefix=${lib_prefix:-$HOME/local}

# PETSc directory to use
export PETSC_DIR=${petsc_dir:-${lib_prefix}/petsc}

# Build Python bindings?
python=${python:-YES}

# Installation prefix
install_dir=${install_dir:-$HOME/local/pism}

git config --global --add safe.directory ${source_dir}

${CMAKE_PREFIX}/bin/cmake -S ${source_dir} \
               -B ${build_dir} \
               -DCMAKE_BUILD_TYPE=Debug \
               -DCMAKE_CXX_FLAGS="-Werror" \
               -DCMAKE_C_FLAGS="-Werror" \
               -DCMAKE_EXE_LINKER_FLAGS="-fuse-ld=lld" \
               -DCMAKE_INSTALL_PREFIX=${install_dir} \
               -DCMAKE_MODULE_LINKER_FLAGS="-fuse-ld=lld" \
               -DCMAKE_PREFIX_PATH="${lib_prefix}" \
               -DCMAKE_SHARED_LINKER_FLAGS="-fuse-ld=lld" \
               -DPism_BUILD_ICEBIN=YES \
               -DPism_BUILD_PYTHON_BINDINGS=${python} \
               -DPism_PEDANTIC_WARNINGS=YES \
               -DPism_USE_PARALLEL_NETCDF4=YES \
               -DPism_USE_PNETCDF=YES \
               -DPism_USE_PROJ=YES \
               -DPism_USE_YAC=YES

make --no-print-directory -C ${build_dir} -j ${N} ${target}
