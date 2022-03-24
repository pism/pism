#!/bin/bash

module purge
module load pism/stable1.0
export PETSC_DIR=/home/albrecht/software/petsc-3.9.1
#export PETSC_ARCH=linux-intel-64bit

export PISM_INSTALL_PREFIX=${PWD}
echo $PISM_INSTALL_PREFIX
export PATH=$HDF5ROOT/lib/:$PATH
export PATH=/p/system/packages/intel/parallel_studio_xe_2018_update1/compilers_and_libraries_2018.1.163/linux/compiler/lib/intel64:$PATH

export PART='' #srun --partition=broadwell --exclusive --ntasks=1'

mkdir build
cd build
$PART cmake -DCMAKE_VERBOSE_MAKEFILE:BOOL=ON \
      -DCMAKE_INSTALL_PREFIX=$PISM_INSTALL_PREFIX \
      -DPism_BUILD_EXTRA_EXECS:BOOL=OFF \
      -DPETSC_EXECUTABLE_RUNS=ON \
      -DPism_USE_JANSSON=NO \
      -DPism_USE_PARALLEL_NETCDF4=YES \
      -DPism_USE_PROJ4=YES \
      -DMPIEXEC=/p/system/slurm/bin/srun \
      -DHDF5_C_LIBRARY_z=$ZLIBROOT/lib/libz.so \
      -DHDF5_C_LIBRARY_mkl_rt=$MKLROOT/lib/intel64 \
      --debug-trycompile ../.
$PART make install
