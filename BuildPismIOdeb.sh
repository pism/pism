#!/bin/bash

module purge
module add cmake/3.17.1-gcc-9.1.0
#module add git/2.19.0
module add gcc/9.1.0-gcc-7.1.0
#module add gcc/6.4.0
module add openmpi/2.0.2p1_hpcx-gcc64

export CC=mpicc 
export CXX=mpicxx
export FC=mpifort
export F77=mpifort
export CFLAGS="-g -O0 -march=native -mtune=native"
export CXXFLAGS="-g -O0 -march=native -mtune=native"
export PETSC_ARCH=arch-linux2-c-debug-omp2
export PETSC_DIR=/home/dkrz/k202136/petsc
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/sw/rhel6-x64/hdf5/hdf5-1.8.20-parallel-openmpi2-gcc64/lib
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/work/bk0993/k202136/netcdf-c-4.7.4-gcc64-openmpi2/lib
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/work/bk0993/k202136/pnetcdf-1.12.1-gcc64-openmpi2/lib
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/work/bk0993/k202136/pio-2.5.2-gcc-openmpi2/lib
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/work/bk0993/k202136/cdi-1.9.9-gcc-openmpi2-parallel/lib
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/work/bk0993/k202136/yaxt-0.9.0-gcc-openmpi2/lib
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/sw/rhel6-x64/util/udunits-2.2.26-gcc64/lib
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/sw/rhel6-x64/numerics/fftw-3.3.7-openmp-gcc64/lib
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/sw/rhel6-x64/numerics/gsl-2.5-gcc48/lib
export UDUNITS2_ROOT=/sw/rhel6-x64/util/udunits-2.2.26-gcc64
export UDUNITS2_INCLUDES=$UDUNITS2_ROOT/include
export UDUNITS2_LIBRARIES=$UDUNITS2_ROOT/lib
export FFTW_ROOT=/sw/rhel6-x64/numerics/fftw-3.3.7-openmp-gcc64
#module add texlive/2016
#module add valgrind/3.13.0-gcc64
#module add nco/4.7.5-gcc64 
#module add cdo/1.9.5-gcc64
#module add ncview/2.1.4-gcc48
#module add netcdf_c/4.3.2-gcc48
export OMP_NUM_THREADS=1

export PATH=/work/bk0993/k202136/pio-2.5.2-gcc-openmpi2/lib:$PATH
export PATH=/sw/rhel6-x64/numerics/gsl-2.5-gcc48/lib:$PATH

export PISM_INSTALL_PREFIX=/home/dkrz/k202136/pismCDI2/pism
export HDF5ROOT=/sw/rhel6-x64/hdf5/hdf5-1.8.20-parallel-openmpi2-gcc64
#export PATH=$HDF5ROOT/lib/:$PATH
export ZLIBROOT=/sw/rhel6-x64/sys/zlib-1.2.8
export CMAKE=cmake
mkdir build
cd build

$CMAKE -DCMAKE_FIND_ROOT_PATH="/work/bk0993/k202136/netcdf-c-4.7.4-gcc64-openmpi2;/work/bk0993/k202136/pnetcdf-1.12.1-gcc-openmpi2;/work/bk0993/k202136/pio-2.5.2-gcc-openmpi2;/work/bk0993/k202136/cdi-1.9.9-gcc-openmpi2-parallel;/work/bk0993/k202136/yaxt-0.9.0-gcc-openmpi2" -DCMAKE_VERBOSE_MAKEFILE:BOOL=ON -DCMAKE_INSTALL_PREFIX=$PISM_INSTALL_PREFIX -DPism_BUILD_EXTRA_EXECS:BOOL=OFF -DPETSC_EXECUTABLE_RUNS=ON -DPism_USE_JANSSON=NO -DPism_USE_PIO=YES -DPism_USE_PARALLEL_NETCDF4=YES -DPism_USE_PNETCDF=YES -DPism_USE_YAXT=YES -DPism_USE_CDI=YES -DPism_USE_CDIPIO=YES -DMPIEXEC=srun -DHDF5_C_LIBRARY_z=$ZLIBROOT/lib/libz.so -DHDF5_C_LIBRARY_hdf5=$HDF5ROOT/lib/libhdf5.so.10 -DHDF5_C_LIBRARY_hdf5_hl=$HDF5ROOT/lib/libhdf5_hl.so.10 --debug-trycompile ../.

make install
