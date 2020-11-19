#!/bin/bash

module purge
module add cmake/3.11.4
module add git/2.19.0
module add intel/18.0.4
module add intelmpi/2018.5.288
#module add inteltools/2018


export CC=mpiicc 
export CXX=mpiicpc
export FC=mpiifort
export F77=mpiifort
export CFLAGS="-fp-model precise -g -xHost -mtune=broadwell -gxx-name=/sw/rhel6-x64/gcc/gcc-6.4.0/bin/g++ -gcc-name=/sw/rhel6-x64/gcc/gcc-6.4.0/bin/gcc"
export CXXFLAGS="-fp-model precise -g -xHost -mtune=broadwell -gxx-name=/sw/rhel6-x64/gcc/gcc-6.4.0/bin/g++ -gcc-name=/sw/rhel6-x64/gcc/gcc-6.4.0/bin/gcc"
#export PETSC_ARCH=arch-linux2-c-opt
export PETSC_DIR=/sw/rhel6-x64/numerics/PETSc-3.12.2-impi2018-intel18
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/work/bk0993/k202136/netcdf-c-4.7.4-intel18-intelmpi2018/lib
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/work/bk0993/k202136/pnetcdf-1.12.1-intel18-intelmpi2018/lib
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/work/bk0993/k202136/pio-2.5.2-intel18-intelmpi2018/lib
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/work/bk0993/k202136/cdi-1.9.9-intel18-intelmpi2018-parallel/lib
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/work/bk0993/k202136/yaxt-0.9.0-intel18-intelmpi2018/lib
#export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/dkrz/k202136/echam6-PalMod/lib
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/sw/rhel6-x64/util/udunits-2.2.17-intel15/lib
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/sw/rhel6-x64/numerics/fftw-3.3.7-openmp-gcc64/lib
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/sw/rhel6-x64/numerics/gsl-2.5-gcc48/lib
export UDUNITS2_ROOT=/sw/rhel6-x64/util/udunits-2.2.17-intel15
export UDUNITS2_INCLUDES=$UDUNITS2_ROOT/include
export UDUNITS2_LIBRARIES=$UDUNITS2_ROOT/lib
export FFTW_ROOT=/sw/rhel6-x64/numerics/fftw-3.3.7-openmp-gcc64
module add texlive/2016
module add valgrind/3.13.0-gcc64
module add nco/4.7.5-gcc64 
module add cdo/1.9.5-gcc64
module add ncview/2.1.4-gcc48
module add netcdf_c/4.3.2-gcc48
export OMP_NUM_THREADS=1

export LD_LIBRARY_PATH=/sw/rhel6-x64/hdf5/hdf5-1.8.18-parallel-impi2017-intel14/lib:$LD_LIBRARY_PATH
export PISM_INSTALL_PREFIX=/home/dkrz/k202136/pismCDI2/pism
export HDF5ROOT=/sw/rhel6-x64/hdf5/hdf5-1.8.18-parallel-impi2017-intel14
export LD_LIBRARY_PATH=/sw/rhel6-x64/gcc/gcc-6.4.0/lib:/sw/rhel6-x64/gcc/gcc-6.4.0/lib64:$LD_LIBRARY_PATH
#export PATH=$HDF5ROOT/lib/:$PATH
export ZLIBROOT=/sw/rhel6-x64/sys/zlib-1.2.8
export CMAKE=cmake
mkdir build
cd build

$CMAKE -DCMAKE_FIND_ROOT_PATH="/work/bk0993/k202136/netcdf-c-4.7.4-intel18-intelmpi2018;/work/bk0993/k202136/pnetcdf-1.12.1-intel18-intelmpi2018;/work/bk0993/k202136/pio-2.5.2-intel18-intelmpi2018;/work/bk0993/k202136/cdi-1.9.9-intel18-intelmpi2018-parallel;/work/bk0993/k202136/yaxt-0.9.0-intel18-intelmpi2018" -DCMAKE_VERBOSE_MAKEFILE:BOOL=ON -DCMAKE_INSTALL_PREFIX=$PISM_INSTALL_PREFIX -DPism_BUILD_EXTRA_EXECS:BOOL=OFF -DPETSC_EXECUTABLE_RUNS=ON -DPism_USE_JANSSON=NO -DPism_USE_PIO=YES -DPism_USE_PARALLEL_NETCDF4=YES -DPism_USE_PNETCDF=YES -DPism_USE_YAXT=YES -DPism_USE_CDI=YES -DPism_USE_CDIPIO=YES -DMPIEXEC=srun -DHDF5_C_LIBRARY_z=$ZLIBROOT/lib/libz.so -DHDF5_C_LIBRARY_hdf5=$HDF5ROOT/lib/libhdf5.so -DHDF5_C_LIBRARY_hdf5_hl=$HDF5ROOT/lib/libhdf5_hl.so --debug-trycompile ../.

make install


#$CMAKE -DCMAKE_FIND_ROOT_PATH="/work/bk0993/k202136/netcdf-c-4.7.4-intel18-intelmpi2018;/work/bk0993/k202136/pnetcdf-1.12.1-intel18-intelmpi2018;/work/bk0993/k202136/pio-2.5.2-intel18-intelmpi2018" -DCMAKE_VERBOSE_MAKEFILE:BOOL=ON -DCMAKE_INSTALL_PREFIX=$PISM_INSTALL_PREFIX -DPism_BUILD_EXTRA_EXECS:BOOL=OFF -DPETSC_EXECUTABLE_RUNS=ON -DPism_USE_JANSSON=NO -DPism_USE_PIO=YES -DPism_USE_PARALLEL_NETCDF4=YES -DPism_USE_PNETCDF=YES -DMPIEXEC=srun -DHDF5_C_LIBRARY_z=$ZLIBROOT/lib/libz.so -DHDF5_C_LIBRARY_hdf5=$HDF5ROOT/lib/libhdf5.so.10 -DHDF5_C_LIBRARY_hdf5_hl=$HDF5ROOT/lib/libhdf5_hl.so.10 --debug-trycompile ../.

#make install
