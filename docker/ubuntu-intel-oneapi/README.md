### PISM container using Intel oneAPI compilers, MPI library and MKL

Based on Ubuntu 24.04 LTS.

Uses FFTW, GSL, libfyaml, PROJ and UDUNITS from Ubuntu packages.

HDF5, NetCDF (with parallel I/O support), YAXT, YAC and PETSc are built from sources using
Intel's C/C++ and Fortran compilers.

PETSc is built with 32 bit indexes, MUMPS, ScaLAPACK and Intel MKL as the BLAS/LAPACK
implementation. See `petsc.sh` for details.

PISM is built with PROJ, YAXT and YAC to be able to use flexible interpolation of gridded
inputs.
