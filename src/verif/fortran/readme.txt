The Fortran 90 sources in this subdirectory,
   testFG.f90
   simpleFG.f90
are *not* built as part of PISM.  Instead, PISM uses the C sources 
   pism/src/exact/exactTestFG.h
   pism/src/exact/exactTestFG.c
   pism/src/exact/simpleFG.c
for computation of the exact solutions corresponding to "pismv -test F"
and "pismv -test G".

The F90 sources are added here to help other ice sheet modellers use the
thermomechanically coupled SIA exact solutions from 
   E. Bueler, J. Brown, and C. Lingle (2007).  "Exact solutions to the 
   thermomechanically coupled shallow ice approximation: effective tools
   for verification", J. Glaciol., to appear.
   [http://www.dms.uaf.edu/~bueler/exactcoupRESUB_TWOCOLA4.pdf]
for verification of other ice sheet simulations.
   
