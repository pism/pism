This directory contains scripts testing PISM's behavior in the "ice
slab warming reversibility" experiment created by Thomas Kleiner.

This setup is described in *Enthalpy benchmark experiments for
numerical ice sheet models* by **T. Kleiner, M. Rückamp, J. H.
Bondzio, and A. Humbert**, see *experiment A* in
www.the-cryosphere.net/9/217/2015/

In short: a 1000 m thick slab of ice rests on a flat bed and is
subject to Dirichlet boundary conditions at the top surface. The
bottom surface boundary conditions are set by PISM depending on basal
conditions.

The top surface B.C. as a function of time is a "boxcar" function
switching from "cold" to "warm" and then back to "cold" conditions.

The steady state corresponding to the initial as well as final "cold"
surface B.C. has a linear-in-depth enthalpy distribution that can be
found analytically.

We expect PISM to be able to "cool down" after a warm period,
returning to the original model state.
