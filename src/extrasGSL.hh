// Copyright (C) 2004-2007 Jed Brown and Ed Bueler
//
// This file is part of Pism.
//
// Pism is free software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the Free Software
// Foundation; either version 2 of the License, or (at your option) any later
// version.
//
// Pism is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
// details.
//
// You should have received a copy of the GNU General Public License
// along with Pism; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

/* 
   extrasGSL.[hh|cc] is a collection of procedures which use the GSL library for numerical 
   computations.  These are useful to the ice model but are not dependent on any members of the
   class IceModel.  A number of the procedures here mimic Matlab functionality: 
      interp1_linear, dblquad_cubature.
   
   The cubature program by Steven Johnson is used for double integration, rather than the iterated
   use of a single-variable adaptive rule as in Matlab's actual dblquad.  cubature.[h|c] by Johnson
   is based on HIntLib by Rudolf Schuerer as well on GSL itself.  Note cubature is distributed with GPL.
   
   If the GSL library is not installed (in particular if one compiles with environment variable
   WITH_GSL=0), many of these codes return null results (of possibly reasonable kinds).

ELB 7/18/06
*/

double interp1_linear(const double x[], const double Y[], int N, double xi);

// get 'integrand' type
#include "cubature.h"
double dblquad_cubature(integrand f, const double ax, const double bx, const double ay, const double by,
                        double reqRelError, void *fdata);


struct ge_params { double dx; double dy; int p; int q; };

double ge_integrand(unsigned ndimMUSTBETWO, const double* xiANDeta, void* paramsIN);

struct vd_params { double t; double R0; double rk; double rho; double grav; double D; double eta; };

double viscDiscIntegrand (double kap, void * paramsIN);

double viscDisc(double t, double H0, double R0, double r, 
                double rho, double grav, double D, double eta);
