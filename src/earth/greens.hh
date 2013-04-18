// Copyright (C) 2007--2009 Ed Bueler
//
// This file is part of PISM.
//
// PISM is free software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the Free Software
// Foundation; either version 3 of the License, or (at your option) any later
// version.
//
// PISM is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
// details.
//
// You should have received a copy of the GNU General Public License
// along with PISM; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

#ifndef __greens_hh
#define __greens_hh

//! Parameters used to access elastic Green's function from the \ref Farrell earth model.
struct ge_params {
   double dx, dy;
   int    p, q; 
};

//! The integrand in defining the elastic Green's function from the \ref Farrell earth model.
/*!
For G^E(r), the Green's function of spherical layered elastic earth model.  From data in 
\ref LingleClark.  See also \ref BLKfastearth.
 */
double ge_integrand(unsigned ndimMUSTBETWO, const double* xiANDeta, void* paramsIN);


//! Parameters used to describe the response of the viscous half-space model to a disc load.
struct vd_params {
   double t, R0, rk, rho, grav, D, eta;
};

//! Integrand defining the response of the viscous half-space model to a disc load.
/*!
For the solution of the disc load case of the viscous half-space model, see 
appendix B of \ref BLK2006earth.  See also \ref LingleClark and \ref BLKfastearth.
 */
double vd_integrand (double kap, void * paramsIN);

//! Actually compute the response of the viscous half-space model in \ref LingleClark, to a disc load.
double viscDisc(double t, double H0, double R0, double r, 
                double rho, double grav, double D, double eta);

#endif	/* __greens_hh */

