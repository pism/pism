// Copyright (C) 2011 Ed Bueler
//
// This file is part of PISM.
//
// PISM is free software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the Free Software
// Foundation; either version 2 of the License, or (at your option) any later
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

#include <petsc.h>  // for PetscErrorPrintf, etc.
#include "pism_const.hh"
#include "enthalpyConverter.hh"


EnthalpyConverter::EnthalpyConverter(const NCConfigVariable &config) {
  beta  = config.get("beta_CC");                                 // K Pa-1
  c_i   = config.get("ice_specific_heat_capacity");              // J kg-1 K-1
  g     = config.get("standard_gravity");			 // m s-2
  L     = config.get("water_latent_heat_fusion");                // J kg-1
  p_air = config.get("surface_pressure");                        // Pa
  rho_i = config.get("ice_density");                             // kg m-3
  T_melting = config.get("water_melting_point_temperature");       // K  
  T_tol = config.get("cold_mode_is_temperate_ice_tolerance");    // K 
  T_0   = config.get("enthalpy_converter_reference_temperature");// K  

  do_cold_ice_methods  = config.get_flag("do_cold_ice_methods");
}


//! Simple view of state of EnthalpyConverter.  viewer==NULL sends to stdout.
PetscErrorCode EnthalpyConverter::viewConstants(PetscViewer viewer) const {
  PetscErrorCode ierr;

  PetscTruth iascii;
  if (!viewer) {
    ierr = PetscViewerASCIIGetStdout(PETSC_COMM_WORLD,&viewer); CHKERRQ(ierr);
  }
  ierr = PetscTypeCompare((PetscObject)viewer,PETSC_VIEWER_ASCII,&iascii); CHKERRQ(ierr);
  if (!iascii) { SETERRQ(1,"Only ASCII viewer for EnthalpyConverter\n"); }

  ierr = PetscViewerASCIIPrintf(viewer,
    "\n<showing EnthalpyConverter constants:\n"); CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,
      "   beta  = %12.5e (K Pa-1)\n",    beta); CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,
      "   c_i   = %12.5f (J kg-1 K-1)\n",c_i); CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,
      "   g     = %12.5f (m s-2)\n",     g); CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,
      "   L     = %12.5e (J kg-1)\n",    L); CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,
      "   p_air = %12.5e (Pa)\n",        p_air); CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,
      "   rho_i = %12.5f (kg m-3)\n",    rho_i); CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,
      "   T_melting = %12.5f (K)\n",      T_melting); CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,
      "   T_tol = %12.5f (K)\n",         T_tol); CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,
      "   T_0   = %12.5f (K)\n",         T_0); CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,
      "   do_cold_ice_methods = %s\n", (do_cold_ice_methods) ? "true" : "false");
      CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,
      ">\n"); CHKERRQ(ierr);
  return 0;
}


//! Get pressure in ice from depth below surface using the hydrostatic assumption.
/*! If \f$d\f$ is the depth then
     \f[ p = p_{\text{air}}  + \rho_i g d. \f]
Frequently \f$d\f$ is computed from the thickess minus a level in the ice, 
something like "H[i][j] - z[k]".  The input depth to this routine is allowed to
be negative, representing a position above the surface of the ice.
 */ 
double EnthalpyConverter::getPressureFromDepth(double depth) const {
  if (depth <= 0.0) { // at or above surface of ice
    return p_air;
  } else {
    return p_air + rho_i * g * depth;
  }
}


//! Get melting temperature from pressure p.
/*!
     \f[ T_m(p) = T_{melting} - \beta p. \f]
 */ 
double EnthalpyConverter::getMeltingTemp(double p) const {
  return T_melting - beta * p;
}


//! Get enthalpy E_s(p) at cold-temperate transition point from pressure p.
/*! Returns 
     \f[ E_s(p) = c_i (T_m(p) - T_0), \f]
 */
double EnthalpyConverter::getEnthalpyCTS(double p) const {
  return c_i * (getMeltingTemp(p) - T_0);
}


//! Get enthalpies E_s(p) and E_l(p) (endpoints of temperate ice enthalpy range) from pressure p.
/*! Ice at enthalpy \f$E\f$ is temperate if \f$E_s(p) < E < E_l(p)\f$:
     \f[ E_s(p) = c_i (T_m(p) - T_0), \f]
     \f[ E_l(p) = E_s(p) + L. \f]
 */
PetscErrorCode EnthalpyConverter::getEnthalpyInterval(
                       double p, double &E_s, double &E_l) const {
  E_s = getEnthalpyCTS(p);
  E_l = E_s + L;
  return 0;
}


//! Computes the ratio CTS = E / E_s(p).  The cold-temperate transition surface (CTS) is the level surface CTS = 1.
/*!
If \f$E\f$ and \f$E_s\f$ are the ice mixture enthalpy and the CTS enthalpy,
respectively, then
  \f[ CTS(E,p) = \frac{E}{E_s(p)}.\f]
The level set CTS = 1 is the actual CTS.  The value computed here is greater
than one for temperate ice and less than one for cold ice.  The output of this
routine is normally only used for postprocessing.
*/
double EnthalpyConverter::getCTS(double E, double p) const {
  const double E_s = getEnthalpyCTS(p);
  return E / E_s;
}


//! Determines if E >= E_s(p), that is, if the ice is at the pressure-melting point.
bool EnthalpyConverter::isTemperate(double E, double p) const {
  if (do_cold_ice_methods) {
      double T_pa;
      getPATemp(E, p, T_pa);
      return (T_pa >= T_melting - T_tol);
  } else
      return (E >= getEnthalpyCTS(p));
}


//! Get absolute (not pressure-adjusted) ice temperature (K) from enthalpy and pressure.
/*! From \ref AschwandenBuelerKhroulevBlatter,
     \f[ T= T(E,p) = \begin{cases} 
                       c_i^{-1} E + T_0,  &  E < E_s(p), \\
                       T_m(p),            &  E_s(p) \le E < E_l(p).
                     \end{cases} \f]

We do not allow liquid water (%i.e. water fraction \f$\omega=1.0\f$) so we
return an error code of 1 if \f$E \ge E_l(p)\f$.
 */
PetscErrorCode EnthalpyConverter::getAbsTemp(double E, double p, double &T) const {
  double E_s, E_l;
  PetscErrorCode ierr = getEnthalpyInterval(p, E_s, E_l); CHKERRQ(ierr);
  if (E >= E_l) { // enthalpy equals or exceeds that of liquid water
    T = getMeltingTemp(p);
    return 1;
  }
  if (E < E_s) {
    T = (E / c_i) + T_0;
  } else {
    T = getMeltingTemp(p);
  }
  return 0;
}


//! Get pressure-adjusted ice temperature, in Kelvin, from enthalpy and pressure.
/*!
The pressure-adjusted temperature is:
     \f[ T_{pa}(E,p) = T(E,p) - T_m(p) + T_{melting}. \f]
 */
PetscErrorCode EnthalpyConverter::getPATemp(double E, double p, double &T_pa) const {
  PetscErrorCode ierr;
  double T = 0;	 // initialized to avoid a compiler warning
  ierr = getAbsTemp(E,p,T); CHKERRQ(ierr);
  T_pa = T - getMeltingTemp(p) + T_melting;
  return 0;
}


//! Get liquid water fraction from enthalpy and pressure.
/*!
From \ref AschwandenBuelerKhroulevBlatter,
   \f[ \omega(E,p) = \begin{cases}  0.0,            & E \le E_s(p), \\
                                    (E-E_s(p)) / L, & E_s(p) < E < E_l(p).
                     \end{cases} \f]

We do not allow liquid water (i.e. water fraction \f$\omega=1.0\f$) so we return
error code 1 if \f$E \ge E_l(p)\f$, but we still compute \f$\omega=1.0\f$.
 */
PetscErrorCode EnthalpyConverter::getWaterFraction(double E, double p, double &omega) const {
  double E_s, E_l;
  PetscErrorCode ierr = getEnthalpyInterval(p, E_s, E_l); CHKERRQ(ierr);
  if (E >= E_l) {
    omega = 1.0;
    return 1;
  }
  if (E <= E_s) {
    omega = 0.0;
  } else {
    omega = (E - E_s) / L;
  }
  return 0;
}


//! Is the ice mixture actually liquid water?
bool EnthalpyConverter::isLiquified(double E, double p) const {
  double E_s, E_l;
  getEnthalpyInterval(p, E_s, E_l);
  return (E >= E_l);
}


//! Compute enthalpy from absolute temperature, liquid water fraction, and pressure.
/*! This is an inverse function to the functions \f$T(E,p)\f$ and
\f$\omega(E,p)\f$ [\ref AschwandenBuelerKhroulevBlatter].  It returns:
  \f[E(T,\omega,p) =
       \begin{cases}
         c_i (T - T_0),     & T < T_m(p) \quad\text{and}\quad \omega = 0, \\
         E_s(p) + \omega L, & T = T_m(p) \quad\text{and}\quad \omega \ge 0.
       \end{cases} \f]
Certain cases are not allowed and return errors:
- \f$T<=0\f$ (error code 1)
- \f$\omega < 0\f$ or \f$\omega > 1\f$ (error code 2)
- \f$T>T_m(p)\f$ (error code 3)
- \f$T<T_m(p)\f$ and \f$\omega > 0\f$ (error code 4)
These inequalities may be violated in the sixth digit or so, however.
 */
PetscErrorCode EnthalpyConverter::getEnth(
                  double T, double omega, double p, double &E) const {
  const double T_m = getMeltingTemp(p);
  if (T <= 0.0) {
    SETERRQ1(1,"\n\nT = %f <= 0 is not a valid absolute temperature\n\n",T);
  }
  if ((omega < 0.0 - 1.0e-6) || (1.0 + 1.0e-6 < omega)) {
    SETERRQ1(2,"\n\nwater fraction omega=%f not in range [0,1]\n\n",omega);
  }
  if (T > T_m + 1.0e-6) {
    SETERRQ2(3,"T=%f exceeds T_m=%f; not allowed\n\n",T,T_m);
  }
  if ((T < T_m - 1.0e-6) && (omega > 0.0 + 1.0e-6)) {
    SETERRQ3(4,"T < T_m AND omega > 0 is contradictory\n\n",T,T_m,omega);
  }
  if (T < T_m) {
    E = c_i * (T - T_0);
  } else {
    E = getEnthalpyCTS(p) + omega * L;
  }
  return 0;
}


//! Compute enthalpy more permissively than getEnth().
/*! Computes enthalpy from absolute temperature, liquid water fraction, and
pressure.  Use this form of getEnth() when outside sources (e.g. information
from a coupler) might generate a temperature above the pressure melting point or
generate cold ice with a positive water fraction.

Treats temperatures above pressure-melting point as \e at the pressure-melting
point.  Interprets contradictory case of \f$T < T_m(p)\f$ and \f$\omega > 0\f$
as cold ice, ignoring the water fraction (\f$\omega\f$) value.

Checks if \f$T <= 0\f$ and returns error code 1 if so.

Computes:
  \f[E = \begin{cases}
            E(T,0.0,p),         & T < T_m(p) \quad \text{and} \quad \omega \ge 0, \\
            E(T_m(p),\omega,p), & T \ge T_m(p) \quad \text{and} \quad \omega \ge 0, 
         \end{cases} \f]
but ensures \f$0\le \omega \le 1\f$ in second case.  Calls getEnth() for
\f$E(T,\omega,p)\f$.
 */
PetscErrorCode EnthalpyConverter::getEnthPermissive(
                  double T, double omega, double p, double &E) const {
  PetscErrorCode ierr;
  if (T <= 0.0) {
    SETERRQ1(1,"\n\nT = %f <= 0 is not a valid absolute temperature\n\n",T);
  }
  const double T_m = getMeltingTemp(p);
  if (T < T_m) {
    ierr = getEnth(T, 0.0, p, E); CHKERRQ(ierr);
  } else { // T >= T_m(p) replaced with T = T_m(p)
    ierr = getEnth(T_m, PetscMax(0.0,PetscMin(omega,1.0)), p, E); CHKERRQ(ierr);
  }
  return 0;
}


//! Returns enthalpy for temperate ice with a given liquid fraction.
/*! Computes
  \f[E = E_s(p) + \omega L.\f]
Only the following case returns an error:
- \f$\omega < 0\f$ or \f$\omega > 1\f$ (error code 2)
These inequalities may be violated in the sixth digit or so, however.
 */
PetscErrorCode EnthalpyConverter::getEnthAtWaterFraction(
                        double omega, double p, double &E) const {
  if ((omega < 0.0 - 1.0e-6) || (1.0 + 1.0e-6 < omega)) {
    SETERRQ1(2,"\n\nwater fraction omega=%f not in range [0,1]\n\n",omega);
  }
  E = getEnthalpyCTS(p) + omega * L;
  return 0;
}

