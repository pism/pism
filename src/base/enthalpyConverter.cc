// Copyright (C) 2009 Andreas Aschwanden and Ed Bueler
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

#include <petsc.h>  // for PetscPrintf, PetscEnd
#include "enthalpyConverter.hh"

/*
*PREVIOUSLY* (r719 and earlier): From draft of \ref AschwandenBlatter2009,
     \f[ H_s(p) = -L + H_l(p), \f]
     \f[ H_l(p) = c_w T_m(p). \f]
NOW THE NORMALIZATION IS DIFFERENT
*/


EnthalpyConverter::EnthalpyConverter(NCConfigVariable *config) {
  if (config == NULL) {
    PetscPrintf(PETSC_COMM_WORLD,
      "\n\n\n  EnthalpyConverter ERROR in constructor:\n"
       "    config == NULL ... \n\n");
    PetscEnd();
  }
  beta  = config->get("beta_CC");                                 // K Pa-1
  c_b   = config->get("bedrock_thermal_specific_heat_capacity");  // J kg-1 K-1
  c_i   = config->get("ice_specific_heat_capacity");              // J kg-1 K-1
  g     = config->get("earth_gravity");                           // m s-2
  L     = config->get("water_latent_heat_fusion");                // J kg-1
  p_air = config->get("surface_pressure");                        // Pa
  rho_i = config->get("ice_density");                             // kg m-3
  T_0   = config->get("water_melting_temperature");               // K  
}


//! Simple view of state of EnthalpyConverter.  Give NULL to get stdout.
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
      "   c_b   = %12.5f (J kg-1 K-1)\n",c_b); CHKERRQ(ierr);
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
      "   T_0   = %12.5f (K)\n",         T_0); CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,
      ">\n"); CHKERRQ(ierr);
  return 0;
}


//! Get pressure in ice from depth below surface using the hydrostatic assumption.
/*! If \f$d\f$ is the depth then
     \f[ p = p_{\text{air}}  + \rho_i g d. \f]
Frequently \f$d\f$ is computed from the thickess minus a level in the ice, 
something like "H[i][j] - z[k]"
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
     \f[ T_m(p) = T_0 - \beta p. \f]
 */ 
double EnthalpyConverter::getMeltingTemp(double p) const {
  return T_0 - beta * p;
}


//! Get enthalpy E_s(p) at cold-temperate transition point, which normalizes enthalpy, from pressure p.
/*! Returns 
     \f[ E_s(p) = c_i T_m(p), \f]
 */
double EnthalpyConverter::getEnthalpyCTS(double p) const {
  return c_i * getMeltingTemp(p);
}


//! Get enthalpies E_s(p) and E_l(p) (endpoints of temperate ice enthalpy range) from pressure p.
/*! Ice at enthalpy \f$E\f$ is temperate if \f$E_s(p) < E < E_l(p)\f$:
     \f[ E_s(p) = c_i T_m(p), \f]
     \f[ E_l(p) = H_s(p) + L. \f]
 */
void EnthalpyConverter::getEnthalpyInterval(double p, double &E_s, double &E_l) const {
  E_s = getEnthalpyCTS(p);
  E_l = E_s + L;
}


//! Get absolute ice temperature (K) from enthalpy and pressure.
/*! From \ref AschwandenBlatter2009, equation (12)
     \f[ T= T(E,p) = \begin{cases} 
                       c_i^{-1} (E-E_s(p)) + T_m(p), & E < E_s(p), \\
                       T_m(p), &                       E_s(p) \le E < E_l(p).
                     \end{cases} \f]
But the first case simplifies if we expand \f$E_s\f$:
     \f[ c_i^{-1} (E-E_s(p)) + T_m(p) = c_i^{-1} (E-c_i T_m(p)) + T_m(p) = c_i^{-1} E.\f]
We do not allow liquid water (i.e. water fraction \f$\omega=1.0\f$) so we fail if
\f$E \ge E_l(p)\f$.
 */
double EnthalpyConverter::getAbsTemp(double E, double p) const {
  double E_s, E_l;
  getEnthalpyInterval(p, E_s, E_l);
  if (E < E_s) {
    return E / c_i;
  } else if (E < E_l) { // two cases in (12)
    return getMeltingTemp(p);
  } else {
    PetscPrintf(PETSC_COMM_WORLD,
      "\n\n\n  EnthalpyConverter ERROR in getAbsTemp():\n"
            "    enthalpy E=%f J kg-1 equals or exceeds that of liquid water (E_l=%f); ending ... \n\n",
      E,E_l);
    PetscEnd();
    return 0.0;
  }
}


//! Get pressure-adjusted ice temperature (K) from enthalpy and pressure.
/*! See getAbsTemp(), which computes \f$T(E,p)\f$, and getMeltingTemp(), which
computes \f$T_m(p)\f$.  We define the pressure-adjusted temperature to be:
     \f[ T_{pa} = T_{pa}(E,p) = T(E,p) - T_m(p) + T_0. \f]
 */
double EnthalpyConverter::getPATemp(double E, double p) const {
  return getAbsTemp(E,p) - getMeltingTemp(p) + T_0;
}


//! Get liquid water fraction from enthalpy and pressure.
/*! From \ref AschwandenBlatter2009, equation (12),
   \f[ \omega = \omega(E,p) = \begin{cases}
                                 0.0,            & E \le E_s(p), \\
                                 (E-E_s(p)) / L, & E_s(p) < E < E_l(p).
                              \end{cases} \f]

We do not allow liquid water (i.e. water fraction \f$\omega=1.0\f$) so we fail if
\f$E \ge E_l(p)\f$.
 */
double EnthalpyConverter::getWaterFraction(double E, double p) const {
  double E_s, E_l;
  getEnthalpyInterval(p, E_s, E_l);
  if (E <= E_s) { // two cases in (12)
    return 0.0;
  } else if (E < E_l) {
    return (E - E_s) / L;
  } else {
    PetscPrintf(PETSC_COMM_WORLD,
      "\n\n\n  EnthalpyConverter ERROR in getWaterFraction():\n"
            "    enthalpy equals or exceeds that of liquid water; ending ... \n\n");
    PetscEnd();
    return 1.0;
  }
}


//! Compute enthalpy from absolute temperature, liquid water fraction, and pressure.
/*! This is an inverse function to the functions \f$T(E,p)\f$ and \f$\omega(E,p)\f$.
It returns this enthalpy value:
  \f[E(T,\omega,p) = \begin{cases}
                        E_s(p) + c_i (T-T_m(p)), & T \le T_m(p) \quad\text{and}\quad \omega = 0, \\
                        E_s(p) + \omega L, &       T = T_m(p) \quad\text{and}\quad \omega \ge 0.
                     \end{cases} \f]

Certain cases are not allowed and are stopped:
- \f$T>T_m(p)\f$ is not allowed,
- \f$T<T_m(p)\f$ and \f$\omega > 0\f$ is not allowed,
- \f$\omega < 0\f$ is not allowed, and
- \f$\omega > 1\f$ is not allowed.

Because of these not-allowed cases, the following expression is also valid:
  \f[E(T,\omega,p) = E_s(p) + c_i (T-T_m(p)) + \omega L.\f]
 */
double EnthalpyConverter::getEnth(double T, double omega, double p) const {
  if ((omega < 0.0 - 1.0e-6) || (1.0 + 1.0e-6 < omega)) {
    PetscPrintf(PETSC_COMM_WORLD,
      "\n\n\n  EnthalpyConverter ERROR in getEnth(): water fraction omega=%f not in range [0,1]; ending ... \n\n",
      omega);
    PetscEnd();
  }
  const double T_m = getMeltingTemp(p);
  if (T > T_m + 1.0e-6) {
    PetscPrintf(PETSC_COMM_WORLD,
      "\n\n\n  EnthalpyConverter ERROR in getEnth(): T=%f exceeds T_m=%f so we have liquid water; ending ... \n\n",
      T,T_m);
    PetscEnd();
  }
  if ((T < T_m - 1.0e-6) && (omega > 0.0 + 1.0e-6)) {
    PetscPrintf(PETSC_COMM_WORLD,
      "\n\n\n  EnthalpyConverter ERROR in getEnth(): T < T_m AND omega > 0 is contradictory;\n"
            "     here T=%f, T_m=%f, omega=%f; ending ... \n\n", T, T_m, omega);
    PetscEnd();
  }
  const double E_s = getEnthalpyCTS(p);
  return E_s + c_i * (T - T_m) + omega * L;
}


//! Compute enthalpy more permissively in terms of meaning of input temperature and water fraction.  Compare getEnth().
/*! Computes enthalpy from absolute temperature, liquid water fraction, and pressure as before.
Use this form of getEnth() when outside sources (e.g. information from a coupler) might generate
a temperature above the pressure melting point or cold ice with a positive water fraction.

Treats temperatures above pressure-melting point as \e at the pressure-melting point.
Interprets contradictory case of \f$T < T_m(p)\f$ and \f$\omega > 0\f$ \e as cold ice,
ignoring water fraction \f$\omega > 0\f$.  Computes:
  \f[E_{\text{permissive}}(T,\omega,p)
       = \begin{cases}
            E(T,0.0,p),         & T < T_m(p) \quad \text{and} \quad \omega \ge 0, \\
            E(T_m(p),\omega,p), & T \ge T_m(p) \quad \text{and} \quad \omega \ge 0, 
         \end{cases} \f]
Calls getEnth() for \f$E(T,\omega,p)\f$.
 */
double EnthalpyConverter::getEnthPermissive(double T, double omega, double p) const {
  const double T_m = getMeltingTemp(p);
  if (T <= T_m) {
    return getEnth(T, 0.0, p);
  } else { // T >= T_m(p) replaced with T = T_m(p)
    return getEnth(T_m, omega, p);
  }
}


//! Compute enthalpy from absolute temperature within the bedrock using conditions at z=0 (bottom of ice and top of the bedrock) to set enthalpy scale.
/*! If \f$c_b\f$ is the specific heat capacity of the bedrock, if \f$E(z=0)\f$ is
the enthalpy at the bottom of the ice, and if \f$T(z=0)\f$ is the temperature at
the top of the bedrock, and if \f$z<0\f$ then we compute and return
  \f[ E_b = E(z=0) + c_b (T_b - T(z=0)) \f]
where \f$T_b\f$ is the temperature at level \f$z<0\f$ in the bedrock.

Input temperatures are assumed to be absolute (not pressure-adjusted).
Generally \f$T > T(z=0)\f$, in which case the resulting enthalpy
value comes out higher than \f$E(z=0)\f$.  Enthalpy scales linearly
with the temperature as we descend into the bedrock.
 */
double EnthalpyConverter::getEnthBedrock(double E_level_zero, double T_level_zero, double Tb) const {
  return E_level_zero + c_b * (Tb - T_level_zero);
}


//! Inverse function from getEnthBedrock().
/*! In same notation as for getEnthBedrock(),
\f[ T_b = \frac{E_b - E(z=0)}{c_b} + T(z=0) \f]
 */
double EnthalpyConverter::getAbsTempBedrock(double E_level_zero, double T_level_zero, double Eb) const {
  return ((Eb - E_level_zero) / c_b) + T_level_zero;
}

