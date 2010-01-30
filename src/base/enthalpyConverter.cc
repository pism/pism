// Copyright (C) 2009-2010 Andreas Aschwanden, Ed Bueler and Constantine Khroulev
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
  T_0   = config.get("water_melting_temperature");               // K  
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
      "   T_0   = %12.5f (K)\n",         T_0); CHKERRQ(ierr);
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
     \f[ T_m(p) = T_0 - \beta p. \f]
 */ 
double EnthalpyConverter::getMeltingTemp(double p) const {
  return T_0 - beta * p;
}


//! Get enthalpy E_s(p) at cold-temperate transition point, which normalizes enthalpy, from pressure p.
/*! Returns 
     \f[ E_s(p) = c_i T_m(p), \f]
In particular,
     \f[ E_s( p_{\text{air}}) = c_i (T_0 - \beta p_{\text{air}}) = 548743.22\, \frac{\text{J}}{\text{kg}} \f]
is the enthalpy for surface ice.

Note that \f$p_{\text{air}} = 10^5\f$ Pa at temperature 
\f$T_0=273.15\f$ K, with standard choices \f$c_i=2009\f$ J kg-1 K-1 and
\f$\beta=7.53\times 10^{-8}\f$ K Pa-1.
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


//! Computes the ratio CTS = E / E_s(p).  The cold-temperate transition surface (CTS) is the level surface CTS = 1
/*!
If \f$E\f$ and \f$E_s\f$ are the enthalpy and the CTS enthalpy, respectively, then
  \f[ CTS = \frac{E}{E_s}.\f]
where CTS = 1 is the CTS contour line.  Thus CTS is greater than one for temperate ice
and less than one for cold ice.  Presumably only used for postprocessing.
*/ 
double EnthalpyConverter::getCTS(double E, double p) const {
  const double E_s = getEnthalpyCTS(p);
  return E / E_s;
}


//! Determines if E >= E_s(p), that is, if the ice is at the pressure-melting point.
bool EnthalpyConverter::isTemperate(double E, double p) const {
  const double E_s = getEnthalpyCTS(p);
  return (E >= E_s);
}


//! Get absolute ice temperature (K) from enthalpy and pressure.
/*! From \ref AschwandenBlatter, equation (12)
     \f[ T= T(E,p) = \begin{cases} 
                       c_i^{-1} (E-E_s(p)) + T_m(p), & E < E_s(p), \\
                       T_m(p), &                       E_s(p) \le E < E_l(p).
                     \end{cases} \f]
But the first case simplifies if we expand \f$E_s\f$:
     \f[ c_i^{-1} (E-E_s(p)) + T_m(p) = c_i^{-1} (E-c_i T_m(p)) + T_m(p) = c_i^{-1} E.\f]

We do not allow liquid water (i.e. water fraction \f$\omega=1.0\f$) so we fail if
\f$E \ge E_l(p)\f$.
 */
PetscErrorCode EnthalpyConverter::getAbsTemp(double E, double p, double &T) const {
  double E_s, E_l;
  getEnthalpyInterval(p, E_s, E_l);
  if (E < E_s) {
    T = E / c_i;
  } else if (E < E_l) { // two cases in (12)
    T = getMeltingTemp(p);
  } else {
    SETERRQ2(3,
      "\n\nenthalpy E=%f equals or exceeds that of liquid water (E_l=%f)\n\n",E,E_l);
  }
  return 0;
}


//! Get pressure-adjusted ice temperature (K) from enthalpy and pressure.
/*!
Calls getAbsTemp(), which computes \f$T(E,p)\f$, and getMeltingTemp(), which
computes \f$T_m(p)\f$.

We define the pressure-adjusted temperature to be:
     \f[ T_{pa} = T_{pa}(E,p) = T(E,p) - T_m(p) + T_0. \f]
 */
PetscErrorCode EnthalpyConverter::getPATemp(double E, double p, double &T_pa) const {
  PetscErrorCode ierr;
  const double T_m = getMeltingTemp(p);
  double T;
  ierr = getAbsTemp(E,p,T); CHKERRQ(ierr);
  T_pa = T - T_m + T_0;
  return 0;
}


//! Get liquid water fraction from enthalpy and pressure.
/*!
From \ref AschwandenBlatter, equation (12),
   \f[ \omega = \omega(E,p) = \begin{cases}
                                 0.0,            & E \le E_s(p), \\
                                 (E-E_s(p)) / L, & E_s(p) < E < E_l(p).
                              \end{cases} \f]

We do not allow liquid water (i.e. water fraction \f$\omega=1.0\f$) so we fail if
\f$E \ge E_l(p)\f$.
 */
PetscErrorCode EnthalpyConverter::getWaterFraction(double E, double p, double &omega) const {
  double E_s, E_l;
  getEnthalpyInterval(p, E_s, E_l);
  if (E <= E_s) { // two cases in (12)
    omega = 0.0;
  } else if (E < E_l) {
    omega = (E - E_s) / L;
  } else {
    SETERRQ1(2,"\n\nenthalpy E=%f equals or exceeds that of liquid water\n\n",E);
  }
  return 0;
}

//! Get liquid water fraction from enthalpy and pressure, but return omega=1 if high enthalpy.
/*! 
Same as getWaterFraction(), except if E > E_l(p) still succeeds and returns omega=1.

Useful in allowing the computation of (appropriately) high basal melt rates.
 */
double EnthalpyConverter::getWaterFractionLimited(double E, double p) const {
  double E_s, E_l;
  getEnthalpyInterval(p, E_s, E_l);
  if (E <= E_s) { // two cases in (12)
    return 0.0;
  } else if (E < E_l) {
    return (E - E_s) / L;
  } else {
    return 1.0;
  }
}


//! Compute enthalpy from absolute temperature, liquid water fraction, and pressure.
/*! This is an inverse function to the functions \f$T(E,p)\f$ and \f$\omega(E,p)\f$.
It returns this enthalpy value:
  \f[E(T,\omega,p) =
       \begin{cases}
         E_s(p) + c_i (T-T_m(p)), & T \le T_m(p) \quad\text{and}\quad \omega = 0, \\
         E_s(p) + \omega L, &       T = T_m(p) \quad\text{and}\quad \omega \ge 0.
       \end{cases} \f]

Certain cases are not allowed and return errors:
- \f$T>T_m(p)\f$ is not allowed,
- \f$T<T_m(p)\f$ and \f$\omega > 0\f$ is not allowed,
- \f$\omega < 0\f$ or \f$\omega > 1\f$ is not allowed.

Because of these not-allowed cases, the following expression is also valid:
  \f[E(T,\omega,p) = E_s(p) + c_i (T-T_m(p)) + \omega L.\f]
 */
PetscErrorCode EnthalpyConverter::getEnth(
                  double T, double omega, double p, double &E) const {
  if (T <= 0.0) {
    SETERRQ1(1,"\n\nT = %f <= 0 is not a valid absolute temperature\n\n",T);
  }
  if ((omega < 0.0 - 1.0e-6) || (1.0 + 1.0e-6 < omega)) {
    SETERRQ1(2,"\n\nwater fraction omega=%f not in range [0,1]\n\n",omega);
  }
  const double T_m = getMeltingTemp(p);
  if (T > T_m + 1.0e-6) {
    SETERRQ2(3,"T=%f exceeds T_m=%f so we have liquid water; not allowed\n\n",T,T_m);
  }
  if ((T < T_m - 1.0e-6) && (omega > 0.0 + 1.0e-6)) {
    SETERRQ3(4,"T < T_m AND omega > 0 is contradictory\n\n",T,T_m,omega);
  }
  const double E_s = getEnthalpyCTS(p);
  E = E_s + c_i * (T - T_m) + omega * L;
  return 0;
}


//! Compute enthalpy more permissively from temperature and water fraction.  Compare getEnth().
/*! Computes enthalpy from absolute temperature, liquid water fraction, and pressure as before.
Use this form of getEnth() when outside sources (e.g. information from a coupler) might generate
a temperature above the pressure melting point or cold ice with a positive water fraction.

Treats temperatures above pressure-melting point as \e at the pressure-melting point.
Interprets contradictory case of \f$T < T_m(p)\f$ and \f$\omega > 0\f$ \e as cold ice,
ignoring water fraction \f$\omega > 0\f$.  

Computes:
  \f[E = \begin{cases}
            E(T,0.0,p),         & T < T_m(p) \quad \text{and} \quad \omega \ge 0, \\
            E(T_m(p),\omega,p), & T \ge T_m(p) \quad \text{and} \quad \omega \ge 0, 
         \end{cases} \f]
but ensures \f$0\le \omega \le 1\f$ in second case.  Calls getEnth() for
\f$E(T,\omega,p)\f$.

Checks that \f$T > 0\f$ and returns error if not.
 */
PetscErrorCode EnthalpyConverter::getEnthPermissive(
                  double T, double omega, double p, double &E) const {
  PetscErrorCode ierr;
  if (T <= 0.0) {
    SETERRQ1(1,"\n\nT = %f <= 0 is not a valid absolute temperature\n\n",T);
  }
  const double T_m = getMeltingTemp(p);
  if (T <= T_m) {
    ierr = getEnth(T, 0.0, p, E); CHKERRQ(ierr);
  } else { // T >= T_m(p) replaced with T = T_m(p)
    ierr = getEnth(T_m, PetscMax(0.0,PetscMin(omega,1.0)), p, E); CHKERRQ(ierr);
  }
  return 0;
}


//! Returns enthalpy for temperate ice with a given liquid fraction.
/*! Computes
  \f[E = E_s(p) + \omega L.\f]
 */
PetscErrorCode EnthalpyConverter::getEnthAtWaterFraction(
                        double omega, double p, double &E) const {
  if ((omega < 0.0 - 1.0e-6) || (1.0 + 1.0e-6 < omega)) {
    SETERRQ1(2,"\n\nwater fraction omega=%f not in range [0,1]\n\n",omega);
  }
  const PetscScalar E_s = getEnthalpyCTS(p);
  E = E_s + omega * L;
  return 0;
}

