// Copyright (C) 2009-2014 Andreas Aschwanden, Ed Bueler and Constantine Khroulev
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

#include <petsc.h>  // for PetscErrorPrintf, etc.
#include "pism_const.hh"
#include "enthalpyConverter.hh"
#include "PISMConfig.hh"

#include "error_handling.hh"

namespace pism {

EnthalpyConverter::EnthalpyConverter(const Config &config) {
  beta      = config.get("beta_CC"); // K Pa-1
  c_i       = config.get("ice_specific_heat_capacity"); // J kg-1 K-1
  g         = config.get("standard_gravity"); // m s-2
  L         = config.get("water_latent_heat_fusion"); // J kg-1
  p_air     = config.get("surface_pressure"); // Pa
  rho_i     = config.get("ice_density"); // kg m-3
  T_melting = config.get("water_melting_point_temperature"); // K  
  T_tol     = config.get("cold_mode_is_temperate_ice_tolerance"); // K
  T_0       = config.get("enthalpy_converter_reference_temperature"); // K

  do_cold_ice_methods  = config.get_flag("do_cold_ice_methods");
}

//! Get pressure in ice from depth below surface using the hydrostatic assumption.
/*! If \f$d\f$ is the depth then
     \f[ p = p_{\text{air}}  + \rho_i g d. \f]
Frequently \f$d\f$ is computed from the thickess minus a level in the ice,
something like "H[i][j] - z[k]".  The input depth to this routine is allowed to
be negative, representing a position above the surface of the ice.
 */
double EnthalpyConverter::getPressureFromDepth(double depth) const {
  if (depth > 0.0) {
    return p_air + rho_i * g * depth;
  } else {
    return p_air; // at or above surface of ice
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
void EnthalpyConverter::getEnthalpyInterval(double p, double &E_s, double &E_l) const {
  E_s = getEnthalpyCTS(p);
  E_l = E_s + L;
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
  return E / getEnthalpyCTS(p);
}


//! Determines if E >= E_s(p), that is, if the ice is at the pressure-melting point.
bool EnthalpyConverter::isTemperate(double E, double p) const {
  if (do_cold_ice_methods) {
      return (getPATemp(E, p) >= T_melting - T_tol);
  } else {
      return (E >= getEnthalpyCTS(p));
  }
}


//! Get absolute (not pressure-adjusted) ice temperature (K) from enthalpy and pressure.
/*! From \ref AschwandenBuelerKhroulevBlatter,
     \f[ T= T(E,p) = \begin{cases} 
                       c_i^{-1} E + T_0,  &  E < E_s(p), \\
                       T_m(p),            &  E_s(p) \le E < E_l(p).
                     \end{cases} \f]

We do not allow liquid water (%i.e. water fraction \f$\omega=1.0\f$) so we
throw an exception if \f$E \ge E_l(p)\f$.
 */
double EnthalpyConverter::getAbsTemp(double E, double p) const {
  double E_s, E_l;
  getEnthalpyInterval(p, E_s, E_l);

  if (E >= E_l) {
    throw RuntimeError::formatted("E=%f at p=%f equals or exceeds that of liquid water",
                                  E, p);
  }

  if (E < E_s) {
    return (E / c_i) + T_0;
  } else {
    return getMeltingTemp(p);
  }
}


//! Get pressure-adjusted ice temperature, in Kelvin, from enthalpy and pressure.
/*!
The pressure-adjusted temperature is:
     \f[ T_{pa}(E,p) = T(E,p) - T_m(p) + T_{melting}. \f]
 */
double EnthalpyConverter::getPATemp(double E, double p) const {
  return getAbsTemp(E, p) - getMeltingTemp(p) + T_melting;
}


//! Get liquid water fraction from enthalpy and pressure.
/*!
From \ref AschwandenBuelerKhroulevBlatter,
   \f[ \omega(E,p) = \begin{cases}  0.0,            & E \le E_s(p), \\
                                    (E-E_s(p)) / L, & E_s(p) < E < E_l(p).
                     \end{cases} \f]

We do not allow liquid water (i.e. water fraction \f$\omega=1.0\f$).
 */
double EnthalpyConverter::getWaterFraction(double E, double p) const {
  double E_s, E_l;
  getEnthalpyInterval(p, E_s, E_l);
  if (E >= E_l) {
    throw RuntimeError::formatted("E=%f and p=%f correspond to liquid water", E, p);
  }
  if (E <= E_s) {
    return 0.0;
  } else {
    return (E - E_s) / L;
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
Certain cases are not allowed and throw exceptions:
- \f$T<=0\f$ (error code 1)
- \f$\omega < 0\f$ or \f$\omega > 1\f$ (error code 2)
- \f$T>T_m(p)\f$ (error code 3)
- \f$T<T_m(p)\f$ and \f$\omega > 0\f$ (error code 4)
These inequalities may be violated in the sixth digit or so, however.
 */
double EnthalpyConverter::getEnth(double T, double omega, double p) const {
  const double T_m = getMeltingTemp(p);

#if (PISM_DEBUG==1)
  if (T <= 0.0) {
    throw RuntimeError::formatted("T = %f <= 0 is not a valid absolute temperature",T);
  }
  if ((omega < 0.0 - 1.0e-6) || (1.0 + 1.0e-6 < omega)) {
    throw RuntimeError::formatted("water fraction omega=%f not in range [0,1]",omega);
  }
  if (T > T_m + 1.0e-6) {
    throw RuntimeError::formatted("T=%f exceeds T_m=%f; not allowed",T,T_m);
  }
  if ((T < T_m - 1.0e-6) && (omega > 0.0 + 1.0e-6)) {
    throw RuntimeError::formatted("T < T_m AND omega > 0 is contradictory",T,T_m,omega);
  }
#endif

  if (T < T_m) {
    return c_i * (T - T_0);
  } else {
    return getEnthalpyCTS(p) + omega * L;
  }
}


//! Compute enthalpy more permissively than getEnth().
/*! Computes enthalpy from absolute temperature, liquid water fraction, and
pressure.  Use this form of getEnth() when outside sources (e.g. information
from a coupler) might generate a temperature above the pressure melting point or
generate cold ice with a positive water fraction.

Treats temperatures above pressure-melting point as \e at the pressure-melting
point.  Interprets contradictory case of \f$T < T_m(p)\f$ and \f$\omega > 0\f$
as cold ice, ignoring the water fraction (\f$\omega\f$) value.

Checks if \f$T <= 0\f$ and throws an exception if so.

Computes:
  \f[E = \begin{cases}
            E(T,0.0,p),         & T < T_m(p) \quad \text{and} \quad \omega \ge 0, \\
            E(T_m(p),\omega,p), & T \ge T_m(p) \quad \text{and} \quad \omega \ge 0, 
         \end{cases} \f]
but ensures \f$0\le \omega \le 1\f$ in second case.  Calls getEnth() for
\f$E(T,\omega,p)\f$.
 */
double EnthalpyConverter::getEnthPermissive(double T, double omega, double p) const {
#if (PISM_DEBUG==1)
  if (T <= 0.0) {
    throw RuntimeError::formatted("T = %f <= 0 is not a valid absolute temperature", T);
  }
#endif

  const double T_m = getMeltingTemp(p);
  if (T < T_m) {
    return getEnth(T, 0.0, p);
  } else { // T >= T_m(p) replaced with T = T_m(p)
    return getEnth(T_m, std::max(0.0, std::min(omega, 1.0)), p);
  }
}


//! Returns enthalpy for temperate ice with a given liquid fraction.
/*! Computes
  \f[E = E_s(p) + \omega L.\f]
Only the following case returns an error:
- \f$\omega < 0\f$ or \f$\omega > 1\f$ (error code 2)
These inequalities may be violated in the sixth digit or so, however.
 */
double EnthalpyConverter::getEnthAtWaterFraction(double omega, double p) const {
#if (PISM_DEBUG==1)
  if ((omega < 0.0 - 1.0e-6) || (1.0 + 1.0e-6 < omega)) {
    throw RuntimeError::formatted("water fraction omega=%f not in range [0,1]",omega);
  }
#endif
  return getEnthalpyCTS(p) + omega * L;
}


} // end of namespace pism
