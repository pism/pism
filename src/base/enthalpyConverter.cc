// Copyright (C) 2009-2015 Andreas Aschwanden, Ed Bueler and Constantine Khroulev
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
  m_beta        = config.get("beta_CC"); // K Pa-1
  m_c_i         = config.get("ice_specific_heat_capacity"); // J kg-1 K-1
  m_g           = config.get("standard_gravity"); // m s-2
  m_L           = config.get("water_latent_heat_fusion"); // J kg-1
  m_p_air       = config.get("surface_pressure"); // Pa
  m_rho_i       = config.get("ice_density"); // kg m-3
  m_T_melting   = config.get("water_melting_point_temperature"); // K  
  m_T_tolerance = config.get("cold_mode_is_temperate_ice_tolerance"); // K
  m_T_0         = config.get("enthalpy_converter_reference_temperature"); // K

  m_do_cold_ice_methods  = config.get_flag("do_cold_ice_methods");
}

EnthalpyConverter::~EnthalpyConverter() {
  // empty
}

bool EnthalpyConverter::is_temperate(double E, double pressure) const {
  return this->is_temperate_impl(E, pressure);
}

double EnthalpyConverter::temperature(double E, double pressure) const {
  return this->temperature_impl(E, pressure);
}


//! Get pressure in ice from depth below surface using the hydrostatic assumption.
/*! If \f$d\f$ is the depth then
     \f[ p = p_{\text{air}}  + \rho_i g d. \f]
Frequently \f$d\f$ is computed from the thickess minus a level in the ice,
something like `ice_thickness(i, j) - z[k]`.  The input depth to this routine is allowed to
be negative, representing a position above the surface of the ice.
 */
double EnthalpyConverter::pressure(double depth) const {
  if (depth > 0.0) {
    return m_p_air + m_rho_i * m_g * depth;
  } else {
    return m_p_air; // at or above surface of ice
  }
}

double EnthalpyConverter::c_from_T(double T) const {
  return this->c_from_T_impl(T);
}

double EnthalpyConverter::c_from_T_impl(double /*T*/) const {
  return m_c_i;
}

//! Get melting temperature from pressure p.
/*!
     \f[ T_m(p) = T_{melting} - \beta p. \f]
 */
double EnthalpyConverter::melting_temperature(double pressure) const {
  return this->melting_temperature_impl(pressure);
}

double EnthalpyConverter::melting_temperature_impl(double pressure) const {
  return m_T_melting - m_beta * pressure;
}


//! Get enthalpy E_s(p) at cold-temperate transition point from pressure p.
/*! Returns
     \f[ E_s(p) = c_i (T_m(p) - T_0), \f]
 */
double EnthalpyConverter::enthalpy_cts(double p) const {
  return this->enthalpy_cts_impl(p);
}

double EnthalpyConverter::enthalpy_cts_impl(double p) const {
  return m_c_i * (melting_temperature(p) - m_T_0);
}

//! @brief Compute the maximum allowed value of ice enthalpy
//! (corresponds to @f$ \omega = 1 @f$).
double EnthalpyConverter::enthalpy_liquid(double p) const {
  return enthalpy_cts(p) + m_L;
}

//! Determines if E >= E_s(p), that is, if the ice is at the pressure-melting point.
bool EnthalpyConverter::is_temperate_impl(double E, double p) const {
  if (m_do_cold_ice_methods) {
    return (pressure_adjusted_temperature(E, p) >= m_T_melting - m_T_tolerance);
  } else {
    return (E >= enthalpy_cts(p));
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
double EnthalpyConverter::temperature_impl(double E, double p) const {

#if (PISM_DEBUG==1)
  if (E >= enthalpy_liquid(p)) {
    throw RuntimeError::formatted("E=%f at p=%f equals or exceeds that of liquid water",
                                  E, p);
  }
#endif

  if (E < enthalpy_cts(p)) {
    return (E / m_c_i) + m_T_0;
  } else {
    return melting_temperature(p);
  }
}


//! Get pressure-adjusted ice temperature, in Kelvin, from enthalpy and pressure.
/*!
The pressure-adjusted temperature is:
     \f[ T_{pa}(E,p) = T(E,p) - T_m(p) + T_{melting}. \f]
 */
double EnthalpyConverter::pressure_adjusted_temperature(double E, double pressure) const {
  return temperature(E, pressure) - melting_temperature(pressure) + m_T_melting;
}


//! Get liquid water fraction from enthalpy and pressure.
/*!
  From [@ref AschwandenBuelerKhroulevBlatter],
   @f[
   \omega(E,p) =
   \begin{cases}
     0.0, & E \le E_s(p), \\
     (E-E_s(p)) / L, & E_s(p) < E < E_l(p).
   \end{cases}
   @f]

   We do not allow liquid water (i.e. water fraction @f$ \omega=1.0 @f$).
 */
double EnthalpyConverter::water_fraction(double E, double pressure) const {
  return this->water_fraction_impl(E, pressure);
}

double EnthalpyConverter::water_fraction_impl(double E, double pressure) const {

#if (PISM_DEBUG==1)
  if (E >= enthalpy_liquid(pressure)) {
    throw RuntimeError::formatted("E=%f and pressure=%f correspond to liquid water",
                                  E, pressure);
  }
#endif

  double E_s = enthalpy_cts(pressure);
  if (E <= E_s) {
    return 0.0;
  } else {
    return (E - E_s) / m_L;
  }
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
double EnthalpyConverter::enthalpy(double T, double omega, double p) const {
  return this->enthalpy_impl(T, omega, p);
}

double EnthalpyConverter::enthalpy_impl(double T, double omega, double p) const {
  const double T_m = melting_temperature(p);

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
    throw RuntimeError::formatted("T < T_m AND omega > 0 is contradictory; got T=%f, T_m=%f, omega=%f",
                                  T, T_m, omega);
  }
#endif

  if (T < T_m) {
    return m_c_i * (T - m_T_0);
  } else {
    return enthalpy_cts(p) + omega * m_L;
  }
}


//! Compute enthalpy more permissively than enthalpy().
/*! Computes enthalpy from absolute temperature, liquid water fraction, and
pressure.  Use this form of enthalpy() when outside sources (e.g. information
from a coupler) might generate a temperature above the pressure melting point or
generate cold ice with a positive water fraction.

Treats temperatures above pressure-melting point as \e at the pressure-melting
point.  Interprets contradictory case of \f$T < T_m(p)\f$ and \f$\omega > 0\f$
as cold ice, ignoring the water fraction (\f$\omega\f$) value.

Calls enthalpy(), which validates its inputs. 

Computes:
  @f[
  E =
  \begin{cases}
    E(T,0.0,p),         & T < T_m(p) \quad \text{and} \quad \omega \ge 0, \\
    E(T_m(p),\omega,p), & T \ge T_m(p) \quad \text{and} \quad \omega \ge 0,
  \end{cases}
  @f]
  but ensures @f$ 0 \le \omega \le 1 @f$ in second case.  Calls enthalpy() for
  @f$ E(T,\omega,p) @f$.
 */
double EnthalpyConverter::enthalpy_permissive(double T, double omega, double pressure) const {
  return this->enthalpy_permissive_impl(T, omega, pressure);
}

double EnthalpyConverter::enthalpy_permissive_impl(double T, double omega, double pressure) const {

  const double T_m = melting_temperature(pressure);

  if (T < T_m) {
    return enthalpy(T, 0.0, pressure);
  } else { // T >= T_m(pressure) replaced with T = T_m(pressure)
    return enthalpy(T_m, std::max(0.0, std::min(omega, 1.0)), pressure);
  }
}

ColdEnthalpyConverter::ColdEnthalpyConverter(const Config &config)
  : EnthalpyConverter(config) {
  m_do_cold_ice_methods = true;
}

ColdEnthalpyConverter::~ColdEnthalpyConverter() {
  // empty
}

double ColdEnthalpyConverter::enthalpy_permissive_impl(double T,
                                                       double /*omega*/,
                                                       double /*pressure*/) const {
  return m_c_i * (T - m_T_0);
}
/*! */
double ColdEnthalpyConverter::enthalpy_impl(double T, double /*omega*/,
                                            double /*pressure*/) const {
  return m_c_i * (T - m_T_0);
}

/*! */
double ColdEnthalpyConverter::water_fraction_impl(double /*E*/,
                                                  double /*pressure*/) const {
  return 0.0;
}

/*! */
double ColdEnthalpyConverter::melting_temperature_impl(double /*pressure*/) const {
  return m_T_melting;
}
/*! */
bool ColdEnthalpyConverter::is_temperate_impl(double /*E*/, double /*pressure*/) const {
  return false;
}
/*! */
double ColdEnthalpyConverter::temperature_impl(double E, double /*pressure*/) const {
  return (E / m_c_i) + m_T_0;
}


} // end of namespace pism
