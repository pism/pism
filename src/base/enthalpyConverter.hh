// Copyright (C) 2009-2011, 2013, 2014, 2015 Andreas Aschwanden, Ed Bueler and Constantine Khroulev
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

#ifndef __enthalpyConverter_hh
#define __enthalpyConverter_hh

namespace pism {

class Config;

//! Converts between specific enthalpy and temperature or liquid content.
/*!
  Use this way, for example within IceModel with Config config member:
  \code
  #include "enthalpyConverter.hh"

  EnthalpyConverter EC(&config);  // runs constructor; do after initialization of Config config
  ...
  for (...) {
  ...
  E_s = EC.enthalpy_cts(p);
  ... etc ...
  }   
  \endcode

  The three methods that get the enthalpy from temperatures and liquid
  fractions, namely enthalpy(), enthalpy_permissive() are more strict
  about error checking. They throw RuntimeError if their arguments are
  invalid.

  This class is documented by [\ref AschwandenBuelerKhroulevBlatter].
*/
class EnthalpyConverter {
public:
  EnthalpyConverter(const Config &config);
  virtual ~EnthalpyConverter();
  
  bool is_temperate(double E, double pressure) const;

  double temperature(double E, double pressure) const;
  double melting_temperature(double pressure) const;
  double pressure_adjusted_temperature(double E, double pressure) const;

  double water_fraction(double E, double pressure) const;

  double enthalpy(double T, double omega, double pressure) const;
  double enthalpy_cts(double pressure) const;
  double enthalpy_permissive(double T, double omega, double pressure) const;

  double c_from_T(double T) const;

  double pressure(double depth) const;

protected:
  virtual double enthalpy_permissive_impl(double T, double omega, double pressure) const;
  virtual double enthalpy_cts_impl(double pressure) const;
  virtual double c_from_T_impl(double T) const;
  virtual double enthalpy_impl(double T, double omega, double pressure) const;
  virtual double water_fraction_impl(double E, double pressure) const;
  virtual double melting_temperature_impl(double pressure) const;
  virtual double temperature_impl(double E, double pressure) const;
  virtual bool is_temperate_impl(double E, double pressure) const;

  void enthalpy_interval(double pressure, double &E_s, double &E_l) const;

  double m_T_melting, m_L, m_c_i, m_rho_i, m_g, m_p_air, m_beta, m_T_tolerance;
  double m_T_0;
  bool m_do_cold_ice_methods;
};


//! An EnthalpyConverter for use in verification tests.
/*!
  Treats ice at any temperature as cold (= zero liquid fraction).  Makes absolute
  temperature (in K) and enthalpy proportional:  \f$E = c_i (T - T_0)\f$. 

  The pressure dependence of the pressure-melting temperature is neglected.

  Note: Any instance of FlowLaw uses an EnthalpyConverter; this is
  the one used in cold mode verification code.
*/
class ColdEnthalpyConverter : public EnthalpyConverter {
public:
  ColdEnthalpyConverter(const Config &config);
  virtual ~ColdEnthalpyConverter();

protected:
  double enthalpy_permissive_impl(double T, double /*omega*/, double /*pressure*/) const;
  double enthalpy_impl(double T, double /*omega*/, double /*pressure*/) const;
  double water_fraction_impl(double /*E*/, double /*pressure*/) const;
  double melting_temperature_impl(double /*pressure*/) const;
  bool is_temperate_impl(double /*E*/, double /*pressure*/) const;
  double temperature_impl(double E, double /*pressure*/) const;
};

/*
class KirchoffEnthalpyConverter : public EnthalpyConverter {
public:
  KirchoffEnthalpyConverter(const Config &config);
  virtual ~KirchoffEnthalpyConverter();
private:
  double L(double T_pm);
};
*/

} // end of namespace pism

#endif // __enthalpyConverter_hh

