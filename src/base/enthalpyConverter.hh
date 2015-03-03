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
  
  bool is_temperate(double E, double P) const;

  double temperature(double E, double P) const;
  double melting_temperature(double P) const;
  double pressure_adjusted_temperature(double E, double P) const;

  double water_fraction(double E, double P) const;

  double enthalpy(double T, double omega, double P) const;
  double enthalpy_cts(double P) const;
  double enthalpy_liquid(double P) const;
  double enthalpy_permissive(double T, double omega, double P) const;

  double c_from_T(double T) const;

  //! @brief Latent heat of fusion of water as a function of
  //! pressure-melting temperature.
  double L(double T_m) const;

  double pressure(double depth) const;

protected:
  virtual double enthalpy_permissive_impl(double T, double omega, double P) const;
  virtual double enthalpy_cts_impl(double P) const;
  virtual double c_from_T_impl(double T) const;
  virtual double L_impl(double T_pm) const;
  virtual double enthalpy_impl(double T, double omega, double P) const;
  virtual double water_fraction_impl(double E, double P) const;
  virtual double melting_temperature_impl(double P) const;
  virtual double temperature_impl(double E, double P) const;
  virtual double enthalpy_liquid_impl(double P) const;

  virtual bool is_temperate_impl(double E, double P) const;

  //! melting temperature of pure water at atmospheric pressure
  double m_T_melting;
  //! latent heat of fusion of water at atmospheric pressure
  double m_L;
  //! specific heat capacity of ice
  double m_c_i;
  //! density of ice
  double m_rho_i;
  //! acceleration due to gravity
  double m_g;
  //! atmospheric pressure
  double m_p_air;
  //! beta in the Clausius-Clapeyron relation (@f$ \diff{T_m}{p} = - \beta @f$).
  double m_beta;
  //! temperature tolerance used in `is_temperate()` in cold ice mode
  double m_T_tolerance;
  //! reference temperature in the definition of ice enthalpy
  double m_T_0;
  //! @brief if cold ice methods are selected, use `is_temperate()`
  //! check based on temperature, not enthalpy
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
  double enthalpy_permissive_impl(double T, double omega, double P) const;
  double enthalpy_impl(double T, double omega, double P) const;
  double water_fraction_impl(double E, double P) const;
  double melting_temperature_impl(double P) const;
  bool is_temperate_impl(double E, double P) const;
  double temperature_impl(double E, double P) const;
};

//! @brief An enthalpy converter including pressure-dependence of the
//! latent heat of fusion of water. (To be used by Glint2.)
class KirchhoffEnthalpyConverter : public EnthalpyConverter {
public:
  KirchhoffEnthalpyConverter(const Config &config);
  virtual ~KirchhoffEnthalpyConverter();
protected:
  double L_impl(double T_m) const;
private:
  //! specific heat capacity of pure water
  double m_c_w;
};


} // end of namespace pism

#endif // __enthalpyConverter_hh

