// Copyright (C) 2009-2011, 2013, 2014, 2015, 2016, 2021, 2024 Andreas Aschwanden, Ed Bueler and Constantine Khroulev
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

#include <vector>
#include <memory>

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
  virtual ~EnthalpyConverter() = default;

  typedef std::shared_ptr<EnthalpyConverter> Ptr;

  bool is_temperate(double E, double P) const;
  bool is_temperate_relaxed(double E, double P) const;

  double temperature(double E, double P) const;
  double melting_temperature(double P) const;
  double pressure_adjusted_temperature(double E, double P) const;

  double water_fraction(double E, double P) const;

  double enthalpy(double T, double omega, double P) const;
  double enthalpy_cts(double P) const;
  double enthalpy_liquid(double P) const;
  double enthalpy_permissive(double T, double omega, double P) const;

  double c() const;
  double L(double T_pm) const;

  double pressure(double depth) const;
  void pressure(const std::vector<double> &depth,
                unsigned int ks, std::vector<double> &result) const;
protected:
  void validate_E_P(double E, double P) const;
  void validate_T_omega_P(double T, double omega, double P) const;

  double temperature_cold(double E) const;
  double enthalpy_cold(double T) const;

  //! melting temperature of pure water at atmospheric pressure
  double m_T_melting;
  //! latent heat of fusion of water at atmospheric pressure
  double m_L;
  //! specific heat capacity of ice
  double m_c_i;
  //! specific heat capacity of pure water
  double m_c_w;
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
  bool m_cold_mode;
};


//! An EnthalpyConverter for use in verification tests.

/*! Treats ice at any temperature below 10^6 kelvin as cold (= zero liquid fraction).

  The pressure dependence of the pressure-melting temperature is neglected.c;

  Note: Any instance of FlowLaw uses an EnthalpyConverter; this is
  the one used in cold mode verification code.


  This is the special enthalpy converter that is used in
  temperature-based verification tests only.

  In these tests ice temperatures in an exact solution may exceed the
  pressure-melting temperature, but we still want to pretend that this
  ice is "cold" to ensure that the map from enthalpy to temperature is
  one-to-one. (Normally enthalpy is mapped to the (temperature, water
  fraction) pair; here water fraction is zero, so enthalpy <-->
  (temperature, 0.)

  So, I had to pick a threshold (melting) temperature that is above
  all ice temperatures.  10^6K was chosen.
*/
class ColdEnthalpyConverter : public EnthalpyConverter {
public:
  ColdEnthalpyConverter(const Config &config);
  virtual ~ColdEnthalpyConverter() = default;
};

} // end of namespace pism

#endif // __enthalpyConverter_hh

