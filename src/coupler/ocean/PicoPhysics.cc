/* Copyright (C) 2018 PISM Authors
 *
 * This file is part of PISM.
 *
 * PISM is free software; you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation; either version 3 of the License, or (at your option) any later
 * version.
 *
 * PISM is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License
 * along with PISM; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 */
#include <cmath> // sqrt

#include "PicoPhysics.hh"

#include "pism/util/ConfigInterface.hh"

namespace pism {
namespace ocean {

PicoPhysics::PicoPhysics(const Config &config) {

  // threshold between deep ocean and continental shelf
  m_continental_shelf_depth = config.get_double("ocean.pico.continental_shelf_depth");

  // value for ocean temperature around Antarctica if no other data available (cold
  // conditions)
  m_T_dummy = -1.5 + config.get_double("constants.fresh_water.melting_point_temperature");

  // value for ocean salinity around Antarctica if no other data available (cold
  // conditions)
  m_S_dummy = 34.7;

  m_earth_grav        = config.get_double("constants.standard_gravity");
  m_ice_density       = config.get_double("constants.ice.density");
  m_sea_water_density = config.get_double("constants.sea_water.density");
  m_rho_star          = 1033;                                // kg/m^3
  m_nu                = m_ice_density / m_sea_water_density; // no unit

  //Joule / kg
  m_latentHeat = config.get_double("constants.fresh_water.latent_heat_of_fusion");

  // J/(K*kg), specific heat capacity of ocean mixed layer
  m_c_p_ocean = 3974.0;

  m_lambda = m_latentHeat / m_c_p_ocean; // °C, NOTE K vs °C

  // Values for linearized potential freezing point (from Xylar Asay-Davis, should be in
  // Asay-Davis et al 2016, but not correct in there )
  m_a_pot = -0.0572;         // K/psu
  m_b_pot = 0.0788 + 273.15; // K
  m_c_pot = 7.77e-4;         // K/dbar

  // in-situ pressure melting point from Jenkins et al. 2010 paper
  m_a_in_situ = -0.0573;         // K/p_in_situu
  m_b_in_situ = 0.0832 + 273.15; // K
  m_c_in_situ = 7.53e-4;         // K/dbar

  // in-situ pressure melting point from Olbers & Hellmer 2010 paper
  // as          = -0.057;       // K/psu
  // bs          = 0.0832 + 273.15;       // K
  // cs          = 7.64e-4;      // K/dbar

  m_alpha = 7.5e-5; // 1/K
  m_beta  = 7.7e-4; // 1/psu

  // m s-1, best fit value in paper
  m_gamma_T = config.get_double("ocean.pico.heat_exchange_coefficent");

  // m6 kg−1 s−1, best fit value in paper
  m_overturning_coeff = config.get_double("ocean.pico.overturning_coefficent");

  // for shelf cells where normal box model is not calculated, used in
  // calculate_basal_melt_missing_cells(), compare POConstantPIK m/s, thermal exchange
  // velocity for Beckmann-Goosse parameterization this is the same meltFactor as in
  // POConstantPIK
  m_meltFactor = config.get_double("ocean.pik_melt_factor");
}


double PicoPhysics::pressure(double ice_thickness) const {
  // pressure in dbar, 1dbar = 10000 Pa = 1e4 kg m-1 s-2
  return m_ice_density * m_earth_grav * ice_thickness * 1e-4;
}

TocBox1 PicoPhysics::Toc_box1(double area, double T_star, double Soc_box0, double Toc_box0) const {

  TocBox1 result = { false, 0.0 };

  const double
    g1 = area * m_gamma_T,
    s1 = Soc_box0 / (m_nu * m_lambda),
    p  = p_coeff(g1, s1),
    q  = p * T_star;

  // This can only happen if T_star > 0.25*p, in particular T_star > 0 which can only
  // happen for values of Toc_box0 close to the local pressure melting point
  double D = 0.25 * p * p - q;
  if (D < 0.0) {
    D             = 0.0;
    result.failed = true;
  }

  // temperature for box 1, p-q formula
  // equation A12 in the PICO paper.
  result.value = Toc_box0 - (-0.5 * p + sqrt(D));

  return result;
}

double PicoPhysics::Toc(double area, double temperature, double T_star, double overturning, double salinity) const {

  double g1 = area * m_gamma_T;
  double g2 = g1 / (m_nu * m_lambda);

  // temperature for Box i > 1
  return temperature + g1 * T_star / (overturning + g1 - g2 * m_a_pot * salinity); // K
}

double PicoPhysics::Soc_box1(double Toc_box0, double Soc_box0, double Toc) const {

  return Soc_box0 - (Soc_box0 / (m_nu * m_lambda)) * (Toc_box0 - Toc);
}

double PicoPhysics::Soc(double salinity, double temperature_in_boundary, double Toc) const {

  return salinity - salinity * (temperature_in_boundary - Toc) / (m_nu * m_lambda); // psu;
}


//! equation 5 in the PICO paper.
//! calculate pressure melting point from potential temperature
double PicoPhysics::theta_pm(double salinity, double pressure) const {
  // using coefficients for potential temperature
  return m_a_pot * salinity + m_b_pot - m_c_pot * pressure;
}

//! equation 5 in the PICO paper.
//! calculate pressure melting point from in-situ temperature
double PicoPhysics::T_pm(double salinity, double pressure) const {
  // using coefficients for in-situ temperature
  return m_a_in_situ * salinity + m_b_in_situ - m_c_in_situ * pressure;
}

//! equation 8 in the PICO paper.
double PicoPhysics::melt_rate(double pm_point, double Toc) const {
  // in m/s
  return m_gamma_T / (m_nu * m_lambda) * (Toc - pm_point);
}

//! Beckmann & Goosse meltrate
double PicoPhysics::melt_rate_beckmann_goosse(double pot_pm_point, double Toc) const {
  // in W/m^2
  double heat_flux = m_meltFactor * m_sea_water_density * m_c_p_ocean * m_gamma_T * (Toc - pot_pm_point);
  // in m s-1
  return heat_flux / (m_latentHeat * m_ice_density);
}

//! equation 3 in the PICO paper. See also equation 4.
double PicoPhysics::overturning(double Soc_box0, double Soc, double Toc_box0, double Toc) const {
  // in m^3/s
  return m_overturning_coeff * m_rho_star * (m_beta * (Soc_box0 - Soc) - m_alpha * (Toc_box0 - Toc));
}

//! See equation A6 and lines before in PICO paper
double PicoPhysics::T_star(double salinity, double temperature, double pressure) const {
  // in Kelvin
  // FIXME: check that this always stays negative.
  // positive values are unphysical as colder temperatures
  // than pressure melting point would be ice, but we are in the ocean.
  // this should not occur as set_ocean_input_fields(...) sets
  // sets too cold temperatures to pressure melting point + 0.001
  return theta_pm(salinity, pressure) - temperature;
}

//! calculate p coefficent for solving the quadratic temperature equation
//! trough the p-q formula. See equation A12 in the PICO paper.
//! is only used once in Toc_box1(...)
double PicoPhysics::p_coeff(double g1, double s1) const {
  // in 1 / (1/K) = K
  // inputs g1 and overturning coefficicient are always positive
  // so output is positive if beta*s1 > alpha
  // which is shown in the text following equation A12
  return g1 / (m_overturning_coeff * m_rho_star * (m_beta * s1 - m_alpha));
}

double PicoPhysics::gamma_T() const {
  return m_gamma_T;
}

double PicoPhysics::overturning_coeff() const {
  return m_overturning_coeff;
}

double PicoPhysics::T_dummy() const {
  return m_T_dummy;
}

double PicoPhysics::S_dummy() const {
  return m_S_dummy;
}

double PicoPhysics::ice_density() const {
  return m_ice_density;
}

double PicoPhysics::continental_shelf_depth() const {
  return m_continental_shelf_depth;
}

} // end of namespace ocean
} // end of namespace pism
