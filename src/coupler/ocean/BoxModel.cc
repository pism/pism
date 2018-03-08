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
#include "Pico.hh"

namespace pism {
namespace ocean {


BoxModel::BoxModel(const Config &config) {

  // threshold between deep ocean and continental shelf
  continental_shelf_depth = config.get_double("ocean.pico.continental_shelf_depth");

  // value for ocean temperature around Antarctica if no other data available (cold conditions)
  T_dummy = -1.5 + config.get_double("constants.fresh_water.melting_point_temperature");
  S_dummy = 34.7; // value for ocean salinity around Antarctica if no other data available (cold conditions)

  earth_grav = config.get_double("constants.standard_gravity");
  rhoi       = config.get_double("constants.ice.density");
  rhow       = config.get_double("constants.sea_water.density");
  rho_star   = 1033;        // kg/m^3
  nu         = rhoi / rhow; // no unit

  latentHeat = config.get_double("constants.fresh_water.latent_heat_of_fusion"); //Joule / kg
  c_p_ocean  = 3974.0;                 // J/(K*kg), specific heat capacity of ocean mixed layer
  lambda     = latentHeat / c_p_ocean; // °C, NOTE K vs °C

  // Valus for linearized potential freezing point (from Xylar Asay-Davis, should be in Asay-Davis et al 2016, but not correct in there )
  a_pot = -0.0572;         // K/psu
  b_pot = 0.0788 + 273.15; // K
  c_pot = 7.77e-4;         // K/dbar

  // in-situ pressure melting point from Jenkins et al. 2010 paper
  a_in_situ = -0.0573;         // K/p_in_situu
  b_in_situ = 0.0832 + 273.15; // K
  c_in_situ = 7.53e-4;         // K/dbar

  // in-situ pressure melting point from Olbers & Hellmer 2010 paper
  // as          = -0.057;       // K/psu
  // bs          = 0.0832 + 273.15;       // K
  // cs          = 7.64e-4;      // K/dbar

  alpha = 7.5e-5; // 1/K
  beta  = 7.7e-4; // 1/psu

  // m s-1, best fit value in paper
  gamma_T = config.get_double("ocean.pico.heat_exchange_coefficent");
  // m6 kg−1 s−1, best fit value in paper
  overturning_coeff = config.get_double("ocean.pico.overturning_coefficent");

  // for shelf cells where normal box model is not calculated,
  // used in calculate_basal_melt_missing_cells(), compare POConstantPIK
  // m/s, thermal exchange velocity for Beckmann-Goose parameterization
  // this is the same meltFactor as in POConstantPIK
  meltFactor = config.get_double("ocean.pik_melt_factor");
}


double BoxModel::pressure(double ice_thickness) const {
      // pressure in dbar, 1dbar = 10000 Pa = 1e4 kg m-1 s-2
      return this->rhoi * this->earth_grav * ice_thickness * 1e-4;
}


double BoxModel::Toc_other_boxes(double area, double temp_in_boundary,
                                 double T_star, double overturning,
                                 double salinity_in_boundary) const {

  double g1 = area* this->gamma_T;
  double g2 = g1 / (this->nu * this->lambda);

  // temperature for Box i > 1
  return temp_in_boundary + g1 * T_star / (overturning + g1 - g2 *
                                           this->a_pot * salinity_in_boundary); // K
}

double BoxModel::Soc_box1(double Toc_box0, double Soc_box0, double Toc) const {

    return Soc_box0 - (Soc_box0 / (this->nu * this->lambda)) * (Toc_box0 - Toc);
}

double BoxModel::Soc_other_boxes(double salinity_in_boundary, double temperature_in_boundary, double Toc) const {

    return salinity_in_boundary - salinity_in_boundary * (temperature_in_boundary -
      Toc) / (this->nu * this->lambda); // psu;
}


//! equation 5 in the PICO paper.
//! calculate pressure melting point from potential temperature
double BoxModel::pot_pressure_melting(double salinity, double pressure) const {
    // using coefficients for potential temperature
    return this->a_pot * salinity + this->b_pot - this->c_pot * pressure;
}

//! equation 5 in the PICO paper.
//! calculate pressure melting point from in-situ temperature
double BoxModel::pressure_melting(double salinity, double pressure) const {
    // using coefficients for potential temperature
    return this->a_in_situ * salinity + this->b_in_situ - this->c_in_situ * pressure;
}

//! equation 8 in the PICO paper.
double BoxModel::bmelt_rate(double pm_point, double Toc) const {
    // in m/s
    return this->gamma_T / (this->nu * this->lambda) * (Toc - pm_point);
}

//! Beckmann & Goose meltrate
double BoxModel::bmelt_rate_beckm_goose(double Toc, double pot_pm_point) const {
      // in W/m^2
      double heatflux = this->meltFactor * this->rhow * this->c_p_ocean * this->gamma_T *
                        (Toc - pot_pm_point); // in W/m^2
      // in m s-1
      return heatflux / (this->latentHeat * this->rhoi);
}

//! equation 3 in the PICO paper. See also equation 4.
double BoxModel::overturning(double Soc_box0,
                             double Soc, double Toc_box0, double Toc) const {
 // in m^3/s
    return this->overturning_coeff * this->rho_star * (this->beta * (Soc_box0 - Soc) - this->alpha * (
        Toc_box0 - Toc));
}

  //! See equation A6 and lines before in PICO paper
double BoxModel::T_star(double salinity, double temperature, double pressure) const {
       // in Kelvin
      // FIXME: check that this stays always negative.
      // positive values are unphysical as colder temperatures
      // than pressure melting point would be ice, but we are in the ocean.
      // this should not occur as set_ocean_input_fields(...) sets
      // sets too cold temperatures to pressure melting point + 0.001
      return this->pot_pressure_melting(salinity, pressure) - temperature;
}

//! calculate p coefficent for solving the quadratic temperature equation
//! trough the p-q formula. See equation A12 in the PICO paper.
//! is only used once in f_Toc_box1(...)
double BoxModel::p_coeff(double g1, double s1) const {
    // in 1 / (1/K) = K
    // inputs g1 and overturning coefficicient are always positive
    // so output is positive if beta*s1 > alpha
    // which is shown in the text following equation A12
    return g1 / (this->overturning_coeff * this->rho_star * (this->beta * s1 - this->alpha));
}

//! calculate q coefficent for solving the quadratic temperature equation
//! trough the p-q formula. See equation A12 in the PICO paper.
//! is only used once in f_Toc_box1(...)
double BoxModel::q_coeff(double g1, double s1, double T_star) const {
    // in K / (1/K) = K^2
    return T_star * this->p_coeff(g1, s1);
}


double BoxModel::Toc_box1(double area, double T_star, double Soc_box0, double Toc_box0,
             bool *success) const  {

  double g1 = area * this->gamma_T;
  double s1 = Soc_box0 / (this->nu * this->lambda);

  // These are the coefficients for solving the quadratic temperature equation
  // trough the p-q formula.
  double p_coeff = this->p_coeff(g1, s1);
  double q_coeff = this->q_coeff(g1, s1, T_star);

  // This can only happen if T_star > 0.25*p_coeff, in particular T_star > 0
  // which can only happen for values of Toc_box0 close to the local pressure melting point
  if ((0.25 * p_coeff * p_coeff - q_coeff) < 0.0) {
    q_coeff = 0.25 * p_coeff * p_coeff;

    if (success) {
      *success = false;
    }
  }

  if (success) {
    *success = true;
  }

  // temperature for box 1, p-q formula
  // equation A12 in the PICO paper.
  return Toc_box0 - (-0.5 * p_coeff + sqrt(0.25 * PetscSqr(p_coeff) - q_coeff));
}

} // end of namespace ocean
} // end of namespace pism
