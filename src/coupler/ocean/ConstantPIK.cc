// Copyright (C) 2008-2019, 2021, 2022, 2023, 2024 Ed Bueler, Constantine Khroulev, Ricarda Winkelmann,
// Gudfinna Adalgeirsdottir, Andy Aschwanden and Torsten Albrecht
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

#include "pism/coupler/ocean/ConstantPIK.hh"
#include "pism/util/ConfigInterface.hh"
#include "pism/util/Grid.hh"
#include "pism/util/MaxTimestep.hh"
#include "pism/geometry/Geometry.hh"

namespace pism {
namespace ocean {

PIK::PIK(std::shared_ptr<const Grid> g)
  : CompleteOceanModel(g) {
  // empty
}

void PIK::init_impl(const Geometry &geometry) {
  (void) geometry;

  m_log->message(2,
                 "* Initializing the constant (PIK) ocean model...\n");

  double
    ice_density   = m_config->get_number("constants.ice.density"),
    water_density = m_config->get_number("constants.sea_water.density"),
    g             = m_config->get_number("constants.standard_gravity");

  compute_average_water_column_pressure(geometry, ice_density, water_density, g,
                                           *m_water_column_pressure);
}

MaxTimestep PIK::max_timestep_impl(double t) const {
  (void) t;
  return MaxTimestep("ocean PIK");
}

void PIK::update_impl(const Geometry &geometry, double t, double dt) {
  (void) t;
  (void) dt;

  const array::Scalar &H = geometry.ice_thickness;

  // Set shelf base temperature to the melting temperature at the base (depth within the
  // ice equal to ice thickness).
  melting_point_temperature(H, *m_shelf_base_temperature);

  mass_flux(H, *m_shelf_base_mass_flux);

  double
    ice_density   = m_config->get_number("constants.ice.density"),
    water_density = m_config->get_number("constants.sea_water.density"),
    g             = m_config->get_number("constants.standard_gravity");

  compute_average_water_column_pressure(geometry, ice_density, water_density, g,
                                           *m_water_column_pressure);
}

/*!
 * Compute melting temperature at a given depth within the ice.
 */
void PIK::melting_point_temperature(const array::Scalar &depth,
                                    array::Scalar &result) const {
  const double
    T0          = m_config->get_number("constants.fresh_water.melting_point_temperature"), // K
    beta_CC     = m_config->get_number("constants.ice.beta_Clausius_Clapeyron"),
    g           = m_config->get_number("constants.standard_gravity"),
    ice_density = m_config->get_number("constants.ice.density");

  array::AccessScope list{&depth, &result};

  for (auto p = m_grid->points(); p; p.next()) {
    const int i = p.i(), j = p.j();
    const double pressure = ice_density * g * depth(i,j); // FIXME task #7297
    // result is set to melting point at depth
    result(i,j) = T0 - beta_CC * pressure;
  }
}

//! \brief Computes mass flux in [kg m-2 s-1].
/*!
 * Assumes that mass flux is proportional to the shelf-base heat flux.
 */
void PIK::mass_flux(const array::Scalar &ice_thickness, array::Scalar &result) const {
  const double
    melt_factor       = m_config->get_number("ocean.pik_melt_factor"),
    L                 = m_config->get_number("constants.fresh_water.latent_heat_of_fusion"),
    sea_water_density = m_config->get_number("constants.sea_water.density"),
    ice_density       = m_config->get_number("constants.ice.density"),
    c_p_ocean         = 3974.0, // J/(K*kg), specific heat capacity of ocean mixed layer
    gamma_T           = 1e-4,   // m/s, thermal exchange velocity
    ocean_salinity    = 35.0,   // g/kg
    T_ocean           = units::convert(m_sys, -1.7, "degree_Celsius", "kelvin"); //Default in PISM-PIK

  //FIXME: gamma_T should be a function of the friction velocity, not a const

  array::AccessScope list{&ice_thickness, &result};

  for (auto p = m_grid->points(); p; p.next()) {
    const int i = p.i(), j = p.j();

    // compute T_f(i, j) according to beckmann_goosse03, which has the
    // meaning of the freezing temperature of the ocean water directly
    // under the shelf, (of salinity 35psu) [this is related to the
    // Pressure Melting Temperature, see beckmann_goosse03 eq. 2 for
    // details]
    double
      shelfbaseelev = - (ice_density / sea_water_density) * ice_thickness(i,j),
      T_f           = 273.15 + (0.0939 -0.057 * ocean_salinity + 7.64e-4 * shelfbaseelev);
    // add 273.15 to convert from Celsius to kelvin

    // compute ocean_heat_flux according to beckmann_goosse03
    // positive, if T_oc > T_ice ==> heat flux FROM ocean TO ice
    double ocean_heat_flux = melt_factor * sea_water_density * c_p_ocean * gamma_T * (T_ocean - T_f); // in W/m^2

    // TODO: T_ocean -> field!

    // shelfbmassflux is positive if ice is freezing on; here it is always negative:
    // same sign as ocean_heat_flux (positive if massflux FROM ice TO ocean)
    result(i,j) = ocean_heat_flux / (L * ice_density); // m s-1

    // convert from [m s-1] to [kg m-2 s-1]:
    result(i,j) *= ice_density;
  }
}

} // end of namespace ocean
} // end of namespace pism
