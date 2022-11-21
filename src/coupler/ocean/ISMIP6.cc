// Copyright (C) 2008-2019, 2022 Ed Bueler, Constantine Khroulev, Ricarda Winkelmann,
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

#include "ISMIP6.hh"
#include "pism/util/pism_options.hh"
#include "pism/util/io/io_helpers.hh"
#include "pism/util/MaxTimestep.hh"
#include "pism/util/pism_utilities.hh"

#include "pism/util/ConfigInterface.hh"
#include "pism/util/IceGrid.hh"
#include "pism/util/Mask.hh"
#include "pism/util/Vars.hh"

#include "pism/util/Time.hh"
#include "pism/geometry/Geometry.hh"

#include "pism/coupler/util/options.hh"

namespace pism {
namespace ocean {


ISMIP6::ISMIP6(IceGrid::ConstPtr g)
  :  CompleteOceanModel(g) {

  m_shelf_base_temperature = allocate_shelf_base_temperature(g);
  m_shelf_base_mass_flux   = allocate_shelf_base_mass_flux(g);

  ForcingOptions opt(*m_grid->ctx(), "ocean.ismip6");

  {
    unsigned int buffer_size = static_cast<unsigned int>(m_config->get_number("input.forcing.buffer_size"));

    File file(m_grid->com, opt.filename, PISM_NETCDF3, PISM_READONLY);

    m_shelfbtemp = std::make_shared<array::Forcing>(m_grid,
                                                    file,
                                                    "shelfbtemp",
                                                    "", // no standard name
                                                    buffer_size,
                                                    opt.periodic,
                                                    LINEAR);

   m_salinity_ocean = std::make_shared<array::Forcing>(m_grid,
                                                 file,
                                                  "salinity_ocean",
                                                  "", // no standard name
                                                   buffer_size,
                                                  opt.periodic,
                                                  LINEAR);
  }

  m_shelfbtemp->set_attrs("climate_forcing",
                          "absolute temperature at ice shelf base",
                          "Kelvin", "Kelvin", "", 0);
  m_salinity_ocean->set_attrs("climate_forcing",
                              "ocean salinity",
                              "g/kg", "g/kg", "", 0);
}

void ISMIP6::init_impl(const Geometry &geometry) {

  m_log->message(2,
             "* Initializing the ISMIP6 ocean reading base of the shelf temperature\n");

  ForcingOptions opt(*m_grid->ctx(), "ocean.ismip6");

  m_shelfbtemp->init(opt.filename, opt.periodic);
  m_salinity_ocean->init(opt.filename, opt.periodic);

  // read time-independent data right away:
  if (m_shelfbtemp->buffer_size() == 1 and m_salinity_ocean->buffer_size() == 1) {
    update(geometry, m_grid->ctx()->time()->current(), 0); // dt is irrelevant
  }
}

void ISMIP6::update_impl(const Geometry &geometry, double t, double dt) {
  (void) t;
  (void) dt;

  //m_shelfbtemp->copy_from(m_input_model->shelfbtemp());

  m_shelfbtemp->update(t, dt);  // FLO
  m_salinity_ocean->update(t, dt);  // FLO

  m_shelfbtemp->average(t, dt);  // FLO
  m_salinity_ocean->average(t, dt);  // FLO

  m_shelf_base_temperature->copy_from(*m_shelfbtemp);
  // FLO m_shelf_base_mass_flux->copy_from(*m_shelfbmassflux);

  // Set shelf base temperature to the melting temperature at the base (depth within the
  // ice equal to ice thickness).
  // FLO melting_point_temperature(H, *m_shelf_base_temperature);

  mass_flux(geometry.ice_thickness,
            *m_shelfbtemp,
            *m_salinity_ocean,
            *m_shelf_base_mass_flux); // call to ISMIP6 quadratic parametrisation

  double
    ice_density   = m_config->get_number("constants.ice.density"),
    water_density = m_config->get_number("constants.sea_water.density"),
    g             = m_config->get_number("constants.standard_gravity");

  compute_average_water_column_pressure(geometry, ice_density, water_density, g,
                                        *m_water_column_pressure);
}


MaxTimestep ISMIP6::max_timestep_impl(double t) const {
  (void) t;

  return MaxTimestep("ocean ismip6");
}


const array::Scalar& ISMIP6::shelf_base_temperature_impl() const {
  return *m_shelf_base_temperature;
}

const array::Scalar& ISMIP6::shelf_base_mass_flux_impl() const {
  return *m_shelf_base_mass_flux;
}

//! \brief Computes mass flux in [kg m-2 s-1].
/*!
 * Assumes that mass flux is proportional to the shelf-base heat flux.
 */
void ISMIP6::mass_flux(const array::Scalar &ice_thickness, const array::Scalar &shelfbtemp, const array::Scalar &salinity_ocean, array::Scalar &result) const {
  const double
    //melt_factor       = m_config->get_number("ocean.pik_melt_factor"),
    L                 = m_config->get_number("constants.fresh_water.latent_heat_of_fusion"),
    sea_water_density = m_config->get_number("constants.sea_water.density"),
    ice_density       = m_config->get_number("constants.ice.density"),
    c_p_ocean         = 3974.0, // J/(K*kg), specific heat capacity of ocean mixed layer
    gamma_0           = 3.52e-4;   // m/s, 11100 m/yr from local_MeanAnt method (Jourdain et al. 2020)
    //ocean_salinity    = 35.0;   // g/kg //********************* TO DO: implement also variable salinity in T_f //
    // FLO T_ocean           = units::convert(m_sys, -1.7, "Celsius", "Kelvin"); //Default in PISM-PIK

  //FIXME: gamma_T should be a function of the friction velocity, not a const

  array::AccessScope list{&ice_thickness, &shelfbtemp, &salinity_ocean, &result};

  // NOLINTBEGIN(readability-magic-numbers)
  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    // compute T_f(i, j) according to beckmann_goosse03, which has the
    // meaning of the freezing temperature of the ocean water directly
    // under the shelf, (of salinity 35psu) [this is related to the
    // Pressure Melting Temperature, see beckmann_goosse03 eq. 2 for
    // details]
    double
      shelfbaseelev = - (ice_density / sea_water_density) * ice_thickness(i,j),
      T_f           = 273.15 + (0.0832 -0.0575 * salinity_ocean(i,j) + 7.59e-4 * shelfbaseelev);
    // add 273.15 to convert from Celsius to Kelvin

    // compute ocean_heat_flux according to beckmann_goosse03
    // positive, if T_oc > T_ice ==> heat flux FROM ocean TO ice
    // FLO double ocean_heat_flux = melt_factor * sea_water_density * c_p_ocean * gamma_T * (T_ocean - T_f); // in W/m^2

    double
      thermal_forcing = shelfbtemp(i,j) - T_f,
      ocean_heat_flux = pow((sea_water_density * c_p_ocean) / (L * ice_density),2)  * gamma_0 * pow(fmax(thermal_forcing,0),2); // in W/m^2

    // TODO: T_ocean -> field!

    // shelfbmassflux is positive if ice is freezing on; here it is always negative:
    // same sign as ocean_heat_flux (positive if massflux FROM ice TO ocean)
    // FLO result(i,j) = ocean_heat_flux / (L * ice_density); // m s-1
    result(i,j) =  ocean_heat_flux; // m s-1

    // convert from [m s-1] to [kg m-2 s-1]:
    result(i,j) *= ice_density;
  }
  // NOLINTEND(readability-magic-numbers)
}

} // end of namespace ocean
} // end of namespace pism
