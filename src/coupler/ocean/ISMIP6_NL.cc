// Copyright (C) 2008-2020, 2022, 2023 Ed Bueler, Constantine Khroulev, Ricarda Winkelmann,
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

#include "ISMIP6_NL.hh"

#include "pism/coupler/util/options.hh"
#include "pism/geometry/Geometry.hh"
#include "pism/util/ConfigInterface.hh"
#include "pism/util/Grid.hh"
#include "pism/util/MaxTimestep.hh"
#include "pism/util/Time.hh"
#include "pism/util/pism_utilities.hh"

/*!
 *
 * https://gmd.copernicus.org/articles/12/2255/2019/ (Favier)
 *
 * https://tc.copernicus.org/articles/14/3111/2020/ (Jourdain)
 *
 */

namespace pism {
namespace ocean {

ISMIP6nl::ISMIP6nl(std::shared_ptr<const Grid> g)
  :  CompleteOceanModel(g),
     m_basin_mask(m_grid, "basin_mask"),
     m_thermal_forcing(m_grid, "thermal_forcing") {

  m_shelf_base_temperature = allocate_shelf_base_temperature(g);
  m_shelf_base_mass_flux   = allocate_shelf_base_mass_flux(g);

  ForcingOptions opt(*m_grid->ctx(), "ocean.ismip6nl");

  {
    unsigned int buffer_size = static_cast<unsigned int>(m_config->get_number("input.forcing.buffer_size"));

    File file(m_grid->com, opt.filename, io::PISM_NETCDF3, io::PISM_READONLY);

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

  m_shelfbtemp->metadata(0).long_name("absolute temperature at ice shelf base").units("Kelvin");
  m_salinity_ocean->metadata(0).long_name("ocean salinity").units("g/kg");
  m_basin_mask.metadata(0).long_name("mask of drainage basins");
  m_n_basins = 0;
}

void ISMIP6nl::init_impl(const Geometry &geometry) {

  m_log->message(2,
                 "* Initializing the ISMIP6 ocean reading base of the shelf temperature\n");

  ForcingOptions opt(*m_grid->ctx(), "ocean.ismip6nl");

  m_shelfbtemp->init(opt.filename, opt.periodic);
  m_salinity_ocean->init(opt.filename, opt.periodic);

  m_basin_mask.regrid(opt.filename, io::Default::Nil());
  m_n_basins = static_cast<int>(array::max(m_basin_mask)) + 1; // Basins id starts at 0 in the input file

  // read time-independent data right away:
  if (m_shelfbtemp->buffer_size() == 1 and m_salinity_ocean->buffer_size() == 1) {
    update(geometry, m_grid->ctx()->time()->current(), 0); // dt is irrelevant
  }
}

void ISMIP6nl::update_impl(const Geometry &geometry, double t, double dt) {
  (void) t;
  (void) dt;

  m_shelfbtemp->update(t, dt);  // FLO
  m_salinity_ocean->update(t, dt);  // FLO

  m_shelfbtemp->average(t, dt);  // FLO
  m_salinity_ocean->average(t, dt);  // FLO

  m_shelf_base_temperature->copy_from(*m_shelfbtemp);

  compute_thermal_forcing(geometry.ice_thickness,
                          *m_shelfbtemp, *m_salinity_ocean, m_thermal_forcing);

  std::vector<double> basin_TF(m_n_basins);
  compute_avg_thermal_forcing(geometry.cell_type,
                              m_basin_mask,
                              m_thermal_forcing,
                              basin_TF); // per basin

  mass_flux(m_thermal_forcing,
            m_basin_mask,
            basin_TF,
            *m_shelf_base_mass_flux); // call to ISMIP6 quadratic parametrisation

  double
    ice_density   = m_config->get_number("constants.ice.density"),
    water_density = m_config->get_number("constants.sea_water.density"),
    g             = m_config->get_number("constants.standard_gravity");

  compute_average_water_column_pressure(geometry, ice_density, water_density, g,
                                        *m_water_column_pressure);
}

MaxTimestep ISMIP6nl::max_timestep_impl(double t) const {
  (void) t;
  return {"ocean ismip6nl"};
}

const array::Scalar& ISMIP6nl::shelf_base_temperature_impl() const {
  return *m_shelf_base_temperature;
}

const array::Scalar& ISMIP6nl::shelf_base_mass_flux_impl() const {
  return *m_shelf_base_mass_flux;
}

void ISMIP6nl::compute_thermal_forcing(const array::Scalar &ice_thickness,
                                       const array::Scalar &shelfbtemp,
                                       const array::Scalar &salinity_ocean,
                                       array::Scalar &result) {

  const double
    sea_water_density = m_config->get_number("constants.sea_water.density"),
    ice_density       = m_config->get_number("constants.ice.density");

  array::AccessScope list{&ice_thickness, &shelfbtemp, &salinity_ocean, &result};

  // NOLINTBEGIN(readability-magic-numbers)
  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    // compute T_f(i, j) according to Favier et al. XXXX and Jourdain et al. XXXX
    // UNESCO Seawater freezing temperature
    double
      shelfbaseelev = - (ice_density / sea_water_density) * ice_thickness(i,j),
      T_f           = 273.15 + (0.0832 -0.0575 * salinity_ocean(i,j) + 7.59e-4 * shelfbaseelev);
    // add 273.15 to convert from Celsius to Kelvin

    result(i,j) = shelfbtemp(i,j) - T_f;
  }
  // NOLINTEND(readability-magic-numbers)
}


// This routine should calculate the averaged thermal-forcing
// over the sum of ice shelves points contained in each drainage basins
void ISMIP6nl::compute_avg_thermal_forcing(const array::CellType &cell_type,
                                           const array::Scalar &basin_mask,
                                           const array::Scalar &thermal_forcing,
                                           std::vector<double> &result) {

  std::vector<int> basin_size_local(m_n_basins, 0);
  std::vector<double> thermal_forcing_sum(m_n_basins, 0.0);

  array::AccessScope list{&cell_type, &thermal_forcing, &basin_mask};

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    if (cell_type.floating_ice(i, j)) {
      int basin_id = basin_mask.as_int(i, j);

      basin_size_local[basin_id]++;
      thermal_forcing_sum[basin_id] += thermal_forcing(i, j);
    }
  }

  std::vector<int> basin_size(m_n_basins, 0);
  GlobalSum(m_grid->com, basin_size_local.data(), basin_size.data(), m_n_basins);

  result.resize(m_n_basins);
  GlobalSum(m_grid->com, thermal_forcing_sum.data(), result.data(), m_n_basins);

  // Divide by the size of a basin to compute the average:
  for (int basin_id = 0; basin_id < m_n_basins; basin_id++) {
    if (basin_size[basin_id] > 0) {
      result[basin_id] /= basin_size[basin_id];
    } else {
      result[basin_id] = 0.0;
    }
  }
}


//! \brief Computes mass flux in [kg m-2 s-1].
/*!
 * NON-LOCAL ISMIP6 Parameterisation (see Jourdain et al., 2019)
 */
void ISMIP6nl::mass_flux(const array::Scalar &thermal_forcing,
                         const array::Scalar &basin_mask,
                         std::vector<double> &basin_TF,
                         array::Scalar &result) {
  const double
    L                 = m_config->get_number("constants.fresh_water.latent_heat_of_fusion"),
    sea_water_density = m_config->get_number("constants.sea_water.density"),
    ice_density       = m_config->get_number("constants.ice.density"),
    c_p_ocean         = 3974.0, // J/(K*kg), specific heat capacity of ocean mixed layer
    gamma_0           = 3.52e-4;   // m/s, 11100 m/yr from local_MeanAnt method (Jourdain et al. 2020)

  array::AccessScope list{&basin_mask, &thermal_forcing, &result};

  // Calculate thermal forcing
  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    auto basin_id = basin_mask.as_int(i, j);

    double
      ocean_heat_flux = pow((sea_water_density * c_p_ocean) / (L * ice_density),2)  * gamma_0 * thermal_forcing(i,j)*basin_TF[basin_id]; // in W/m^2

    result(i,j) = ocean_heat_flux; // m s-1

    // convert from [m s-1] to [kg m-2 s-1]:
    result(i,j) *= ice_density;
  }
}

} // end of namespace ocean
} // end of namespace pism
