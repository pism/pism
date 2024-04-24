// Copyright (C) 2018, 2019, 2021, 2022, 2023, 2024 Andy Aschwanden and Constantine Khroulev
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

#include "pism/coupler/frontalmelt/DischargeRouting.hh"

#include "pism/util/Grid.hh"
#include "pism/geometry/Geometry.hh"
#include "pism/coupler/util/options.hh"
#include "pism/coupler/frontalmelt/FrontalMeltPhysics.hh"
#include "pism/util/array/Forcing.hh"

namespace pism {
namespace frontalmelt {
  
DischargeRouting::DischargeRouting(std::shared_ptr<const Grid> grid)
  : FrontalMelt(grid, nullptr),
    m_frontal_melt_rate(grid, "frontal_melt_rate") {

  m_frontal_melt_rate.metadata(0)
      .long_name("frontal melt rate")
      .units("m s^-1")
      .output_units("m day^-1");

  m_log->message(2, "* Initializing the frontal melt model\n"
                    "  using the Rignot/Xu parameterization\n"
                    "  and routing of subglacial discharge\n");

  m_theta_ocean = array::Forcing::Constant(grid, "theta_ocean", 0.0);
}

void DischargeRouting::init_impl(const Geometry &geometry) {
  (void)geometry;

  ForcingOptions opt(*m_grid->ctx(), "frontal_melt.routing");

  {
    unsigned int buffer_size = m_config->get_number("input.forcing.buffer_size");

    File file(m_grid->com, opt.filename, io::PISM_NETCDF3, io::PISM_READONLY);

    m_theta_ocean = std::make_shared<array::Forcing>(m_grid, file, "theta_ocean",
                                                     "", // no standard name
                                                     buffer_size, opt.periodic, LINEAR);
  }

  m_theta_ocean->metadata(0)
      .long_name("potential temperature of the adjacent ocean")
      .units("degree_Celsius");

  m_theta_ocean->init(opt.filename, opt.periodic);
}

/*!
 * Initialize potential temperature from an array instead of an input
 * file (for testing).
 */
void DischargeRouting::initialize(const array::Scalar &theta) {
  m_theta_ocean->copy_from(theta);
}

void DischargeRouting::update_impl(const FrontalMeltInputs &inputs, double t, double dt) {

  m_theta_ocean->update(t, dt);

  FrontalMeltPhysics physics(*m_config);

  const auto &cell_type                    = inputs.geometry->cell_type;
  const array::Scalar &bed_elevation       = inputs.geometry->bed_elevation;
  const array::Scalar &ice_thickness       = inputs.geometry->ice_thickness;
  const array::Scalar &sea_level_elevation = inputs.geometry->sea_level_elevation;
  const array::Scalar &water_flux          = *inputs.subglacial_water_flux;

  array::AccessScope list{ &ice_thickness,       &bed_elevation, &cell_type,
                           &sea_level_elevation, &water_flux,    m_theta_ocean.get(),
                           &m_frontal_melt_rate };

  double seconds_per_day = 86400, grid_spacing = 0.5 * (m_grid->dx() + m_grid->dy());

  for (auto p = m_grid->points(); p; p.next()) {
    const int i = p.i(), j = p.j();

    if (cell_type.icy(i, j)) {
      // Assume for now that thermal forcing is equal to theta_ocean. Also, thermal
      // forcing is generally not available at the grounding line.
      double TF = (*m_theta_ocean)(i, j);

      double water_depth          = std::max(sea_level_elevation(i, j) - bed_elevation(i, j), 0.0),
             submerged_front_area = water_depth * grid_spacing;

      // Convert subglacial water flux (m^2/s) to an "effective subglacial freshwater
      // velocity" or flux per unit area of ice front in m/day (see Xu et al 2013, section
      // 2, paragraph 11).
      //
      // [flux] = m^2 / s, so
      // [flux * grid_spacing] = m^3 / s, so
      // [flux * grid_spacing / submerged_front_area] = m / s, and
      // [flux * grid_spacing  * (s / day) / submerged_front_area] = m / day
      double Q_sg = water_flux(i, j) * grid_spacing;
      double q_sg = Q_sg / submerged_front_area * seconds_per_day;

      m_frontal_melt_rate(i, j) = physics.frontal_melt_from_undercutting(water_depth, q_sg, TF);
      // convert from m / day to m / s
      m_frontal_melt_rate(i, j) /= seconds_per_day;
    } else {
      m_frontal_melt_rate(i, j) = 0.0;
    }
  } // end of the loop over grid points

  // Set frontal melt rate *near* grounded termini to the average of grounded icy
  // neighbors: front retreat code uses values at these locations (the rest is for
  // visualization).

  m_frontal_melt_rate.update_ghosts();

  for (auto p = m_grid->points(); p; p.next()) {
    const int i = p.i(), j = p.j();

    if (apply(cell_type, i, j) and cell_type.ice_free(i, j)) {

      auto R = m_frontal_melt_rate.star(i, j);
      auto M = cell_type.star_int(i, j);

      int N        = 0;
      double R_sum = 0.0;
      for (auto d : { North, East, South, West }) {
        if (mask::grounded_ice(M[d]) or (m_include_floating_ice and mask::icy(M[d]))) {
          R_sum += R[d];
          N++;
        }
      }

      if (N > 0) {
        m_frontal_melt_rate(i, j) = R_sum / N;
      }
    }
  }
}

const array::Scalar &DischargeRouting::frontal_melt_rate_impl() const {
  return m_frontal_melt_rate;
}

MaxTimestep DischargeRouting::max_timestep_impl(double t) const {

  auto dt = m_theta_ocean->max_timestep(t);

  if (dt.finite()) {
    return {dt.value(), "frontal_melt routing"};
  }

  return {"frontal_melt routing"};
}

} // end of namespace frontalmelt
} // end of namespace pism
