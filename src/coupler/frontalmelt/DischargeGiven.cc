// Copyright (C) 2018, 2019, 2021 Andy Aschwanden and Constantine Khroulev
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

#include "DischargeGiven.hh"

#include "pism/util/IceGrid.hh"
#include "pism/geometry/Geometry.hh"
#include "pism/coupler/util/options.hh"
#include "FrontalMeltPhysics.hh"

namespace pism {
namespace frontalmelt {
  
DischargeGiven::DischargeGiven(IceGrid::ConstPtr grid)
  : FrontalMelt(grid, nullptr) {

  m_frontal_melt_rate = allocate_frontal_melt_rate(grid, 1);

  m_log->message(2,
                 "* Initializing the frontal melt model\n"
                 "  UAF-UT\n");

  m_theta_ocean = IceModelVec2T::Constant(grid, "theta_ocean", 0.0);

  m_subglacial_discharge = IceModelVec2T::Constant(grid, "subglacial_discharge", 0.0);
}

DischargeGiven::~DischargeGiven() {
  // empty
}

void DischargeGiven::init_impl(const Geometry &geometry) {
  (void) geometry;

  ForcingOptions opt(*m_grid->ctx(), "frontal_melt.discharge_given");

  {
    unsigned int buffer_size = m_config->get_number("input.forcing.buffer_size");
    unsigned int evaluations_per_year = m_config->get_number("input.forcing.evaluations_per_year");

    File file(m_grid->com, opt.filename, PISM_NETCDF3, PISM_READONLY);

    m_theta_ocean = IceModelVec2T::ForcingField(m_grid,
                                                file,
                                                "theta_ocean",
                                                "", // no standard name
                                                buffer_size,
                                                evaluations_per_year,
                                                opt.periodic);

    m_subglacial_discharge = IceModelVec2T::ForcingField(m_grid,
                                                file,
                                                "subglacial_discharge",
                                                "", // no standard name
                                                buffer_size,
                                                evaluations_per_year,
                                                opt.periodic);
  }

  m_theta_ocean->set_attrs("climate_forcing",
                           "potential temperature of the adjacent ocean",
                           "Celsius", "Celsius", "", 0);

  m_theta_ocean->init(opt.filename, opt.periodic);

  m_subglacial_discharge->set_attrs("climate_forcing",
                                    "subglacial discharge",
                                    "kg m-2 s-1", "kg m-2 year-1", "", 0);

  m_subglacial_discharge->init(opt.filename, opt.periodic);
}

/*!
 * Initialize potential temperature from IceModelVecs instead of an input
 * file (for testing).
 */
void DischargeGiven::initialize(const IceModelVec2S &theta, const IceModelVec2S &sgl) {
  m_theta_ocean->copy_from(theta);
  m_subglacial_discharge->copy_from(sgl);
}

void DischargeGiven::update_impl(const FrontalMeltInputs &inputs, double t, double dt) {

  m_theta_ocean->update(t, dt);
  m_subglacial_discharge->update(t, dt);

  m_theta_ocean->average(t, dt);
  m_subglacial_discharge->average(t, dt);

  FrontalMeltPhysics physics(*m_config);

  const IceModelVec2CellType &cell_type           = inputs.geometry->cell_type;
  const IceModelVec2S        &bed_elevation       = inputs.geometry->bed_elevation;
  const IceModelVec2S        &sea_level_elevation = inputs.geometry->sea_level_elevation;

  IceModelVec::AccessList list
    {&bed_elevation, &cell_type, &sea_level_elevation,
     m_theta_ocean.get(), m_subglacial_discharge.get(), m_frontal_melt_rate.get()};

  double
    water_density   = m_config->get_number("constants.fresh_water.density"),
    seconds_per_day = 86400;

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    if (cell_type.icy(i, j)) {
      // Assume for now that thermal forcing is equal to theta_ocean. Also, thermal
      // forcing is generally not available at the grounding line.
      double TF = (*m_theta_ocean)(i, j);

      // Convert subglacial discharge (kg /(m^2 * s)) to an "effective subglacial freshwater
      // velocity" or volume flux per unit area of ice front in m/day (see Xu et al 2013,
      // section 2, paragraph 11).
      //
      // [flux] = kg / (m^2 * s), so
      // [flux / water_density] = m / s, so
      // [flux / water_density  * (s / day) ] = m / day
      double q_sg = ((*m_subglacial_discharge)(i, j) / water_density) * seconds_per_day;

      double water_depth = sea_level_elevation(i, j) - bed_elevation(i, j);

      (*m_frontal_melt_rate)(i, j) = physics.frontal_melt_from_undercutting(water_depth, q_sg, TF);
      // convert from m / day to m / s
      (*m_frontal_melt_rate)(i, j) /= seconds_per_day;
    } else {
      (*m_frontal_melt_rate)(i, j) = 0.0;
    }
  } // end of the loop over grid points

  // Set frontal melt rate *near* grounded termini to the average of grounded icy
  // neighbors: front retreat code uses values at these locations (the rest is for
  // visualization).

  m_frontal_melt_rate->update_ghosts();

  const Direction dirs[] = {North, East, South, West};

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    if (apply(cell_type, i, j) and cell_type.ice_free(i, j)) {

      auto R = m_frontal_melt_rate->star(i, j);
      auto M = cell_type.int_star(i, j);

      int N = 0;
      double R_sum = 0.0;
      for (int n = 0; n < 4; ++n) {
        Direction direction = dirs[n];
        if (mask::grounded_ice(M[direction]) or
            (m_include_floating_ice and mask::icy(M[direction]))) {
          R_sum += R[direction];
          N++;
        }
      }

      if (N > 0) {
        (*m_frontal_melt_rate)(i, j) = R_sum / N;
      }
    }
  }
}

const IceModelVec2S& DischargeGiven::frontal_melt_rate_impl() const {
  return *m_frontal_melt_rate;
}

MaxTimestep DischargeGiven::max_timestep_impl(double t) const {

  auto dt = std::min(m_theta_ocean->max_timestep(t),
                     m_subglacial_discharge->max_timestep(t));

  if (dt.finite()) {
    return MaxTimestep(dt.value(), "frontal_melt discharge_given");
  } else {
    return MaxTimestep("frontal_melt discharge_given");
  }
}

} // end of namespace frontalmelt
} // end of namespace pism
