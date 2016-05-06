/* Copyright (C) 2016 PISM Authors
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
#include "CalvingFrontRetreat.hh"

#include "base/util/iceModelVec.hh"
#include "base/util/IceModelVec2CellType.hh"
#include "remove_narrow_tongues.hh"
#include "base/util/MaxTimestep.hh"
#include "base/util/PISMVars.hh"
#include "base/util/pism_utilities.hh"

namespace pism {

CalvingFrontRetreat::CalvingFrontRetreat(IceGrid::ConstPtr g)
  : Component(g) {

  m_thk_loss.create(m_grid, "temporary_storage", WITH_GHOSTS, 1);

  m_horizontal_calving_rate.create(m_grid, "horizontal_calving_rate", WITHOUT_GHOSTS);
  m_horizontal_calving_rate.set_attrs("diagnostic", "calving rate", "m second-1", "");
  m_horizontal_calving_rate.set_time_independent(false);
  m_horizontal_calving_rate.metadata().set_string("glaciological_units", "m year-1");
  m_horizontal_calving_rate.write_in_glaciological_units = true;

  m_restrict_timestep = m_config->get_boolean("cfl_eigen_calving");
}

CalvingFrontRetreat::~CalvingFrontRetreat() {

}

/**
 * @brief Compute the maximum time-step length allowed by the CFL
 * condition applied to the calving rate.
 *
 * Note: this code uses the mask variable obtained from the Vars
 * dictionary. This is not the same mask that is used in the update()
 * call, since max_timestep() is called *before* the mass-continuity
 * time-step.
 *
 * @return 0 on success
 */
MaxTimestep CalvingFrontRetreat::max_timestep() {

  if (not m_restrict_timestep) {
    return MaxTimestep();
  }

  // About 9 hours which corresponds to 10000 km year-1 on a 10 km grid
  double dt_min = units::convert(m_sys, 0.001, "years", "seconds");

  double
    calving_rate_max  = 0.0,
    calving_rate_mean = 0.0;
  int N_calving_cells = 0;

  const IceModelVec2CellType &mask = *m_grid->variables().get_2d_cell_type("mask");

  compute_calving_rate(mask, m_horizontal_calving_rate);

  IceModelVec::AccessList list;
  list.add(m_horizontal_calving_rate);

  for (Points pt(*m_grid); pt; pt.next()) {
    const int i = pt.i(), j = pt.j();

    const double C = m_horizontal_calving_rate(i, j);

    if (C > 0.0) {
      N_calving_cells   += 1;
      calving_rate_mean += C;
      calving_rate_max   = std::max(C, calving_rate_max);
    }
  }

  N_calving_cells   = GlobalSum(m_grid->com, N_calving_cells);
  calving_rate_mean = GlobalSum(m_grid->com, calving_rate_mean);
  calving_rate_max  = GlobalMax(m_grid->com, calving_rate_max);

  if (N_calving_cells > 0.0) {
    calving_rate_mean /= N_calving_cells;
  } else {
    calving_rate_mean = 0.0;
  }

  using units::convert;

  double denom = calving_rate_max / m_grid->dx();
  const double epsilon = convert(m_sys, 0.001 / (m_grid->dx() + m_grid->dy()), "seconds", "years");

  double dt = 1.0 / (denom + epsilon);

  m_log->message(3,
                 "  calving: maximum rate = %.2f m/year gives dt=%.5f years\n"
                 "           mean rate    = %.2f m/year over %d cells\n",
                 convert(m_sys, calving_rate_max, "m second-1", "m year-1"),
                 convert(m_sys, dt, "seconds", "years"),
                 convert(m_sys, calving_rate_mean, "m second-1", "m year-1"),
                 N_calving_cells);

  return MaxTimestep(std::max(dt, dt_min));
}

//! Update ice thickness, ice volume in partially-filled cells (Href),
//! and cell type mask by applying the 2D horizontal calving rate.

/*
  FIXME: we don't really need to call remove_narrow_tongues here: it is necessary when we use a
  calving parameterization which uses strain rates (eigen-calving), but it may not be appropriate
  with a frontal melt parameterization.
 */
void CalvingFrontRetreat::update(double dt,
                                 double sea_level,
                                 const IceModelVec2S &bed_topography,
                                 IceModelVec2CellType &mask,
                                 IceModelVec2S &Href,
                                 IceModelVec2S &ice_thickness) {

  compute_calving_rate(mask, m_horizontal_calving_rate);

  const double dx = m_grid->dx();

  m_thk_loss.set(0.0);

  IceModelVec::AccessList list;
  list.add(ice_thickness);
  list.add(mask);
  list.add(Href);
  list.add(m_thk_loss);
  list.add(m_horizontal_calving_rate);

  // Prepare to loop over neighbors: directions
  const Direction dirs[] = {North, East, South, West};

  // Apply the computed horizontal calving rate:
  for (Points pt(*m_grid); pt; pt.next()) {
    const int i = pt.i(), j = pt.j();

    // Find partially filled or empty grid boxes on the icefree ocean, which
    // have floating ice neighbors after the mass continuity step
    if (mask.ice_free_ocean(i, j) and mask.next_to_floating_ice(i, j)) {

      // Compute the number of floating neighbors and the neighbor-averaged ice thickness:
      double H_average = 0.0;
      int N_floating_neighbors = 0;
      {
        StarStencil<int> M_star = mask.int_star(i, j);
        StarStencil<double> H_star = ice_thickness.star(i, j);

        for (int n = 0; n < 4; ++n) {
          Direction direction = dirs[n];
          if (mask::floating_ice(M_star[direction])) {
            H_average += H_star[direction];
            N_floating_neighbors += 1;
          }
        }

        if (N_floating_neighbors > 0) {
          H_average /= N_floating_neighbors;
        }
      }

      // Calculate mass loss with respect to the associated ice thickness and the grid size:
      double calving_rate = m_horizontal_calving_rate(i, j) * H_average / dx; // in m/s

      // Apply calving rate at partially filled or empty grid cells
      if (calving_rate > 0.0) {
        Href(i, j) -= calving_rate * dt; // in m

        if (Href(i, j) < 0.0) {
          // Partially filled grid cell became ice-free

          m_thk_loss(i, j) = -Href(i, j); // in m, corresponds to additional ice loss
          Href(i, j)       = 0.0;

          // additional mass loss will be distributed among
          // N_floating_neighbors:
          if (N_floating_neighbors > 0) {
            m_thk_loss(i, j) /= N_floating_neighbors;
          }
        }
      }

    } // end of "if (ice_free_ocean and next_to_floating)"
  }

  m_thk_loss.update_ghosts();

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();
    double thk_loss_ij = 0.0;

    if (mask.floating_ice(i, j) and
        (m_thk_loss(i + 1, j) > 0.0 or m_thk_loss(i - 1, j) > 0.0 or
         m_thk_loss(i, j + 1) > 0.0 or m_thk_loss(i, j - 1) > 0.0)) {

      thk_loss_ij = (m_thk_loss(i + 1, j) + m_thk_loss(i - 1, j) +
                     m_thk_loss(i, j + 1) + m_thk_loss(i, j - 1));     // in m/s

      // Note std::max: we do not account for further calving
      // ice-inwards! Alternatively CFL criterion for time stepping
      // could be adjusted to maximum of calving rate
      Href(i, j) = std::max(ice_thickness(i, j) - thk_loss_ij, 0.0); // in m

      ice_thickness(i, j) = 0.0;
    }
  }

  // need to update ghosts of thickness to compute mask in place
  ice_thickness.update_ghosts();

  // update mask
  GeometryCalculator gc(*m_config);
  gc.set_icefree_thickness(m_config->get_double("mask_icefree_thickness_stress_balance_standard"));
  gc.compute_mask(sea_level, bed_topography, ice_thickness, mask);

  // remove narrow ice tongues
  remove_narrow_tongues(mask, ice_thickness);

  // update mask again
  gc.compute_mask(sea_level, bed_topography, ice_thickness, mask);
}

} // end of namespace pism
