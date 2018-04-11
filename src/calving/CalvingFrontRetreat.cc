/* Copyright (C) 2016, 2017, 2018 PISM Authors
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
#include "remove_narrow_tongues.hh"

#include "pism/util/iceModelVec.hh"
#include "pism/util/IceModelVec2CellType.hh"
#include "pism/util/MaxTimestep.hh"
#include "pism/util/pism_utilities.hh"
#include "pism/geometry/part_grid_threshold_thickness.hh"
#include "pism/geometry/Geometry.hh"

namespace pism {

CalvingInputs::CalvingInputs() {
  geometry = nullptr;

  bc_mask              = nullptr;
  ice_enthalpy         = nullptr;
  ice_velocity         = nullptr;
  shelf_base_mass_flux = nullptr;
}

CalvingFrontRetreat::CalvingFrontRetreat(IceGrid::ConstPtr g, unsigned int mask_stencil_width)
  : Component(g) {

  m_tmp.create(m_grid, "temporary_storage", WITH_GHOSTS, 1);
  m_tmp.set_attrs("internal", "additional mass loss at points near the calving front",
                  "m", "");

  m_horizontal_calving_rate.create(m_grid, "horizontal_calving_rate", WITHOUT_GHOSTS);
  m_horizontal_calving_rate.set_attrs("diagnostic", "calving rate", "m second-1", "land_ice_calving_rate");
  m_horizontal_calving_rate.set_time_independent(false);
  m_horizontal_calving_rate.metadata().set_string("glaciological_units", "m year-1");

  m_mask.create(m_grid, "m_mask", WITH_GHOSTS, mask_stencil_width);
  m_mask.set_attrs("internal", "cell type mask", "", "");

  m_surface_topography.create(m_grid, "m_surface_topography", WITH_GHOSTS, 1);
  m_surface_topography.set_attrs("internal", "surface topography", "m", "surface_altitude");

  m_restrict_timestep = m_config->get_boolean("calving.front_retreat.use_cfl");
}

CalvingFrontRetreat::~CalvingFrontRetreat() {
  // empty
}

/**
 * @brief Compute the maximum time-step length allowed by the CFL
 * condition applied to the calving rate.
 */
MaxTimestep CalvingFrontRetreat::max_timestep(const CalvingInputs &inputs,
                                              double t) const {
  (void) t;

  if (not m_restrict_timestep) {
    return MaxTimestep();
  }

  // About 9 hours which corresponds to 10000 km year-1 on a 10 km grid
  double dt_min = units::convert(m_sys, 0.001, "years", "seconds");

  double
    calving_rate_max  = 0.0,
    calving_rate_mean = 0.0;
  int N_calving_cells = 0;

  IceModelVec2S &horizontal_calving_rate = m_tmp;

  compute_calving_rate(inputs, horizontal_calving_rate);

  IceModelVec::AccessList list(horizontal_calving_rate);

  for (Points pt(*m_grid); pt; pt.next()) {
    const int i = pt.i(), j = pt.j();

    const double C = horizontal_calving_rate(i, j);

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

/*!
 * Adjust the mask near domain boundaries to avoid "wrapping around."
 */
void CalvingFrontRetreat::prepare_mask(const IceModelVec2CellType &input,
                                       IceModelVec2CellType &output) const {

  output.copy_from(input);

  if (m_config->get_boolean("calving.front_retreat.wrap_around")) {
    return;
  }

  IceModelVec::AccessList list(output);

  const int Mx = m_grid->Mx();
  const int My = m_grid->My();

  ParallelSection loop(m_grid->com);
  try {
    for (PointsWithGhosts p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      if (i < 0 or i >= Mx or j < 0 or j >= My) {
        output(i, j) = MASK_ICE_FREE_OCEAN;
      }
    }
  } catch (...) {
    loop.failed();
  }
  loop.check();
}

/*! Update ice geometry and mask using the computed horizontal calving rate.
 * @param[in] dt time step, seconds
 * @param[in] sea_level sea level elevation, meters
 * @param[in] thickness_bc_mask Dirichlet B.C. mask for the ice thickness
 * @param[in] bed_topography bed elevation, meters
 * @param[in,out] mask cell type mask
 * @param[in,out] Href "area specific volume"
 * @param[in,out] ice_thickness ice thickness
 *
 * FIXME: we don't really need to call remove_narrow_tongues here: it is necessary when we use a
 * calving parameterization which uses strain rates (eigen-calving), but it may not be appropriate
 * with a frontal melt parameterization.
 */
void CalvingFrontRetreat::update(double dt,
                                 const CalvingInputs &inputs,
                                 IceModelVec2CellType &mask,
                                 IceModelVec2S &Href,
                                 IceModelVec2S &ice_thickness) {

  const IceModelVec2S   &sea_level      = inputs.geometry->sea_level_elevation;
  const IceModelVec2S   &bed_topography = inputs.geometry->bed_elevation;
  const IceModelVec2Int &bc_mask        = *inputs.bc_mask;

  GeometryCalculator gc(*m_config);
  gc.compute_surface(sea_level, bed_topography, ice_thickness, m_surface_topography);

  // use mask with a wide stencil to compute the calving rate
  compute_calving_rate(inputs, m_horizontal_calving_rate);

  const double dx = m_grid->dx();

  m_tmp.set(0.0);

  IceModelVec::AccessList list{&ice_thickness, &bc_mask,
      &bed_topography, &sea_level, &mask, &Href, &m_tmp, &m_horizontal_calving_rate,
      &m_surface_topography};

  // Prepare to loop over neighbors: directions
  const Direction dirs[] = {North, East, South, West};

  // Step 1: Apply the computed horizontal calving rate:
  for (Points pt(*m_grid); pt; pt.next()) {
    const int i = pt.i(), j = pt.j();

    if (bc_mask(i, j) > 0.5) {
      // don't modify cells marked as Dirichlet B.C. locations
      continue;
    }

    const double rate = m_horizontal_calving_rate(i, j);

    if (mask.ice_free(i, j) and rate > 0.0) {
      // apply calving rate at the margin (i.e. to partially-filled cells) only

      const double Href_old = Href(i, j);

      // Compute the number of floating neighbors and the neighbor-averaged ice thickness:
      double H_threshold = part_grid_threshold_thickness(mask.int_star(i, j),
                                                         ice_thickness.star(i, j),
                                                         m_surface_topography.star(i, j),
                                                         bed_topography(i, j));

      // Calculate mass loss with respect to the associated ice thickness and the grid size:
      const double Href_change = -dt * rate * H_threshold / dx; // in m

      if (Href_old + Href_change >= 0.0) {
        // Href is high enough to absorb the mass loss
        Href(i, j) = Href_old + Href_change;
      } else {
        Href(i, j) = 0.0;
        // Href is below Href_change: need to distribute mass loss to neighboring points

        // Find the number of neighbors to distribute to.
        //
        // We consider floating cells and grounded cells with the base below sea level. In other
        // words, additional mass losses are distributed to shelf calving fronts and grounded marine
        // termini.
        int N = 0;
        {
          auto
            M  = mask.int_star(i, j),
            BC = bc_mask.int_star(i, j);

          auto
            bed = bed_topography.star(i, j),
            sl  = sea_level.star(i, j);

          for (int n = 0; n < 4; ++n) {
            Direction direction = dirs[n];
            int m = M[direction];
            int bc = BC[direction];

            if (bc == 0 and     // distribute to regular (*not* Dirichlet B.C.) neighbors only
                (mask::floating_ice(m) or
                 (mask::grounded_ice(m) and bed[direction] < sl[direction]))) {
              N += 1;
            }
          }
        }

        if (N > 0) {
          m_tmp(i, j) = (Href_old + Href_change) / (double)N;
        } else {
          // No shelf calving front of grounded terminus to distribute to: calving stops here.
          m_tmp(i, j) = 0.0;
        }
      }

    } // end of "if (rate > 0.0)"
  }   // end of loop over grid points

  // Step 2: update ice thickness and Href in neighboring cells if we need to propagate mass losses
  // due to calving front retreat.
  m_tmp.update_ghosts();

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    // Note: this condition has to match the one in step 1 above.
    if (bc_mask.as_int(i, j) == 0 and
        (mask.floating_ice(i, j) or
         (mask.grounded_ice(i, j) and bed_topography(i, j) < sea_level(i, j)))) {

      const double delta_H = (m_tmp(i + 1, j) + m_tmp(i - 1, j) +
                              m_tmp(i, j + 1) + m_tmp(i, j - 1));

      if (delta_H < 0.0) {
        Href(i, j) = ice_thickness(i, j) + delta_H; // in m
        ice_thickness(i, j) = 0.0;
      }

      // Stop calving if the current cell does not have enough ice to absorb the loss.
      if (Href(i, j) < 0.0) {
        Href(i, j) = 0.0;
      }

    }
  }

  // need to update ghosts of thickness to compute mask in place
  ice_thickness.update_ghosts();

  // update mask
  gc.set_icefree_thickness(m_config->get_double("stress_balance.ice_free_thickness_standard"));
  gc.compute_mask(sea_level, bed_topography, ice_thickness, mask);

  // remove narrow ice tongues
  remove_narrow_tongues(mask, ice_thickness);

  // update mask again
  gc.compute_mask(sea_level, bed_topography, ice_thickness, mask);
}

const IceModelVec2S& CalvingFrontRetreat::calving_rate() const {
  return m_horizontal_calving_rate;
}

CalvingRate::CalvingRate(const CalvingFrontRetreat *m,
                         const std::string &name,
                         const std::string &long_name)
  : Diag<CalvingFrontRetreat>(m) {

  /* set metadata: */
  m_vars = {SpatialVariableMetadata(m_sys, name)};

  set_attrs(long_name, "",      // land_ice_calving_rate
            "m second-1", "m year-1", 0);
}

IceModelVec::Ptr CalvingRate::compute_impl() const {

  IceModelVec2S::Ptr result(new IceModelVec2S(m_grid, "", WITHOUT_GHOSTS));
  result->metadata(0) = m_vars[0];

  result->copy_from(model->calving_rate());

  return result;
}

} // end of namespace pism
