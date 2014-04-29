/* Copyright (C) 2013, 2014 PISM Authors
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

#include "PISMEigenCalving.hh"
#include "PISMVars.hh"
#include "PISMStressBalance.hh"
#include "Mask.hh"

namespace pism {

EigenCalving::EigenCalving(IceGrid &g, const Config &conf,
                                   StressBalance *stress_balance)
  : Component(g, conf), m_stencil_width(2), m_mask(NULL),
    m_stress_balance(stress_balance) {
  PetscErrorCode ierr;
  ierr = m_strain_rates.create(grid, "edot", WITH_GHOSTS,
                               m_stencil_width,
                               2); CHKERRCONTINUE(ierr);
  if (ierr != 0) {
    PetscPrintf(grid.com, "PISM ERROR: memory allocation failed.\n");
    PISMEnd();
  }

  ierr = m_thk_loss.create(grid, "temporary_storage", WITH_GHOSTS, 1); CHKERRCONTINUE(ierr);
  if (ierr != 0) {
    PetscPrintf(grid.com, "PISM ERROR: memory allocation failed.\n");
    PISMEnd();
  }

  m_strain_rates.set_name("edot_1", 0);
  m_strain_rates.set_attrs("internal",
                           "major principal component of horizontal strain-rate",
                           "1/s", "", 0);

  m_strain_rates.set_name("edot_2", 1);
  m_strain_rates.set_attrs("internal",
                           "minor principal component of horizontal strain-rate",
                           "1/s", "", 1);

  m_K = config.get("eigen_calving_K");
  m_restrict_timestep = config.get_flag("cfl_eigen_calving");
}

EigenCalving::~EigenCalving() {
  // empty
}

PetscErrorCode EigenCalving::init(Vars &vars) {
  PetscErrorCode ierr;

  ierr = verbPrintf(2, grid.com,
                    "* Initializing the 'eigen-calving' mechanism...\n");
  CHKERRQ(ierr);

  if (PetscAbs(grid.dx - grid.dy) / PetscMin(grid.dx, grid.dy) > 1e-2) {
    PetscPrintf(grid.com,
                "PISM ERROR: -calving eigen_calving using a non-square grid cell is not implemented (yet);\n"
                "            dx = %f, dy = %f, relative difference = %f",
                grid.dx, grid.dy,
                PetscAbs(grid.dx - grid.dy) / PetscMax(grid.dx, grid.dy));
    PISMEnd();
  }

  ierr = m_strain_rates.set(0.0); CHKERRQ(ierr);

  m_mask = dynamic_cast<IceModelVec2Int*>(vars.get("mask"));
  if (m_mask == NULL) SETERRQ(grid.com, 1, "mask is not available");

  return 0;
}

//! \brief Uses principal strain rates to apply "eigencalving" with constant K.
/*!
  See equation (26) in [\ref Winkelmannetal2011].
*/
PetscErrorCode EigenCalving::update(double dt,
                                        IceModelVec2Int &pism_mask,
                                        IceModelVec2S &Href,
                                        IceModelVec2S &ice_thickness) {
  PetscErrorCode ierr;

  // Distance (grid cells) from calving front where strain rate is evaluated
  int offset = m_stencil_width;

  ierr = m_thk_loss.set(0.0); CHKERRQ(ierr);

  ierr = update_strain_rates(); CHKERRQ(ierr);

  MaskQuery mask(pism_mask);

  ierr = ice_thickness.begin_access(); CHKERRQ(ierr);
  ierr = pism_mask.begin_access(); CHKERRQ(ierr);
  ierr = Href.begin_access(); CHKERRQ(ierr);
  ierr = m_strain_rates.begin_access(); CHKERRQ(ierr);
  ierr = m_thk_loss.begin_access(); CHKERRQ(ierr);
  for (int i = grid.xs; i < grid.xs + grid.xm; ++i) {
    for (int j = grid.ys; j < grid.ys + grid.ym; ++j) {

      // Average of strain-rate eigenvalues in adjacent floating grid
      // cells to be used for eigen-calving:
      double
        eigen1    = 0.0,
        eigen2    = 0.0;
      // Neighbor-averaged ice thickness:
      double
        H_average = 0.0;

      int N_floating_neighbors = 0, M = 0;
      // M is the number of cells used to compute averages of strain
      // rate components.

      // Find partially filled or empty grid boxes on the icefree ocean, which
      // have floating ice neighbors after the mass continuity step
      if (mask.ice_free_ocean(i, j) &&
          mask.next_to_floating_ice(i, j)) {

        if (mask.floating_ice(i + 1, j)) {
          N_floating_neighbors += 1;
          H_average += ice_thickness(i + 1, j);
        }

        if (mask.floating_ice(i - 1, j)) {
          N_floating_neighbors += 1;
          H_average += ice_thickness(i - 1, j);
        }

        if (mask.floating_ice(i, j + 1)) {
          N_floating_neighbors += 1;
          H_average += ice_thickness(i, j + 1);
        }

        if (mask.floating_ice(i, j - 1)) {
          N_floating_neighbors += 1;
          H_average += ice_thickness(i, j - 1);
        }

        if (N_floating_neighbors > 0)
          H_average /= N_floating_neighbors;

        for (int p = -1; p < 2; p += 2) {
          int i_offset = p * offset;
          if (mask.floating_ice(i + i_offset, j) &&
              mask.ice_margin(i + i_offset, j) == false) {
            eigen1 += m_strain_rates(i + i_offset, j, 0);
            eigen2 += m_strain_rates(i + i_offset, j, 1);
            M += 1;
          }
        }

        for (int q = -1; q < 2; q += 2) {
          int j_offset = q * offset;
          if (mask.floating_ice(i, j + j_offset) &&
              mask.ice_margin(i , j + j_offset) == false) {
            eigen1 += m_strain_rates(i, j + j_offset, 0);
            eigen2 += m_strain_rates(i, j + j_offset, 1);
            M += 1;
          }
        }

        if (M > 0) {
          eigen1 /= M;
          eigen2 /= M;
        }

        double
          calving_rate_horizontal = 0.0,
          eigenCalvOffset         = 0.0;
        // eigenCalvOffset allows adjusting the transition from
        // compressive to extensive flow regime

        // Calving law
        //
        // eigen1 * eigen2 has units [s^-2] and calving_rate_horizontal
        // [m*s^1] hence, eigen_calving_K has units [m*s]
        if (eigen2 > eigenCalvOffset &&
            eigen1 > 0.0) { // if spreading in all directions
          calving_rate_horizontal = m_K * eigen1 * (eigen2 - eigenCalvOffset);
        }

        // calculate mass loss with respect to the associated ice thickness and the grid size:
        double calving_rate = calving_rate_horizontal * H_average / grid.dx; // in m/s

        // apply calving rate at partially filled or empty grid cells
        if (calving_rate > 0.0) {
          Href(i, j) -= calving_rate * dt; // in m

          if (Href(i, j) < 0.0) {
            // Partially filled grid cell became ice-free

            m_thk_loss(i, j) = -Href(i, j); // in m, corresponds to additional ice loss
            Href(i, j)       = 0.0;

            // additional mass loss will be distributed among
            // N_floating_neighbors:
            if (N_floating_neighbors > 0)
              m_thk_loss(i, j) /= N_floating_neighbors;
          }
        }

      } // end of "if (ice_free_ocean && next_to_floating)"
    } // j-loop
  } // i-loop
  ierr = m_thk_loss.end_access(); CHKERRQ(ierr);

  ierr = m_thk_loss.update_ghosts(); CHKERRQ(ierr);

  ierr = m_thk_loss.begin_access(); CHKERRQ(ierr);
  for (int i = grid.xs; i < grid.xs + grid.xm; ++i) {
    for (int j = grid.ys; j < grid.ys + grid.ym; ++j) {
      double thk_loss_ij = 0.0;

      if (mask.floating_ice(i, j) &&
          (m_thk_loss(i + 1, j) > 0.0 || m_thk_loss(i - 1, j) > 0.0 ||
           m_thk_loss(i, j + 1) > 0.0 || m_thk_loss(i, j - 1) > 0.0)) {

        thk_loss_ij = (m_thk_loss(i + 1, j) + m_thk_loss(i - 1, j) +
                       m_thk_loss(i, j + 1) + m_thk_loss(i, j - 1));     // in m/s

        // Note PetscMax: we do not account for further calving
        // ice-inwards! Alternatively CFL criterion for time stepping
        // could be adjusted to maximum of calving rate
        Href(i, j) = PetscMax(ice_thickness(i, j) - thk_loss_ij, 0.0); // in m

        ice_thickness(i, j) = 0.0;
        pism_mask(i, j) = MASK_ICE_FREE_OCEAN;
      }

    } // j-loop
  } // i-loop
  ierr = ice_thickness.end_access(); CHKERRQ(ierr);
  ierr = pism_mask.end_access(); CHKERRQ(ierr);
  ierr = Href.end_access(); CHKERRQ(ierr);
  ierr = m_strain_rates.end_access(); CHKERRQ(ierr);
  ierr = m_thk_loss.end_access(); CHKERRQ(ierr);

  ierr = pism_mask.update_ghosts(); CHKERRQ(ierr);

  ierr = remove_narrow_tongues(pism_mask, ice_thickness); CHKERRQ(ierr);

  ierr = ice_thickness.update_ghosts(); CHKERRQ(ierr);
  ierr = pism_mask.update_ghosts(); CHKERRQ(ierr);

  return 0;
}


/**
 * @brief Compute the maximum time-step length allowed by the CFL
 * condition applied to the calving rate.
 *
 * @param[in] my_t current time, in seconds
 * @param[out] my_dt set the the maximum allowed time-step, in seconds
 * @param[out] restrict set to "true" if this component has a time-step restriction
 *
 * Note: this code uses the mask variable obtained from the Vars
 * dictionary. This is not the same mask that is used in the update()
 * call, since max_timestep() is called *before* the mass-continuity
 * time-step.
 *
 * @return 0 on success
 */
PetscErrorCode EigenCalving::max_timestep(double /*my_t*/,
                                              double &my_dt, bool &restrict) {
  PetscErrorCode ierr;

  if (m_restrict_timestep == false) {
    restrict = false;
    my_dt    = -1;
    return 0;
  }

  restrict = true;

  // About 9 hours which corresponds to 10000 km/year on a 10 km grid
  double dt_min = grid.convert(0.001, "years", "seconds");

  // Distance (grid cells) from calving front where strain rate is evaluated
  int offset = m_stencil_width,
    i0 = 0, j0 = 0;
  double
    my_calving_rate_max     = 0.0,
    my_calving_rate_mean    = 0.0,
    my_calving_rate_counter = 0.0;

  MaskQuery mask(*m_mask);

  ierr = update_strain_rates(); CHKERRQ(ierr);

  ierr = m_mask->begin_access(); CHKERRQ(ierr);
  ierr = m_strain_rates.begin_access(); CHKERRQ(ierr);

  for (int i = grid.xs; i < grid.xs + grid.xm; ++i) {
    for (int j = grid.ys; j < grid.ys + grid.ym; ++j) {
      // Average of strain-rate eigenvalues in adjacent floating grid cells to
      // be used for eigencalving
      double eigen1 = 0.0, eigen2 = 0.0;
      // Number of cells used in computing eigen1 and eigen2:
      int M = 0;

      // find partially filled or empty grid boxes on the ice-free
      // ocean which have floating ice neighbors
      if ((mask.ice_free_ocean(i, j) &&
           mask.next_to_floating_ice(i, j)) == false)
        continue;

      double
        calving_rate_horizontal = 0.0,
        eigenCalvOffset = 0.0;

      for (int p = -1; p < 2; p += 2) {
        int i_offset = p * offset;
        if (mask.floating_ice(i + i_offset, j) &&
            mask.ice_margin(i + i_offset, j) == false) {
          eigen1 += m_strain_rates(i + i_offset, j, 0);
          eigen2 += m_strain_rates(i + i_offset, j, 1);
          M += 1;
        }
      }

      for (int q = -1; q < 2; q += 2) {
        int j_offset = q * offset;
        if (mask.floating_ice(i, j + j_offset) &&
            mask.ice_margin(i,   j + j_offset) == false) {
          eigen1 += m_strain_rates(i, j + j_offset, 0);
          eigen2 += m_strain_rates(i, j + j_offset, 1);
          M += 1;
        }
      }

      if (M > 0) {
        eigen1 /= M;
        eigen2 /= M;
      }

      // calving law
      if (eigen2 > eigenCalvOffset && eigen1 > 0.0) { // if spreading in all directions
        calving_rate_horizontal = m_K * eigen1 * (eigen2 - eigenCalvOffset);
        my_calving_rate_counter += 1.0;
        my_calving_rate_mean += calving_rate_horizontal;
        if (my_calving_rate_max < calving_rate_horizontal) {
          i0 = i;
          j0 = j;
        }
        my_calving_rate_max = PetscMax(my_calving_rate_max, calving_rate_horizontal);
      } else calving_rate_horizontal = 0.0;

    } // i-loop
  } // j-loop

  ierr = m_mask->end_access(); CHKERRQ(ierr);
  ierr = m_strain_rates.end_access(); CHKERRQ(ierr);

  double calving_rate_max = 0.0, calving_rate_mean = 0.0, calving_rate_counter = 0.0;
  ierr = GlobalSum(&my_calving_rate_mean, &calving_rate_mean, grid.com); CHKERRQ(ierr);
  ierr = GlobalSum(&my_calving_rate_counter, &calving_rate_counter, grid.com); CHKERRQ(ierr);
  ierr = GlobalMax(&my_calving_rate_max, &calving_rate_max, grid.com); CHKERRQ(ierr);

  calving_rate_mean /= calving_rate_counter;

  double denom = calving_rate_max / grid.dx;
  const double epsilon = grid.convert(0.001 / (grid.dx + grid.dy), "seconds", "years");

  my_dt = 1.0 / (denom + epsilon);

  ierr = verbPrintf(2, grid.com,
                    "!!!!! c_rate = %.0f m/year (dt=%.5f a) at point %d, %d with mean_c=%.0f m/year over %.0f cells \n",
                    grid.convert(calving_rate_max, "m/s", "m/year"),
                    grid.convert(my_dt, "seconds", "years"),
                    i0, j0,
                    grid.convert(calving_rate_mean, "m/s", "m/year"),
                    calving_rate_counter); CHKERRQ(ierr);

  my_dt = PetscMax(my_dt, dt_min);

  return 0;
}

void EigenCalving::add_vars_to_output(std::string /*keyword*/, std::set<std::string> &/*result*/) {
  // empty
}

PetscErrorCode EigenCalving::define_variables(std::set<std::string> /*vars*/, const PIO &/*nc*/,
                                                  IO_Type /*nctype*/) {
  // empty
  return 0;
}

PetscErrorCode EigenCalving::write_variables(std::set<std::string> /*vars*/, const PIO& /*nc*/) {
  // empty
  return 0;
}

/**
 * Update the strain rates field.
 *
 * Note: this code uses the mask obtained from the Vars
 * dictionary, because the velocity field used to compute it need not
 * extend past the ice margin corresponding to the *beginning* of the
 * time-step.
 *
 * @return 0 on success
 */
PetscErrorCode EigenCalving::update_strain_rates() {
  PetscErrorCode ierr;

  IceModelVec2V *ssa_velocity;
  ierr = m_stress_balance->get_2D_advective_velocity(ssa_velocity); CHKERRQ(ierr);
  ierr = m_stress_balance->compute_2D_principal_strain_rates(*ssa_velocity,
                                                             *m_mask, m_strain_rates); CHKERRQ(ierr);

  return 0;
}

/** Remove tips of one-cell-wide ice tongues.
 *
 * ice tongues like this one (and equivalent)
 *
 * @code
   O O O
   X X O
   O O O
   @endcode
 *
 * are removed, while ones like this one
 *
 * @code
   X O O
   X X O
   X O O
   @endcode
 * are not.
 *
 * @note We use `pism_mask` (and not ice_thickness) to make decisions.
 * This means that we can update `ice_thickness` in place without
 * introducing a dependence on the grid traversal order.
 *
 * @param[in,out] pism_mask cell type mask
 * @param[in,out] ice_thickness modeled ice thickness
 *
 * @return 0 on success
 */
PetscErrorCode EigenCalving::remove_narrow_tongues(IceModelVec2Int &pism_mask,
                                                       IceModelVec2S &ice_thickness) {
  PetscErrorCode ierr;

  MaskQuery mask(pism_mask);

  ierr = pism_mask.begin_access(); CHKERRQ(ierr);
  ierr = ice_thickness.begin_access(); CHKERRQ(ierr);
  for (int   i = grid.xs; i < grid.xs + grid.xm; ++i) {
    for (int j = grid.ys; j < grid.ys + grid.ym; ++j) {
      if (mask.ice_free(i, j))  // FIXME: it might be better to have access to bedrock elevation b(i,j)
                                // and sea level SL so that the predicate can be
                                //   mask.ice_free(i,j) || (mask.grounded_ice(i,j) && (b(i,j) >= SL)))
        continue;

      const bool
        ice_free_N  = mask.ice_free_ocean(i, j + 1),
        ice_free_E  = mask.ice_free_ocean(i + 1, j),
        ice_free_S  = mask.ice_free_ocean(i, j - 1),
        ice_free_W  = mask.ice_free_ocean(i - 1, j),
        ice_free_NE = mask.ice_free_ocean(i + 1, j + 1),
        ice_free_NW = mask.ice_free_ocean(i - 1, j + 1),
        ice_free_SE = mask.ice_free_ocean(i + 1, j - 1),
        ice_free_SW = mask.ice_free_ocean(i - 1, j - 1);

      if ((ice_free_W == false &&
           ice_free_NW         &&
           ice_free_SW         &&
           ice_free_N          &&
           ice_free_S          &&
           ice_free_E)         ||
          (ice_free_N == false &&
           ice_free_NW         &&
           ice_free_NE         &&
           ice_free_W          &&
           ice_free_E          &&
           ice_free_S)         ||
          (ice_free_E == false &&
           ice_free_NE         &&
           ice_free_SE         &&
           ice_free_W          &&
           ice_free_S          &&
           ice_free_N)         ||
          (ice_free_S == false &&
           ice_free_SW         &&
           ice_free_SE         &&
           ice_free_W          &&
           ice_free_E          &&
           ice_free_N)) {
        pism_mask(i, j) = MASK_ICE_FREE_OCEAN;
        ice_thickness(i, j) = 0.0;
      }
    }
  }
  ierr = ice_thickness.end_access(); CHKERRQ(ierr);
  ierr = pism_mask.end_access(); CHKERRQ(ierr);

  return 0;
}

} // end of namespace pism
