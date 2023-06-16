// Copyright (C) 2004-2017, 2019, 2020, 2021, 2022 Jed Brown, Ed Bueler and Constantine Khroulev
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

#include <algorithm>            // std::sort
#include <cmath>                // std::floor

#include "IceModel.hh"
#include "pism/util/Grid.hh"
#include "pism/util/ConfigInterface.hh"
#include "pism/util/Time.hh"
#include "pism/util/MaxTimestep.hh"
#include "pism/stressbalance/StressBalance.hh"
#include "pism/stressbalance/ShallowStressBalance.hh"
#include "pism/util/Component.hh" // ...->max_timestep()

#include "pism/frontretreat/calving/EigenCalving.hh"
#include "pism/frontretreat/calving/HayhurstCalving.hh"
#include "pism/frontretreat/calving/vonMisesCalving.hh"
#include "pism/frontretreat/FrontRetreat.hh"

#include "pism/energy/EnergyModel.hh"
#include "pism/coupler/OceanModel.hh"
#include "pism/coupler/FrontalMelt.hh"

namespace pism {

//! Compute the maximum time step allowed by the diffusive SIA.
/*!
If maximum diffusivity is positive (i.e. if there is diffusion going on) then
updates dt.

Note adapt_ratio * 2 is multiplied by dx^2/(2*maxD) so dt <= adapt_ratio *
dx^2/maxD (if dx=dy).

Reference: [\ref MortonMayers] pp 62--63.
 */
MaxTimestep IceModel::max_timestep_diffusivity() {
  double D_max = m_stress_balance->max_diffusivity();

  if (D_max > 0.0) {
    const double
      dx = m_grid->dx(),
      dy = m_grid->dy(),
      adaptive_timestepping_ratio = m_config->get_number("time_stepping.adaptive_ratio"),
      grid_factor                 = 1.0 / (dx*dx) + 1.0 / (dy*dy);

    return MaxTimestep(adaptive_timestepping_ratio * 2.0 / (D_max * grid_factor),
                       "diffusivity");
  } else {
    return MaxTimestep(m_config->get_number("time_stepping.maximum_time_step", "seconds"),
                       "max time step");
  }
}

/** @brief Compute the skip counter using "long" (usually determined
 * using the CFL stability criterion) and "short" (typically
 * determined using the diffusivity-based stability criterion) time
 * step lengths.
 *
 *
 * @param[in] input_dt long time-step
 * @param[in] input_dt_diffusivity short time-step
 *
 * @return new skip counter
 */
unsigned int IceModel::skip_counter(double input_dt, double input_dt_diffusivity) {

  if (not m_config->get_flag("time_stepping.skip.enabled")) {
    return 0;
  }

  const int skip_max = static_cast<int>(m_config->get_number("time_stepping.skip.max"));

  if (input_dt_diffusivity > 0.0) {
    const double conservativeFactor = 0.95;
    const double counter = floor(conservativeFactor * (input_dt / input_dt_diffusivity));
    return std::min(static_cast<int>(counter), skip_max);
  } else {
    return skip_max;
  }

  return 0;
}

//! Use various stability criteria to determine the time step for an evolution run.
/*!
The main loop in run() approximates many physical processes.  Several of these approximations,
including the mass continuity and temperature equations in particular, involve stability
criteria.  This procedure builds the length of the next time step by using these criteria and
by incorporating choices made by options (e.g. `-max_dt`) and by derived classes.

@param[in] counter current time-step skipping counter
 */
IceModel::TimesteppingInfo IceModel::max_timestep(unsigned int counter) {

  const double current_time = m_time->current();

  std::vector<MaxTimestep> restrictions;

  // get time-stepping restrictions from sub-models
  for (auto m : m_submodels) {
    restrictions.push_back(m.second->max_timestep(current_time));
  }

  // mechanisms that use a retreat rate
  bool front_retreat = (m_eigen_calving or m_vonmises_calving or
                        m_hayhurst_calving or m_frontal_melt);
  if (front_retreat and m_config->get_flag("geometry.front_retreat.use_cfl")) {
    // at least one of front retreat mechanisms is active *and* PISM is told to use a CFL
    // restriction

    array::Scalar &retreat_rate = *m_work2d[0];
    retreat_rate.set(0.0);

    if (m_eigen_calving) {
      retreat_rate.add(1.0, m_eigen_calving->calving_rate());
    }

    if (m_hayhurst_calving) {
      retreat_rate.add(1.0, m_hayhurst_calving->calving_rate());
    }

    if (m_vonmises_calving) {
      retreat_rate.add(1.0, m_vonmises_calving->calving_rate());
    }

    if (m_frontal_melt) {
      retreat_rate.add(1.0, m_frontal_melt->retreat_rate());
    }

    assert(m_front_retreat);

    restrictions.push_back(m_front_retreat->max_timestep(m_geometry.cell_type,
                                                         m_ice_thickness_bc_mask,
                                                         retreat_rate));
  }

  const char* end = "end of the run";
  const char* max = "max";

  // Always consider the maximum allowed time-step length.
  double max_timestep = m_config->get_number("time_stepping.maximum_time_step", "seconds");
  if (max_timestep > 0.0) {
    restrictions.push_back(MaxTimestep(max_timestep, max));
  }

  // Never go past the end of a run.
  const double time_to_end = m_time->end() - current_time;
  if (time_to_end > 0.0) {
    restrictions.push_back(MaxTimestep(time_to_end, end));
  }

  // reporting
  {
    restrictions.push_back(ts_max_timestep(current_time));
    restrictions.push_back(extras_max_timestep(current_time));
    restrictions.push_back(save_max_timestep(current_time));
  }

  // mass continuity stability criteria
  if (m_config->get_flag("geometry.update.enabled")) {
    auto cfl = m_stress_balance->max_timestep_cfl_2d();

    restrictions.push_back(MaxTimestep(cfl.dt_max.value(), "2D CFL"));
    restrictions.push_back(max_timestep_diffusivity());
  }

  // sort time step restrictions to find the strictest one
  std::sort(restrictions.begin(), restrictions.end());

  // note that restrictions has at least 2 elements
  // the first element is the max time step we can take
  auto dt_max = restrictions[0];
  auto dt_other = restrictions[1];

  TimesteppingInfo result;
  result.dt = dt_max.value();
  result.reason = (dt_max.description() + " (overrides " + dt_other.description() + ")");
  result.skip_counter = 0;

  double resolution = m_config->get_number("time_stepping.resolution", "seconds");

  // Hit all multiples of X years, if requested.
  {
    int year_increment = static_cast<int>(m_config->get_number("time_stepping.hit_multiples"));
    if (year_increment > 0) {
      auto next_time = m_time->increment_date(m_timestep_hit_multiples_last_time,
                                              year_increment);

      if (std::fabs(current_time - next_time) < resolution) {
        // the current time is a multiple of year_increment
        m_timestep_hit_multiples_last_time = current_time;
        next_time = m_time->increment_date(current_time, year_increment);
      }

      auto dt = next_time - current_time;
      assert(dt > resolution);

      if (dt < result.dt) {
        result.dt = dt;
        result.reason = pism::printf("hit multiples of %d years (overrides %s)",
                                     year_increment, dt_max.description().c_str());
      }
    }
  }

  // the "skipping" mechanism
  {
    if (dt_max.description() == "diffusivity" and counter == 0) {
      result.skip_counter = skip_counter(dt_other.value(), dt_max.value());
    } else {
      result.skip_counter = counter;
    }

    // "max" and "end of the run" limit the "big" time-step (in
    // the context of the "skipping" mechanism), so we might need to
    // reset the skip_counter_result to 1.
    if (member(dt_max.description(), {max, end}) and counter > 1) {
      result.skip_counter = 1;
    }
  }

  if (resolution > 0.0) {
    double dt = std::floor(result.dt * resolution) / resolution;

    // Ensure that the resulting time step is never zero. This may happen if the length of
    // the run is not an integer multiple of "resolution".
    if (dt >= resolution) {
      result.dt = dt;
    }
  }

  return result;
}

} // end of namespace pism
