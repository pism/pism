// Copyright (C) 2004-2016 Jed Brown, Ed Bueler and Constantine Khroulev
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

#include <petscvec.h>
#include <sstream>
#include <algorithm>

#include "iceModel.hh"
#include "base/calving/EigenCalving.hh"
#include "base/calving/vonMisesCalving.hh"
#include "base/energy/BedThermalUnit.hh"
#include "base/hydrology/PISMHydrology.hh"
#include "base/stressbalance/PISMStressBalance.hh"
#include "base/util/IceGrid.hh"
#include "base/util/Mask.hh"
#include "base/util/PISMConfigInterface.hh"
#include "base/util/PISMTime.hh"
#include "base/util/error_handling.hh"
#include "coupler/PISMOcean.hh"
#include "coupler/PISMSurface.hh"
#include "base/util/MaxTimestep.hh"
#include "base/util/pism_utilities.hh"
#include "base/age/AgeModel.hh"
#include "base/energy/EnergyModel.hh"

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
      adaptive_timestepping_ratio = m_config->get_double("time_stepping.adaptive_ratio"),
      grid_factor                 = 1.0 / (dx*dx) + 1.0 / (dy*dy);

    return MaxTimestep(adaptive_timestepping_ratio * 2.0 / (D_max * grid_factor),
                       "diffusivity");
  } else {
    return MaxTimestep(m_config->get_double("time_stepping.maximum_time_step", "seconds"),
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

  if (not m_config->get_boolean("time_stepping.skip.enabled")) {
    return 0;
  }

  const unsigned int skip_max = static_cast<int>(m_config->get_double("time_stepping.skip.max"));

  if (input_dt_diffusivity > 0.0) {
    const double conservativeFactor = 0.95;
    const double counter = floor(conservativeFactor * (input_dt / input_dt_diffusivity));
    const unsigned int result = static_cast<unsigned int>(counter);
    return std::min(result, skip_max);
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
by incorporating choices made by options (e.g. <c>-max_dt</c>) and by derived classes.

@param[out] dt_result computed maximum time step
@param[in,out] skip_counter_result time-step skipping counter
 */
void IceModel::max_timestep(double &dt_result, unsigned int &skip_counter_result) {

  const double current_time = m_time->current();

  // get time-stepping restrictions from sub-models
  std::vector<MaxTimestep> restrictions;
  {
    std::map<std::string, const Component*>::const_iterator k = m_submodels.begin();
    for (; k != m_submodels.end(); ++k) {
      const Component_TS* m = dynamic_cast<const Component_TS*>(k->second);
      if (m != NULL) {
        restrictions.push_back(m->max_timestep(current_time));
      }
    }
  }

  // Always consider the maximum allowed time-step length.
  if (m_config->get_double("time_stepping.maximum_time_step") > 0.0) {
    restrictions.push_back(MaxTimestep(m_config->get_double("time_stepping.maximum_time_step",
                                                            "seconds"),
                                       "max"));
  }

  // Never go past the end of a run.
  const double time_to_end = m_time->end() - current_time;
  if (time_to_end > 0.0) {
    restrictions.push_back(MaxTimestep(time_to_end, "end of the run"));
  }

  // reporting
  {
    restrictions.push_back(ts_max_timestep(current_time));
    restrictions.push_back(extras_max_timestep(current_time));
    restrictions.push_back(save_max_timestep(current_time));
  }

  // mass continuity stability criteria
  if (m_config->get_boolean("geometry.update.enabled")) {
    CFLData cfl = m_stress_balance->max_timestep_cfl_2d();

    restrictions.push_back(MaxTimestep(cfl.dt_max.value(), "2D CFL"));
    restrictions.push_back(max_timestep_diffusivity());
  }

  // Hit multiples of X years, if requested.
  {
    const int timestep_hit_multiples = static_cast<int>(m_config->get_double("time_stepping.hit_multiples"));
    if (timestep_hit_multiples > 0) {
      const double epsilon = 1.0; // 1 second tolerance
      double
        next_time = m_timestep_hit_multiples_last_time;

      while (m_time->increment_date(next_time, timestep_hit_multiples) <= current_time + dt_result + epsilon) {
        next_time = m_time->increment_date(next_time, timestep_hit_multiples);
      }

      if (next_time > current_time && next_time <= current_time + dt_result + epsilon) {
        dt_result = next_time - current_time;
        m_timestep_hit_multiples_last_time = next_time;

        std::stringstream str;
        str << "hit multiples of " << timestep_hit_multiples << " years";

        restrictions.push_back(MaxTimestep(next_time - current_time, str.str()));

      }
    }
  }

  // sort time step restrictions to find the strictest one
  std::sort(restrictions.begin(), restrictions.end());

  // note that restrictions has at least 2 elements
  // the first element is the max time step we can take
  MaxTimestep dt_max = restrictions[0];
  MaxTimestep dt_other = restrictions[1];
  dt_result = dt_max.value();
  m_adaptive_timestep_reason = (dt_max.description() +
                                " (overrides " + dt_other.description() + ")");

  // the "skipping" mechanism
  {
    if (dt_max.description() == "diffusivity" and skip_counter_result == 0) {
      skip_counter_result = skip_counter(dt_other.value(), dt_max.value());
    }

    // "max" and "end of the run" limit the "big" time-step (in
    // the context of the "skipping" mechanism), so we might need to
    // reset the skip_counter_result to 1.
    if ((m_adaptive_timestep_reason == "max" ||
         m_adaptive_timestep_reason == "end of the run") &&
        skip_counter_result > 1) {
      skip_counter_result = 1;
    }
  }
}


//! Because of the -skip mechanism it is still possible that we can have CFL violations: count them.
/*!
This applies to the horizontal part of the three-dimensional advection problem
solved by AgeModel and the advection, ice-only part of the problem solved by
temperatureStep().

Communication is needed to determine total CFL violation count over entire grid.
It is handled by temperatureAgeStep(), not here.
*/
unsigned int count_CFL_violations(const IceModelVec3 &u3,
                                  const IceModelVec3 &v3,
                                  const IceModelVec2S &ice_thickness,
                                  double dt) {


  IceGrid::ConstPtr grid = u3.get_grid();

  const double
    CFL_x = grid->dx() / dt,
    CFL_y = grid->dy() / dt;

  IceModelVec::AccessList list;
  list.add(ice_thickness);
  list.add(u3);
  list.add(v3);

  unsigned int CFL_violation_count = 0;
  ParallelSection loop(grid->com);
  try {
    for (Points p(*grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      const int ks = grid->kBelowHeight(ice_thickness(i,j));

      const double
        *u = u3.get_column(i, j),
        *v = v3.get_column(i, j);

      // check horizontal CFL conditions at each point
      for (int k = 0; k <= ks; k++) {
        if (fabs(u[k]) > CFL_x) {
          CFL_violation_count += 1;
        }
        if (fabs(v[k]) > CFL_y) {
          CFL_violation_count += 1;
        }
      }
    }
  } catch (...) {
    loop.failed();
  }
  loop.check();

  return (unsigned int)GlobalMax(grid->com, CFL_violation_count);
}

} // end of namespace pism
