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
#include "base/calving/PISMEigenCalving.hh"
#include "base/energy/bedrockThermalUnit.hh"
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

namespace pism {

//! Compute the maximum velocities for time-stepping and reporting to user.
/*!
Computes the maximum magnitude of the components \f$u,v,w\f$ of the 3D velocity.
Then sets `CFLmaxdt`, the maximum time step allowed under the
Courant-Friedrichs-Lewy (CFL) condition on the
horizontal advection scheme for age and for temperature.

Under BOMBPROOF there is no CFL condition for the vertical advection.
The maximum vertical velocity is computed but it does not affect
`CFLmaxdt`.
 */
double IceModel::max_timestep_cfl_3d() {
  double max_dt = m_config->get_double("maximum_time_step_years", "seconds");

  const IceModelVec3
    &u3 = m_stress_balance->velocity_u(),
    &v3 = m_stress_balance->velocity_v(),
    &w3 = m_stress_balance->velocity_w();

  IceModelVec::AccessList list;
  list.add(m_ice_thickness);
  list.add(u3);
  list.add(v3);
  list.add(w3);
  list.add(m_cell_type);

  // update global max of abs of velocities for CFL; only velocities under surface
  const double
    one_over_dx = 1.0 / m_grid->dx(),
    one_over_dy = 1.0 / m_grid->dy();

  double max_u = 0.0, max_v = 0.0, max_w = 0.0;
  ParallelSection loop(m_grid->com);
  try {
    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      if (m_cell_type.icy(i, j)) {
        const int ks = m_grid->kBelowHeight(m_ice_thickness(i, j));
        const double
          *u = u3.get_column(i, j),
          *v = v3.get_column(i, j),
          *w = w3.get_column(i, j);

        for (int k = 0; k <= ks; ++k) {
          const double
            absu = fabs(u[k]),
            absv = fabs(v[k]);
          max_u = std::max(max_u, absu);
          max_v = std::max(max_v, absv);
          const double denom = fabs(absu * one_over_dx) + fabs(absv * one_over_dy);
          if (denom > 0.0) {
            max_dt = std::min(max_dt, 1.0 / denom);
          }
        }

        for (int k = 0; k <= ks; ++k) {
          max_w = std::max(max_w, fabs(w[k]));
        }
      }
    }
  } catch (...) {
    loop.failed();
  }
  loop.check();

  m_max_u_speed = GlobalMax(m_grid->com, max_u);
  m_max_v_speed = GlobalMax(m_grid->com, max_v);
  m_max_w_speed = GlobalMax(m_grid->com, max_w);

  return GlobalMin(m_grid->com, max_dt);
}


//! Compute the CFL constant associated to first-order upwinding for the sliding contribution to mass continuity.
/*!
  This procedure computes the maximum horizontal speed in the icy
  areas. In particular it computes CFL constant for the upwinding, in
  massContExplicitStep(), which applies to the basal component of mass
  flux.

  That is, because the map-plane mass continuity is advective in the
  sliding case we have a CFL condition.
 */
double IceModel::max_timestep_cfl_2d() {
  double max_dt = m_config->get_double("maximum_time_step_years", "seconds");

  const IceModelVec2V &vel = m_stress_balance->advective_velocity();

  const double
    dx = m_grid->dx(),
    dy = m_grid->dy();

  IceModelVec::AccessList list;
  list.add(vel);
  list.add(m_cell_type);
  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    if (m_cell_type.icy(i, j)) {
      const double denom = fabs(vel(i, j).u) / dx + fabs(vel(i, j).v) / dy;
      if (denom > 0.0) {
        max_dt = std::min(max_dt, 1.0 / denom);
      }
    }
  }

  return GlobalMin(m_grid->com, max_dt);
}


//! Compute the maximum time step allowed by the diffusive SIA.
/*!
If maximum diffusivity is positive (i.e. if there is diffusion going on) then
updates dt.

Note adapt_ratio * 2 is multiplied by dx^2/(2*maxD) so dt <= adapt_ratio *
dx^2/maxD (if dx=dy).

Reference: [\ref MortonMayers] pp 62--63.
 */
double IceModel::max_timestep_diffusivity() {
  double D_max = m_stress_balance->max_diffusivity();

  if (D_max > 0.0) {
    const double
      dx = m_grid->dx(),
      dy = m_grid->dy(),
      adaptive_timestepping_ratio = m_config->get_double("adaptive_timestepping_ratio"),
      grid_factor                 = 1.0 / (dx*dx) + 1.0 / (dy*dy);

    return adaptive_timestepping_ratio * 2.0 / (D_max * grid_factor);
  } else {
    return m_config->get_double("maximum_time_step_years", "seconds");
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

  if (not m_config->get_boolean("do_skip")) {
    return 0;
  }

  const unsigned int skip_max = static_cast<int>(m_config->get_double("skip_max"));

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

  const bool update_3d = (skipCountDown == 0);
  const double current_time = m_time->current();
  const double time_to_end = m_time->end() - current_time;

  // FIXME: we should probably create a std::vector<const Component_TS*>
  // (or similar) and iterate over that instead.
  std::map<std::string, double> dt_restrictions;

  // Always consider the maximum allowed time-step length.
  if (m_config->get_double("maximum_time_step_years") > 0.0) {
    dt_restrictions["max"] = m_config->get_double("maximum_time_step_years", "seconds");
  }

  // Always consider maxdt_temporary.
  //
  // FIXME: maxdt_temporary is used by iceCompModel (and only there).
  // It should probably be removed.
  if (maxdt_temporary > 0.0) {
    dt_restrictions["internal (derived class)"] = maxdt_temporary;
  }

  // Never go past the end of a run.
  if (time_to_end > 0.0) {
    dt_restrictions["end of the run"] = time_to_end;
  }

  //! Always apply the time-step restriction from the
  //! -ts_{times,file,vars} mechanism (the user asked for it).
  MaxTimestep ts_dt = ts_max_timestep(current_time);
  if (ts_dt.is_finite()) {
    dt_restrictions["-ts_... reporting"] = ts_dt.value();
  }

  //! Always apply the time-step restriction from the
  //! -extra_{times,file,vars} mechanism (the user asked for it).
  MaxTimestep extras_dt = extras_max_timestep(current_time);
  if (extras_dt.is_finite()) {
    dt_restrictions["-extra_... reporting"] = extras_dt.value();
  }

  MaxTimestep save_dt = save_max_timestep(current_time);
  if (save_dt.is_finite()) {
    dt_restrictions["-save_... reporting"] = save_dt.value();
  }

  {
    // Query sub-models, which might add time-step restrictions.

    MaxTimestep surface_dt = m_surface->max_timestep(current_time);
    if (surface_dt.is_finite())  {
      dt_restrictions["surface"] = surface_dt.value();
    }

    MaxTimestep ocean_dt = m_ocean->max_timestep(current_time);
    if (ocean_dt.is_finite()) {
      dt_restrictions["ocean"] = ocean_dt.value();
    }

    MaxTimestep hydrology_dt = subglacial_hydrology->max_timestep(current_time);
    if (hydrology_dt.is_finite()) {
      dt_restrictions["hydrology"] = hydrology_dt.value();
    }

    if (btu != NULL) {
      MaxTimestep btu_dt = btu->max_timestep(current_time);
      if (btu_dt.is_finite()) {
        dt_restrictions["BTU"] = btu_dt.value();
      }
    }

    if (eigen_calving != NULL) {
      MaxTimestep eigencalving_dt = eigen_calving->max_timestep();
      if (eigencalving_dt.is_finite()) {
        dt_restrictions["eigencalving"] = eigencalving_dt.value();
      }
    }

    if (m_config->get_boolean("do_energy")) {
      if (update_3d) {
        CFLmaxdt = max_timestep_cfl_3d();
      }
      dt_restrictions["3D CFL"] = CFLmaxdt;
    }

    if (m_config->get_boolean("do_mass_conserve")) {
      CFLmaxdt2D = max_timestep_cfl_2d();

      dt_restrictions["2D CFL"] = CFLmaxdt2D;

      double max_dt_diffusivity = max_timestep_diffusivity();
      dt_restrictions["diffusivity"] = max_dt_diffusivity;
    }
  }

  // find the smallest of the max. time-steps reported by boundary models:
  {
    std::map<std::string, double>::const_iterator j = dt_restrictions.begin();
    if (dt_restrictions["max"] > 0.0) {
      m_adaptive_timestep_reason = "max";
    } else {
      m_adaptive_timestep_reason = j->first;
    }
    dt_result = dt_restrictions[m_adaptive_timestep_reason];
    while (j != dt_restrictions.end()) {
      if (j->second < dt_result) {
        dt_result = j->second;
        m_adaptive_timestep_reason = j->first;
      }
      ++j;
    }
  }

  // Hit multiples of X years, if requested (this has to go last):
  {
    const int timestep_hit_multiples = static_cast<int>(m_config->get_double("timestep_hit_multiples"));
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
        str << "(hit multiples of " << timestep_hit_multiples << " years; overrides " << m_adaptive_timestep_reason << ")";
        m_adaptive_timestep_reason = str.str();
      }
    }
  }


  if (dt_restrictions["diffusivity"] > 0.0) {
    double dt_diffusivity = dt_restrictions["diffusivity"];

    // remove the max. time-step corresponding to the diffusivity
    // stability criterion and find the smallest one again:
    dt_restrictions.erase("diffusivity");
    std::map<std::string, double>::const_iterator j = dt_restrictions.begin();
    double dt_other = j->second;
    std::string other_reason = j->first;
    while (j != dt_restrictions.end()) {
      if (j->second < dt_other) {
        dt_other = j->second;
        other_reason = j->first;
      }
      ++j;
    }

    if (dt_diffusivity < dt_other) {
      m_adaptive_timestep_reason += " (overrides " + other_reason + ")";
    }

    if (skip_counter_result == 0) {
      skip_counter_result = skip_counter(dt_other, dt_diffusivity);
    }
  }

  // "max", "internal (derived class)", and "end of the run" limit the "big" time-step (in
  // the context of the "skipping" mechanism), so we might need to
  // reset the skip_counter_result to 1.
  if ((m_adaptive_timestep_reason == "max" ||
       m_adaptive_timestep_reason == "internal (derived class)" ||
       m_adaptive_timestep_reason == "end of the run") &&
      skip_counter_result > 1) {
    skip_counter_result = 1;
  }
}


//! Because of the -skip mechanism it is still possible that we can have CFL violations: count them.
/*!
This applies to the horizontal part of the three-dimensional advection problem
solved by IceModel::ageStep() and the advection, ice-only part of the problem solved by
temperatureStep().  These methods use a fine vertical grid, and so we consider CFL
violations on that same fine grid. (FIXME: should we actually use the fine grid?)

Communication is needed to determine total CFL violation count over entire grid.
It is handled by temperatureAgeStep(), not here.
*/
unsigned int IceModel::countCFLViolations() {

  const double
    CFL_x = m_grid->dx() / dt_TempAge,
    CFL_y = m_grid->dy() / dt_TempAge;

  const IceModelVec3
    &u3 = m_stress_balance->velocity_u(),
    &v3 = m_stress_balance->velocity_v();

  IceModelVec::AccessList list;
  list.add(m_ice_thickness);
  list.add(u3);
  list.add(v3);

  unsigned int CFL_violation_count = 0;
  ParallelSection loop(m_grid->com);
  try {
    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      const int fks = m_grid->kBelowHeight(m_ice_thickness(i,j));

      const double
        *u = u3.get_column(i, j),
        *v = v3.get_column(i, j);

      // check horizontal CFL conditions at each point
      for (int k = 0; k <= fks; k++) {
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

  return (unsigned int)GlobalMax(m_grid->com, CFL_violation_count);
}

} // end of namespace pism
