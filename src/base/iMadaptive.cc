// Copyright (C) 2004-2015 Jed Brown, Ed Bueler and Constantine Khroulev
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
#include "Mask.hh"
#include "PISMStressBalance.hh"
#include "bedrockThermalUnit.hh"
#include "PISMTime.hh"
#include "PISMEigenCalving.hh"
#include "PISMOcean.hh"
#include "PISMSurface.hh"
#include "PISMHydrology.hh"
#include "IceGrid.hh"

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
void IceModel::max_timestep_cfl_3d(double &dt_result) {
  double maxtimestep = config.get("maximum_time_step_years", "years", "seconds");

  IceModelVec3 *u3 = NULL, *v3 = NULL, *w3 = NULL;
  stress_balance->get_3d_velocity(u3, v3, w3);

  IceModelVec::AccessList list;
  list.add(ice_thickness);
  list.add(*u3);
  list.add(*v3);
  list.add(*w3);
  list.add(vMask);

  MaskQuery mask(vMask);
  double *u = NULL, *v = NULL, *w = NULL;

  // update global max of abs of velocities for CFL; only velocities under surface
  double maxu = 0.0, maxv = 0.0, maxw = 0.0;
  for (Points p(grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    if (mask.icy(i, j)) {
      const int ks = grid.kBelowHeight(ice_thickness(i, j));
      u3->getInternalColumn(i, j, &u);
      v3->getInternalColumn(i, j, &v);
      w3->getInternalColumn(i, j, &w);
      for (int k = 0; k <= ks; ++k) {
        const double
          absu = fabs(u[k]),
          absv = fabs(v[k]);
        maxu = std::max(maxu, absu);
        maxv = std::max(maxv, absv);
        maxw = std::max(maxw, fabs(w[k]));
        const double denom = fabs(absu / grid.dx()) + fabs(absv / grid.dy());
        if (denom > 0.0) {
          maxtimestep = std::min(maxtimestep, 1.0 / denom);
        }
      }
    }
  }


  gmaxu = GlobalMax(grid.com, maxu);
  gmaxv = GlobalMax(grid.com, maxv);
  gmaxw = GlobalMax(grid.com, maxw);

  dt_result = GlobalMin(grid.com, maxtimestep);
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
void IceModel::max_timestep_cfl_2d(double &dt_result) {
  double maxtimestep = config.get("maximum_time_step_years", "years", "seconds");

  MaskQuery mask(vMask);

  IceModelVec2V *vel_advective;
  stress_balance->get_2D_advective_velocity(vel_advective);
  IceModelVec2V &vel = *vel_advective; // a shortcut

  IceModelVec::AccessList list;
  list.add(vel);
  list.add(vMask);
  for (Points p(grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    if (mask.icy(i, j)) {
      const double denom = fabs(vel(i,j).u)/grid.dx() + fabs(vel(i,j).v)/grid.dy();
      if (denom > 0.0) {
        maxtimestep = std::min(maxtimestep, 1.0/denom);
      }
    }
  }

  dt_result = GlobalMin(grid.com, maxtimestep);
}


//! Compute the maximum time step allowed by the diffusive SIA.
/*!
If maximum diffusivity is positive (i.e. if there is diffusion going on) then
updates dt.

Note adapt_ratio * 2 is multiplied by dx^2/(2*maxD) so dt <= adapt_ratio *
dx^2/maxD (if dx=dy).

Reference: [\ref MortonMayers] pp 62--63.
 */
void IceModel::max_timestep_diffusivity(double &dt_result) {
  double D_max = 0.0;
  stress_balance->get_max_diffusivity(D_max);

  if (D_max > 0.0) {
    const double
      adaptive_timestepping_ratio = config.get("adaptive_timestepping_ratio"),
      grid_factor                 = 1.0 / (grid.dx()*grid.dx()) + 1.0 / (grid.dy()*grid.dy());

    dt_result = adaptive_timestepping_ratio * 2.0 / (D_max * grid_factor);
  } else {
    dt_result = config.get("maximum_time_step_years", "years", "seconds");
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

  if (config.get_flag("do_skip") == false) {
    return 0;
  }

  const unsigned int skip_max = static_cast<int>(config.get("skip_max"));

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

  const bool updateAtDepth = (skipCountDown == 0);
  const double time_to_end = grid.time->end() - grid.time->current();

  // FIXME: we should probably create a std::vector<const Component_TS*>
  // (or similar) and iterate over that instead.
  bool restrict_dt = false;
  const double current_time = grid.time->current();
  std::map<std::string, double> dt_restrictions;

  // Always consider the maximum allowed time-step length.
  if (config.get("maximum_time_step_years") > 0.0) {
    dt_restrictions["max"] = config.get("maximum_time_step_years", "years", "seconds");
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
  double ts_dt = 0.0;
  ts_max_timestep(current_time, ts_dt, restrict_dt);
  if (restrict_dt) {
    dt_restrictions["-ts_... reporting"] = ts_dt;
  }

  //! Always apply the time-step restriction from the
  //! -extra_{times,file,vars} mechanism (the user asked for it).
  double extras_dt = 0.0;
  extras_max_timestep(current_time, extras_dt, restrict_dt);
  if (restrict_dt) {
    dt_restrictions["-extra_... reporting"] = extras_dt;
  }

  if (dt_force > 0.0) {
    dt_restrictions["fixed"] = dt_force;
    // If the user asked for fixed time-steps, we're done; proceed to
    // comparing time-step restrictions.
  } else {

    // ... else query sub-models, which might add more time-step
    // restrictions.

    double surface_dt = 0.0;
    surface->max_timestep(current_time, surface_dt, restrict_dt);
    if (restrict_dt)  {
      dt_restrictions["surface"] = surface_dt;
    }

    double ocean_dt = 0.0;
    ocean->max_timestep(current_time, ocean_dt, restrict_dt);
    if (restrict_dt) {
      dt_restrictions["ocean"] = ocean_dt;
    }

    double hydrology_dt = 0.0;
    subglacial_hydrology->max_timestep(current_time, hydrology_dt, restrict_dt);
    if (restrict_dt) {
      dt_restrictions["hydrology"] = hydrology_dt;
    }

    if (btu != NULL) {
      double btu_dt = 0.0;
      btu->max_timestep(current_time, btu_dt, restrict_dt);
      if (restrict_dt) {
        dt_restrictions["BTU"] = btu_dt;
      }
    }

    if (eigen_calving != NULL) {
      double eigencalving_dt = 0.0;
      eigen_calving->max_timestep(current_time, eigencalving_dt, restrict_dt);
      if (restrict_dt) {
        dt_restrictions["eigencalving"] = eigencalving_dt;
      }
    }

    if (config.get_flag("do_energy") == true) {
      if (updateAtDepth == true) {
        // this call sets CFLmaxdt
        max_timestep_cfl_3d(CFLmaxdt);
      }
      dt_restrictions["3D CFL"] = CFLmaxdt;
    }

    if (config.get_flag("do_mass_conserve")) {
      // this call sets CFLmaxdt2D
      max_timestep_cfl_2d(CFLmaxdt2D);

      dt_restrictions["2D CFL"] = CFLmaxdt2D;

      double max_dt_diffusivity = 0.0;
      max_timestep_diffusivity(max_dt_diffusivity);
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
    const int timestep_hit_multiples = static_cast<int>(config.get("timestep_hit_multiples"));
    if (timestep_hit_multiples > 0) {
      const double epsilon = 1.0; // 1 second tolerance
      double
        next_time = timestep_hit_multiples_last_time;

      while (grid.time->increment_date(next_time, timestep_hit_multiples) <= current_time + dt_result + epsilon) {
        next_time = grid.time->increment_date(next_time, timestep_hit_multiples);
      }

      if (next_time > current_time && next_time <= current_time + dt_result + epsilon) {
        dt_result = next_time - current_time;
        timestep_hit_multiples_last_time = next_time;

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
void IceModel::countCFLViolations(double* CFLviol) {

  const double
    CFL_x = grid.dx() / dt_TempAge,
    CFL_y = grid.dy() / dt_TempAge;

  double *u, *v;
  IceModelVec3 *u3, *v3, *dummy;
  stress_balance->get_3d_velocity(u3, v3, dummy);

  IceModelVec::AccessList list;
  list.add(ice_thickness);
  list.add(*u3);
  list.add(*v3);

  for (Points p(grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    const int  fks = grid.kBelowHeight(ice_thickness(i,j));

    u3->getInternalColumn(i,j,&u);
    v3->getInternalColumn(i,j,&v);

    // check horizontal CFL conditions at each point
    for (int k=0; k<=fks; k++) {
      if (fabs(u[k]) > CFL_x) {
        *CFLviol += 1.0;
      }
      if (fabs(v[k]) > CFL_y) {
        *CFLviol += 1.0;
      }
    }
  }
}

} // end of namespace pism
