/* Copyright (C) 2016, 2017 PISM Authors
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

#include "timestepping.hh"
#include "pism/util/IceGrid.hh"
#include "pism/util/iceModelVec.hh"
#include "pism/util/IceModelVec2CellType.hh"
#include "pism/util/pism_utilities.hh"

namespace pism {

CFLData::CFLData() {
  u_max = 0.0;
  v_max = 0.0;
  w_max = 0.0;
}

//! Compute the maximum velocities for time-stepping and reporting to user.
/*!
Computes the maximum magnitude of the components \f$u,v,w\f$ of the 3D velocity.

Under BOMBPROOF there is no CFL condition for the vertical advection.
The maximum vertical velocity is computed but it does not affect the output.
 */
CFLData max_timestep_cfl_3d(const IceModelVec2S &ice_thickness,
                            const IceModelVec2CellType &cell_type,
                            const IceModelVec3 &u3,
                            const IceModelVec3 &v3,
                            const IceModelVec3 &w3) {

  IceGrid::ConstPtr grid = ice_thickness.grid();
  Config::ConstPtr config = grid->ctx()->config();

  double dt_max = config->get_double("time_stepping.maximum_time_step", "seconds");

  IceModelVec::AccessList list{&ice_thickness, &u3, &v3, &w3, &cell_type};

  // update global max of abs of velocities for CFL; only velocities under surface
  const double
    one_over_dx = 1.0 / grid->dx(),
    one_over_dy = 1.0 / grid->dy();

  double u_max = 0.0, v_max = 0.0, w_max = 0.0;
  ParallelSection loop(grid->com);
  try {
    for (Points p(*grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      if (cell_type.icy(i, j)) {
        const int ks = grid->kBelowHeight(ice_thickness(i, j));
        const double
          *u = u3.get_column(i, j),
          *v = v3.get_column(i, j),
          *w = w3.get_column(i, j);

        for (int k = 0; k <= ks; ++k) {
          const double
            u_abs = fabs(u[k]),
            v_abs = fabs(v[k]);
          u_max = std::max(u_max, u_abs);
          v_max = std::max(v_max, v_abs);
          const double denom = fabs(u_abs * one_over_dx) + fabs(v_abs * one_over_dy);
          if (denom > 0.0) {
            dt_max = std::min(dt_max, 1.0 / denom);
          }
        }

        for (int k = 0; k <= ks; ++k) {
          w_max = std::max(w_max, fabs(w[k]));
        }
      }
    }
  } catch (...) {
    loop.failed();
  }
  loop.check();

  CFLData result;

  result.u_max = GlobalMax(grid->com, u_max);
  result.v_max = GlobalMax(grid->com, v_max);
  result.w_max = GlobalMax(grid->com, w_max);
  result.dt_max = MaxTimestep(GlobalMin(grid->com, dt_max));

  return result;
}

//! Compute the CFL constant associated to first-order upwinding for the sliding contribution to mass continuity.
/*!
  This procedure computes the maximum horizontal speed in the icy
  areas. In particular it computes CFL constant for the upwinding, in
  GeometryEvolution::step(), which applies to the basal component of mass
  flux.

  That is, because the map-plane mass continuity is advective in the
  sliding case we have a CFL condition.
 */
CFLData max_timestep_cfl_2d(const IceModelVec2S &ice_thickness,
                            const IceModelVec2CellType &cell_type,
                            const IceModelVec2V &velocity) {

  IceGrid::ConstPtr grid = ice_thickness.grid();
  Config::ConstPtr config = grid->ctx()->config();

  double dt_max = config->get_double("time_stepping.maximum_time_step", "seconds");

  const double
    dx = grid->dx(),
    dy = grid->dy();

  IceModelVec::AccessList list{&velocity, &cell_type};

  double u_max = 0.0, v_max = 0.0;
  for (Points p(*grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    if (cell_type.icy(i, j)) {
      const double
        u_abs = fabs(velocity(i, j).u),
        v_abs = fabs(velocity(i, j).v);

      u_max = std::max(u_max, u_abs);
      v_max = std::max(v_max, v_abs);

      const double denom = u_abs / dx + v_abs / dy;
      if (denom > 0.0) {
        dt_max = std::min(dt_max, 1.0 / denom);
      }
    }
  }

  CFLData result;

  result.u_max = GlobalMax(grid->com, u_max);
  result.v_max = GlobalMax(grid->com, v_max);
  result.w_max = 0.0;
  result.dt_max = MaxTimestep(GlobalMin(grid->com, dt_max));

  return result;
}

} // end of namespace pism
