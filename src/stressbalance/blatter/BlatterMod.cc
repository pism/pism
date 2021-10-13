/* Copyright (C) 2020 PISM Authors
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

#include <algorithm>           // std::min, std::max

#include "BlatterMod.hh"

#include "pism/rheology/FlowLawFactory.hh"

#include "pism/geometry/Geometry.hh"
#include "pism/stressbalance/StressBalance.hh" // Inputs
#include "pism/util/pism_utilities.hh" // GlobalMax

namespace pism {
namespace stressbalance {

BlatterMod::BlatterMod(std::shared_ptr<Blatter> solver)
  : SSB_Modifier(solver->grid()),
    m_solver(solver) {

  rheology::FlowLawFactory ice_factory("stress_balance.blatter.", m_config, m_EC);
  ice_factory.remove(ICE_GOLDSBY_KOHLSTEDT);
  m_flow_law = ice_factory.create();
}

void BlatterMod::init() {
  // empty
}

/*!
 * Post-process ice velocity computed by the Blatter solver.
 *
 * - transfers velocity from the sigma grid onto PISM's vertical grid
 *
 * - estimates the maximum diffusivity used to compute the time step restriction
 */
void BlatterMod::update(const IceModelVec2V &sliding_velocity,
                        const Inputs &inputs,
                        bool full_update) {
  (void) sliding_velocity;
  (void) full_update;

  // transfer velocity onto the PISM vertical grid, setting m_u and m_v
  transfer(inputs.geometry->ice_thickness);

  // estimate max diffusivity to use in adaptive time stepping
  compute_max_diffusivity(m_solver->velocity(),
                          inputs.geometry->ice_thickness,
                          inputs.geometry->ice_surface_elevation);

  m_diffusive_flux.set(0.0);
}

/*!
 * Copy ice velocity from the sigma vertical grid onto PISM's vertical grid.
 *
 * Uses constant extrapolation above the ice surface.
 */
void BlatterMod::transfer(const IceModelVec2S &ice_thickness) {

  auto u_sigma = m_solver->velocity_u_sigma();
  auto v_sigma = m_solver->velocity_v_sigma();

  const auto &zlevels = m_u.levels();
  int Mz = zlevels.size();

  IceModelVec::AccessList list{&m_u, &m_v, u_sigma.get(), v_sigma.get(), &ice_thickness};

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    auto u = m_u.get_column(i, j);
    auto v = m_v.get_column(i, j);

    double H = ice_thickness(i, j);

    if (H > 0.0) {
      for (int k = 0; k < Mz; ++k) {
        double sigma = std::min(zlevels[k] / H, 1.0);

        u[k] = u_sigma->interpolate(i, j, sigma);
        v[k] = v_sigma->interpolate(i, j, sigma);
      }
    } else {
      m_u.set_column(i, j, 0.0);
      m_v.set_column(i, j, 0.0);
    }
  }

  m_u.update_ghosts();
  m_v.update_ghosts();
}

/*!
 * Estimate max SIA-type diffusivity assuming that `Q = -D \nabla h`.
 */
void BlatterMod::compute_max_diffusivity(const IceModelVec2V &velocity,
                                         const IceModelVec2S &ice_thickness,
                                         const IceModelVec2S &surface) {
  double eps = 1e-3;
  double
    dx = m_grid->dx(),
    dy = m_grid->dy();

  IceModelVec::AccessList list{&velocity, &ice_thickness, &surface};

  m_D_max = 0.0;
  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    auto h = surface.star(i, j);
    auto H = ice_thickness(i, j);

    if (H > 0.0) {
      Vector2 grad_h = {(h.e - h.w) / (2.0 * dx),
                        (h.n - h.s) / (2.0 * dy)};

      double D = H * velocity(i, j).magnitude() / (grad_h.magnitude() + eps);

      m_D_max = std::max(D, m_D_max);
    }
  }

  m_D_max = GlobalMax(m_grid->com, m_D_max);
}

} // end of namespace stressbalance
} // end of namespace pism
