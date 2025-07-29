/* Copyright (C) 2020, 2021, 2023 PISM Authors
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

#include "pism/stressbalance/molho/MOLHOMod.hh"

#include "pism/rheology/FlowLawFactory.hh"

#include "pism/geometry/Geometry.hh"
#include "pism/stressbalance/StressBalance.hh" // Inputs
#include "pism/util/pism_utilities.hh" // GlobalMax

namespace pism {
namespace stressbalance {

MOLHOMod::MOLHOMod(std::shared_ptr<MOLHO> solver)
  : SSB_Modifier(solver->grid()),
    m_solver(solver) {

  rheology::FlowLawFactory ice_factory("stress_balance.blatter.", m_config, m_EC);
  ice_factory.remove(ICE_GOLDSBY_KOHLSTEDT);
  m_flow_law = ice_factory.create();
}

void MOLHOMod::init() {
  // empty
}

/*!
 * Post-process ice velocity computed by the MOLHO solver.
 *
 * - estimates the maximum diffusivity used to compute the time step restriction
 */
void MOLHOMod::update(const array::Vector &sliding_velocity,
                        const Inputs &inputs,
                        bool full_update) {
  (void) sliding_velocity;
  (void) full_update;

  // estimate max diffusivity to use in adaptive time stepping
  compute_max_diffusivity(m_solver->velocity(),
                          inputs.geometry->ice_thickness,
                          inputs.geometry->ice_surface_elevation);

  m_diffusive_flux.set(0.0);
}

/*!
 * Estimate max SIA-type diffusivity assuming that `Q = -D \nabla h`.
 */
void MOLHOMod::compute_max_diffusivity(const array::Vector &velocity,
                                         const array::Scalar &ice_thickness,
                                         const array::Scalar1 &surface) {
  const double eps = 1e-3;
  double
    dx = m_grid->dx(),
    dy = m_grid->dy();

  array::AccessScope list{&velocity, &ice_thickness, &surface};

  m_D_max = 0.0;
  for (auto p = m_grid->points(); p; p.next()) {
    const int i = p.i(), j = p.j();

    auto h = surface.star(i, j);
    auto H = ice_thickness(i, j);

    if (H > 0.0) {
      Vector2d grad_h = {(h.e - h.w) / (2.0 * dx),
                        (h.n - h.s) / (2.0 * dy)};

      double D = H * velocity(i, j).magnitude() / (grad_h.magnitude() + eps);

      m_D_max = std::max(D, m_D_max);
    }
  }

  m_D_max = GlobalMax(m_grid->com, m_D_max);
}

} // end of namespace stressbalance
} // end of namespace pism
