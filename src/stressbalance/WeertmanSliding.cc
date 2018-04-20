/* Copyright (C) 2018 PISM Authors
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

#include "WeertmanSliding.hh"

#include "pism/rheology/FlowLawFactory.hh"
#include "pism/geometry/Geometry.hh"
#include "StressBalance.hh"

namespace pism {
namespace stressbalance {

WeertmanSliding::WeertmanSliding(IceGrid::ConstPtr grid)
  : ShallowStressBalance(grid) {
  // Use the SIA flow law.
  rheology::FlowLawFactory ice_factory("stress_balance.sia.", m_config, m_EC);
  m_flow_law = ice_factory.create();
}

WeertmanSliding::~WeertmanSliding() {
  // empty
}

void WeertmanSliding::init_impl() {
  m_log->message(2, "* Initializing Weertman-style basal sliding...\n");
}


/*!
 * Compute sliding velocity using a Weertman-style parameterization from [@ref Tomkin2007],
 * equation 5.
 *
 * @f[ u_s = \frac{2 A_s \beta_c (\rho g H)^{n}}{N - P} \cdot |\nabla h|^{n-1} \cdot \nabla h,
 * @f]
 *
 * where
 *
 * - @f$ A_s @f$ is the sliding parameter,
 * - @f$ \beta_c @f$ is a parameter capturing the effect of the presence of valley walls
 *   (set to 1 in this implementation),
 * - @f$ \rho @f$ is the ice density,
 * - @f$ H @f$ is the ice thickness,
 * - @f$ h @f$ is the ice surface elevation ,
 * - @f$ n @f$ is the flow law exponent (usually 3),
 * - @f$ g @f$ is the acceleration due to gravity,
 * - @f$ N @f$ is the ice overburden pressure,
 * - @f$ P @f$ is the basal water pressure.
 *
 * With these modifications and noting that @f$ N = \rho g H @f$, the formula above
 * becomes
 *
 * @f[ u_s = \frac{2 A_s}{1 - k} \cdot (N |\nabla h|)^{n-1} \cdot \nabla h,
 * @f]
 *
 * where
 * - @f$ N = \rho g H @f$,
 * - @f$ P = k N @f$ (we assume that basal water pressure is a given fraction of overburden)
 *
 * This parameterization is used for areas of grounded ice where the base of the ice is
 * temperate.
 */
void WeertmanSliding::update(const Inputs &inputs, bool full_update) {

  (void) full_update;

  const IceModelVec2S        &H         = inputs.geometry->ice_thickness;
  const IceModelVec2S        &h         = inputs.geometry->ice_surface_elevation;
  const IceModelVec2CellType &cell_type = inputs.geometry->cell_type;
  const IceModelVec3         &enthalpy  = *inputs.enthalpy;

  double n   = m_flow_law->exponent();
  double A_s = m_config->get_double("stress_balance.weertman_sliding.A");
  double k   = m_config->get_double("stress_balance.weertman_sliding.k");

  IceModelVec::AccessList list{&m_velocity, &H, &h, &enthalpy, &cell_type};

  ParallelSection loop(m_grid->com);
  try {
    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      double
        P_o    = m_EC->pressure(H(i, j)),
        E_base = enthalpy(i, j, 0);

      if (not m_EC->is_temperate(E_base, P_o) or cell_type.ocean(i, j)) {
        m_velocity(i, j) = {0.0, 0.0};
        continue;
      }

      // Note: we may need to decide if we should use one-sided FD at ice margins.
      Vector2 grad_h = {h.diff_x_p(i, j), h.diff_y_p(i, j)};

      // Note: this could be optimized by computing this instead
      // 2 * A_s / (1 - k) * pow(P * P * (h_x * h_x + h_y * h_y), (n - 1) / 2) * grad_h;
      // ... but I'm not sure we need to and the current code is cleaner.
      m_velocity(i, j) = 2.0 * A_s / (1.0 - k) * pow(P_o * grad_h.magnitude(), n - 1) * grad_h;
    }
  } catch (...) {
    loop.failed();
  }
  loop.check();

}

} // end of namespace stressbalance
} // end of namespace pism
