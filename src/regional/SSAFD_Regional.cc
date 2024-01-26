/* Copyright (C) 2015, 2016, 2017, 2020, 2023, 2024 PISM Authors
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

#include "pism/regional/SSAFD_Regional.hh"
#include "pism/stressbalance/StressBalance.hh"

namespace pism {

namespace stressbalance {

SSAFD_Regional::SSAFD_Regional(std::shared_ptr<const Grid> g)
  : SSAFD(g, true) {

  m_h_stored = nullptr;
  m_H_stored = nullptr;
}

void SSAFD_Regional::init_impl() {

  SSAFD::init_impl();

  m_log->message(2, "  using the regional version of the SSA solver...\n");

  if (m_config->get_flag("stress_balance.ssa.dirichlet_bc")) {
    m_log->message(2, "  using stored SSA velocities as Dirichlet B.C. in the no_model_strip...\n");
  }
}

void SSAFD_Regional::solve(const Inputs &inputs) {
  m_h_stored = inputs.no_model_surface_elevation;
  m_H_stored = inputs.no_model_ice_thickness;

  SSAFD::solve(inputs);

  m_h_stored = nullptr;
  m_H_stored = nullptr;
}


SSA * SSAFD_RegionalFactory(std::shared_ptr<const Grid> grid) {
  return new SSAFD_Regional(grid);
}

} // end of namespace stressbalance

} // end of namespace pism
