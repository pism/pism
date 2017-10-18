/* Copyright (C) 2015, 2016, 2017 PISM Authors
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

#include "SSAFD_Regional.hh"
#include "pism/util/Vars.hh"
#include "pism/stressbalance/StressBalance.hh"

namespace pism {

namespace stressbalance {

SSAFD_Regional::SSAFD_Regional(IceGrid::ConstPtr g)
  : SSAFD(g) {

  m_h_stored      = NULL;
  m_H_stored      = NULL;
  m_no_model_mask = NULL;
}

SSAFD_Regional::~SSAFD_Regional() {
  // empty
}

void SSAFD_Regional::init() {

  SSAFD::init();

  m_log->message(2, "  using the regional version of the SSA solver...\n");

  if (m_config->get_boolean("stress_balance.ssa.dirichlet_bc")) {
    m_log->message(2, "  using stored SSA velocities as Dirichlet B.C. in the no_model_strip...\n");
  }
}

void SSAFD_Regional::update(const Inputs &inputs, bool full_update) {
  m_h_stored      = inputs.no_model_surface_elevation;
  m_H_stored      = inputs.no_model_ice_thickness;
  m_no_model_mask = inputs.no_model_mask;

  SSA::update(inputs, full_update);

  m_h_stored      = NULL;
  m_H_stored      = NULL;
  m_no_model_mask = NULL;
}

void SSAFD_Regional::compute_driving_stress(const Geometry &geometry, IceModelVec2V &result) const {

  SSAFD::compute_driving_stress(geometry, result);

  const IceModelVec2Int &nmm = *m_no_model_mask;

  const IceModelVec2S
    &usurfstore = *m_h_stored,
    &thkstore   = *m_H_stored;

  IceModelVec::AccessList list{&result, &nmm, &usurfstore, &thkstore};

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    double pressure = m_EC->pressure(thkstore(i,j));
    if (pressure <= 0) {
      pressure = 0;
    }

    if (nmm(i, j) > 0.5 || nmm(i - 1, j) > 0.5 || nmm(i + 1, j) > 0.5) {
      if (i - 1 < 0 || i + 1 > (int)m_grid->Mx() - 1) {
        result(i, j).u = 0;
      } else {
        result(i, j).u = - pressure * usurfstore.diff_x(i,j);
      }
    }

    if (nmm(i, j) > 0.5 || nmm(i, j - 1) > 0.5 || nmm(i, j + 1) > 0.5) {
      if (j - 1 < 0 || j + 1 > (int)m_grid->My() - 1) {
        result(i, j).v = 0;
      } else {
        result(i, j).v = - pressure * usurfstore.diff_y(i,j);
      }
    }
  }
}

} // end of namespace stressbalance

} // end of namespace pism
