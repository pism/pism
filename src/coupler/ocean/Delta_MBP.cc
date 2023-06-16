/* Copyright (C) 2021 PISM Authors
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

#include "Delta_MBP.hh"
#include "pism/util/ScalarForcing.hh"
#include "pism/geometry/Geometry.hh"

namespace pism {
namespace ocean {

Delta_MBP::Delta_MBP(std::shared_ptr<const IceGrid> g, std::shared_ptr<OceanModel> in)
  : OceanModel(g, in) {

  m_forcing.reset(new ScalarForcing(*g->ctx(),
                                    "ocean.delta_MBP",
                                    "delta_MBP",
                                    "Pa", "Pa",
                                    "melange back pressure"));

  m_water_column_pressure = allocate_water_column_pressure(g);
}

Delta_MBP::~Delta_MBP() {
  // empty
}

void Delta_MBP::init_impl(const Geometry &geometry) {

  m_input_model->init(geometry);

  m_log->message(2,
                 "* Initializing melange back pressure forcing using scalar offsets...\n");
}

void Delta_MBP::update_impl(const Geometry &geometry, double t, double dt) {
  m_input_model->update(geometry, t, dt);

  double
    melange_thickness = m_config->get_number("ocean.delta_MBP.melange_thickness"),
    dP_melange        = m_forcing->value(t + 0.5 * dt);

  const auto &P = m_input_model->average_water_column_pressure();
  const auto &H = geometry.ice_thickness;
  auto &P_new = *m_water_column_pressure;
  array::AccessScope list{&P, &H, &P_new};

  for (auto p = m_grid->points(); p; p.next()) {
    const int i = p.i(), j = p.j();

    double dP = H(i, j) > 0.0 ? (melange_thickness * dP_melange) / H(i, j) : 0.0;

    P_new(i, j) = P(i, j) + dP;
  }
}

const array::Scalar& Delta_MBP::average_water_column_pressure_impl() const {
  return *m_water_column_pressure;
}


} // end of namespace ocean
} // end of namespace pism
