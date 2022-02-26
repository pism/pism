/* Copyright (C) 2013, 2014, 2015, 2016, 2017, 2018, 2019, 2021 PISM Authors
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

#include "Frac_MBP.hh"
#include "pism/util/ScalarForcing.hh"
#include "pism/geometry/Geometry.hh"
#include "pism/util/pism_utilities.hh"

namespace pism {
namespace ocean {

Frac_MBP::Frac_MBP(IceGrid::ConstPtr g, std::shared_ptr<OceanModel> in)
  : OceanModel(g, in) {

  m_forcing.reset(new ScalarForcing(*g->ctx(),
                                    "ocean.frac_MBP",
                                    "frac_MBP",
                                    "1", "1",
                                    "melange back pressure fraction"));

  m_water_column_pressure = allocate_water_column_pressure(g);
}

Frac_MBP::~Frac_MBP() {
  // empty
}

void Frac_MBP::init_impl(const Geometry &geometry) {

  m_input_model->init(geometry);

  m_log->message(2,
                 "* Initializing melange back pressure fraction forcing...\n");
}


/*!
 * Modify water column pressure using a scalar factor `\lambda`
 *
 * We assume that the "melange back pressure" `P_m = \lambda (P_i - P_o)`, where `P_i` is
 * the vertically-integrated cryostatic pressure of an ice column and, similarly, `P_o` is
 * the vertically-integrated pressure of a water column and `\lambda \in [0, 1]`.
 *
 * Then the modified integrated water column pressure is
 *
 * `P_o^{*} = P_o + P_m`.
 */
void Frac_MBP::update_impl(const Geometry &geometry, double t, double dt) {
  m_input_model->update(geometry, t, dt);

  m_water_column_pressure->copy_from(m_input_model->average_water_column_pressure());

  double
    lambda      = m_forcing->value(t + 0.5 * dt),
    ice_density = m_config->get_number("constants.ice.density"),
    g           = m_config->get_number("constants.standard_gravity");

  array::Scalar &P_o = *m_water_column_pressure;

  IceModelVec::AccessList list{&P_o, &geometry.ice_thickness};

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    double
      P_i = 0.5 * ice_density * g * geometry.ice_thickness(i, j), // vertical average
      P_m = lambda * (P_i - P_o(i, j));

    P_o(i, j) += P_m;
  }
}

const array::Scalar& Frac_MBP::average_water_column_pressure_impl() const {
  return *m_water_column_pressure;
}


} // end of namespace ocean
} // end of namespace pism
