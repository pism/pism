/* Copyright (C) 2013, 2014, 2015, 2016, 2017, 2018 PISM Authors
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
#include "pism/coupler/util/ScalarForcing.hh"

namespace pism {
namespace ocean {

Frac_MBP::Frac_MBP(IceGrid::ConstPtr g, std::shared_ptr<OceanModel> in)
  : OceanModel(g, in) {

  m_forcing.reset(new ScalarForcing(g->ctx(),
                                    "-ocean_frac_MBP",
                                    "frac_MBP",
                                    "1", "1",
                                    "melange back pressure fraction"));

  m_melange_back_pressure_fraction = allocate_melange_back_pressure(g);
}

Frac_MBP::~Frac_MBP() {
  // empty
}

void Frac_MBP::init_impl(const Geometry &geometry) {

  m_input_model->init(geometry);

  m_log->message(2, "* Initializing melange back pressure fraction forcing...\n");

  m_forcing->init();
}

void Frac_MBP::update_impl(const Geometry &geometry, double t, double dt) {
  m_input_model->update(geometry, t, dt);

  m_forcing->update(t, dt);

  m_melange_back_pressure_fraction->copy_from(m_input_model->melange_back_pressure_fraction());
  m_melange_back_pressure_fraction->scale(m_forcing->value());
}

const IceModelVec2S& Frac_MBP::melange_back_pressure_fraction_impl() const {
  return *m_melange_back_pressure_fraction;
}


} // end of namespace ocean
} // end of namespace pism
