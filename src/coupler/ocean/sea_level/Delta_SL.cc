// Copyright (C) 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018 PISM Authors
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

#include "Delta_SL.hh"
#include "pism/coupler/util/ScalarForcing.hh"

namespace pism {
namespace ocean {
namespace sea_level {

Delta_SL::Delta_SL(IceGrid::ConstPtr grid, std::shared_ptr<SeaLevel> in)
  : SeaLevel(grid, in) {

  m_forcing.reset(new ScalarForcing(grid->ctx(),
                                    "-ocean_delta_sl", "delta_SL",
                                    "m", "m",
                                    "sea level elevation offsets"));
}

Delta_SL::~Delta_SL() {
  // empty
}

void Delta_SL::init_impl(const Geometry &geometry) {

  m_input_model->init(geometry);

  m_log->message(2, "* Initializing scalar sea level forcing...\n");

  m_forcing->init();
}

void Delta_SL::update_impl(const Geometry &geometry, double t, double dt) {
  m_input_model->update(geometry, t, dt);

  m_forcing->update(t, dt);

  m_sea_level.copy_from(m_input_model->elevation());
  m_sea_level.shift(m_forcing->value());
}

} // end of namespace sea_level
} // end of namespace ocean
} // end of namespace pism
