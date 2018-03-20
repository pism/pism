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

/// -ocean_delta_SL_file, ...

Delta_SL::Delta_SL(IceGrid::ConstPtr g, std::shared_ptr<OceanModel> in)
  : OceanModel(g, in) {

  m_forcing.reset(new ScalarForcing(g->ctx(),
                                    "-ocean_delta_SL", "delta_SL",
                                    "m", "m",
                                    "sea level elevation offsets"));

  m_sea_level_elevation = allocate_sea_level_elevation(g);
}

Delta_SL::~Delta_SL() {
  // empty
}

void Delta_SL::init_impl(const Geometry &geometry) {

  m_input_model->init(geometry);

  m_log->message(2, "* Initializing sea level forcing...\n");

  m_forcing->init();
}

void Delta_SL::update_impl(const Geometry &geometry, double t, double dt) {
  m_input_model->update(geometry, t, dt);

  m_forcing->update(t, dt);

  m_sea_level_elevation->copy_from(m_input_model->sea_level_elevation());
  m_sea_level_elevation->shift(m_forcing->value());
}

const IceModelVec2S& Delta_SL::sea_level_elevation_impl() const {
  return *m_sea_level_elevation;
}

} // end of namespace ocean
} // end of namespace pism
