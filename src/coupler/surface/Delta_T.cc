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

#include "Delta_T.hh"
#include "pism/coupler/util/ScalarForcing.hh"

namespace pism {
namespace surface {

/// -surface ...,delta_T (scalar forcing of ice surface temperatures)

Delta_T::Delta_T(IceGrid::ConstPtr g, std::shared_ptr<SurfaceModel> in)
  : SurfaceModel(g, in) {

  m_forcing.reset(new ScalarForcing(g->ctx(),
                                    "-surface_delta_T",
                                    "delta_T",
                                    "Kelvin",
                                    "Kelvin",
                                    "ice-surface temperature offsets"));

  m_temperature = allocate_temperature(g);
}

Delta_T::~Delta_T() {
  // empty
}

void Delta_T::init_impl(const Geometry &geometry) {
  m_input_model->init(geometry);

  m_log->message(2,
                 "* Initializing ice-surface temperature forcing using scalar offsets...\n");

  m_forcing->init();
}

void Delta_T::update_impl(const Geometry &geometry, double t, double dt) {
  m_input_model->update(geometry, t, dt);
  m_forcing->update(t, dt);

  m_temperature->copy_from(m_input_model->temperature());
  m_temperature->shift(m_forcing->value());
}

const IceModelVec2S &Delta_T::temperature_impl() const {
  return *m_temperature;
}

} // end of namespace surface
} // end of namespace pism
