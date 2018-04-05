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

#include "Delta_SMB.hh"
#include "pism/coupler/util/ScalarForcing.hh"

namespace pism {
namespace ocean {

Delta_SMB::Delta_SMB(IceGrid::ConstPtr g, std::shared_ptr<OceanModel> in)
  : OceanModel(g, in) {

  m_forcing.reset(new ScalarForcing(g->ctx(),
                                    "-ocean_delta_mass_flux",
                                    "delta_mass_flux",
                                    "kg m-2 second-1",
                                    "kg m-2 year-1",
                                    "ice-shelf-base mass flux offsets"));

  m_shelf_base_mass_flux = allocate_shelf_base_mass_flux(g);
}

Delta_SMB::~Delta_SMB() {
  // empty
}

void Delta_SMB::init_impl(const Geometry &geometry) {
  m_input_model->init(geometry);

  m_log->message(2,
                 "* Initializing ice shelf base mass flux forcing using scalar offsets...\n");

  m_forcing->init();
}

void Delta_SMB::update_impl(const Geometry &geometry, double t, double dt) {
  m_input_model->update(geometry, t, dt);

  m_forcing->update(t, dt);

  m_shelf_base_mass_flux->copy_from(m_input_model->shelf_base_mass_flux());
  m_shelf_base_mass_flux->shift(m_forcing->value());
}

const IceModelVec2S& Delta_SMB::shelf_base_mass_flux_impl() const {
  return *m_shelf_base_mass_flux;
}

} // end of namespace ocean
} // end of namespace pism
