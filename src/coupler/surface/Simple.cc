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


#include "Simple.hh"

#include "pism/coupler/AtmosphereModel.hh"

namespace pism {
namespace surface {

Simple::Simple(IceGrid::ConstPtr g, std::shared_ptr<atmosphere::AtmosphereModel> atmosphere)
  : SurfaceModel(g, atmosphere) {

  m_temperature = allocate_temperature(g);
  m_mass_flux   = allocate_mass_flux(g);
}

void Simple::init_impl(const Geometry &geometry) {

  m_atmosphere->init(geometry);

  m_log->message(2,
             "* Initializing the simplest PISM surface (snow) processes model Simple.\n"
             "  It passes atmospheric state directly to upper ice fluid surface:\n"
             "    surface mass balance          := precipitation,\n"
             "    ice upper surface temperature := 2m air temperature.\n");
}

void Simple::update_impl(const Geometry &geometry, double t, double dt) {
  if (m_atmosphere) {
    m_atmosphere->update(geometry, t, dt);
  }

  m_mass_flux->copy_from(m_atmosphere->mean_precipitation());
  m_temperature->copy_from(m_atmosphere->mean_annual_temp());
}

const IceModelVec2S &Simple::mass_flux_impl() const {
  return *m_mass_flux;
}

const IceModelVec2S &Simple::temperature_impl() const {
  return *m_temperature;
}

} // end of namespace surface
} // end of namespace pism
