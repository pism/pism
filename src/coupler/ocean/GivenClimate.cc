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

#include "GivenClimate.hh"
#include "pism/util/IceGrid.hh"

namespace pism {
namespace ocean {

Given::Given(IceGrid::ConstPtr g)
  : PGivenClimate<OceanModel>(g, nullptr) {

  m_sea_level_elevation    = allocate_sea_level_elevation(g);
  m_shelf_base_temperature = allocate_shelf_base_temperature(g);
  m_shelf_base_mass_flux   = allocate_shelf_base_mass_flux(g);

  m_filename = process_options("-ocean_given");

  {
    unsigned int buffer_size = m_config->get_double("climate_forcing.buffer_size");
    PIO file(m_grid->com, "netcdf3", m_filename, PISM_READONLY);

    for (auto name : {"shelfbtemp", "shelfbmassflux"}) {
      m_fields[name] = allocate(file,
                                m_sys,
                                name,
                                "", // no standard name
                                buffer_size,
                                m_bc_period > 0);

    }
  }

  m_shelfbtemp     = m_fields["shelfbtemp"].get();
  m_shelfbmassflux = m_fields["shelfbmassflux"].get();

  m_shelfbtemp->set_attrs("climate_forcing",
                        "absolute temperature at ice shelf base",
                        "Kelvin", "");
  m_shelfbmassflux->set_attrs("climate_forcing",
                            "ice mass flux from ice shelf base (positive flux is loss from ice shelf)",
                            "kg m-2 s-1", "");
  m_shelfbmassflux->metadata().set_string("glaciological_units", "kg m-2 year-1");
}

Given::~Given() {
  // empty
}

void Given::init_impl(const Geometry &geometry) {

  m_log->message(2,
             "* Initializing the ocean model reading base of the shelf temperature\n"
             "  and sub-shelf mass flux from a file...\n");

  m_shelfbtemp->init(m_filename, m_bc_period, m_bc_reference_time);
  m_shelfbmassflux->init(m_filename, m_bc_period, m_bc_reference_time);

  // read time-independent data right away:
  if (m_shelfbtemp->n_records() == 1 && m_shelfbmassflux->n_records() == 1) {
    update(geometry, m_grid->ctx()->time()->current(), 0); // dt is irrelevant
  }
}

void Given::update_impl(const Geometry &geometry, double t, double dt) {
  update_internal(geometry, t, dt);

  m_shelfbmassflux->average(t, dt);
  m_shelfbtemp->average(t, dt);

  m_shelf_base_temperature->copy_from(*m_shelfbtemp);
  m_shelf_base_mass_flux->copy_from(*m_shelfbmassflux);
}

const IceModelVec2S& Given::shelf_base_temperature_impl() const {
  return *m_shelf_base_temperature;
}

const IceModelVec2S& Given::shelf_base_mass_flux_impl() const {
  return *m_shelf_base_mass_flux;
}

const IceModelVec2S& Given::sea_level_elevation_impl() const {
  return *m_sea_level_elevation;
}

} // end of namespace ocean
} // end of namespace pism
