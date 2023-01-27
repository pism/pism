// Copyright (C) 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019, 2021, 2022 PISM Authors
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

#include "pism/util/Time.hh"
#include "pism/util/IceGrid.hh"
#include "pism/coupler/util/options.hh"

namespace pism {
namespace surface {

Given::Given(IceGrid::ConstPtr grid, std::shared_ptr<atmosphere::AtmosphereModel> input)
  : SurfaceModel(grid)
{
  (void) input;

  ForcingOptions opt(*m_grid->ctx(), "surface.given");

  {
    unsigned int buffer_size = m_config->get_number("input.forcing.buffer_size");

    File file(m_grid->com, opt.filename, PISM_GUESS, PISM_READONLY);

    m_temperature = std::make_shared<array::Forcing>(m_grid,
                                                file,
                                                "ice_surface_temp",
                                                "", // no standard name
                                                buffer_size,
                                                opt.periodic,
                                                LINEAR);

    m_mass_flux = std::make_shared<array::Forcing>(m_grid,
                                              file,
                                              "climatic_mass_balance",
                                              "land_ice_surface_specific_mass_balance_flux",
                                              buffer_size,
                                              opt.periodic);
  }

  m_temperature->set_attrs("climate_forcing",
                           "temperature of the ice at the ice surface but below firn processes",
                           "Kelvin", "Kelvin", "", 0);
  m_temperature->metadata()["valid_range"] = {0.0, 323.15}; // [0C, 50C]

  const double smb_max = m_config->get_number("surface.given.smb_max", "kg m-2 second-1");

  m_mass_flux->set_attrs("climate_forcing",
                         "surface mass balance (accumulation/ablation) rate",
                         "kg m-2 s-1", "kg m-2 year-1",
                         "land_ice_surface_specific_mass_balance_flux", 0);

  m_mass_flux->metadata()["valid_range"] = {-smb_max, smb_max};
}

void Given::init_impl(const Geometry &geometry) {

  m_log->message(2,
                 "* Initializing the surface model reading temperature at the top of the ice\n"
                 "  and ice surface mass flux from a file...\n");

  ForcingOptions opt(*m_grid->ctx(), "surface.given");

  m_temperature->init(opt.filename, opt.periodic);
  m_mass_flux->init(opt.filename, opt.periodic);

  // read time-independent data right away:
  if (m_temperature->buffer_size() == 1 && m_mass_flux->buffer_size() == 1) {
    update(geometry, m_grid->ctx()->time()->current(), 0); // dt is irrelevant
  }
}

void Given::update_impl(const Geometry &geometry, double t, double dt) {
  (void) geometry;

  m_mass_flux->update(t, dt);
  m_temperature->update(t, dt);

  m_mass_flux->average(t, dt);
  m_temperature->average(t, dt);

  dummy_accumulation(*m_mass_flux, *m_accumulation);
  dummy_melt(*m_mass_flux, *m_melt);
  dummy_runoff(*m_mass_flux, *m_runoff);

}

const array::Scalar &Given::mass_flux_impl() const {
  return *m_mass_flux;
}

const array::Scalar &Given::temperature_impl() const {
  return *m_temperature;
}

const array::Scalar &Given::accumulation_impl() const {
  return *m_accumulation;
}

const array::Scalar &Given::melt_impl() const {
  return *m_melt;
}

const array::Scalar &Given::runoff_impl() const {
  return *m_runoff;
}

void Given::define_model_state_impl(const File &output) const {
  m_mass_flux->define(output);
  m_temperature->define(output);
}

void Given::write_model_state_impl(const File &output) const {
  m_mass_flux->write(output);
  m_temperature->write(output);
}

} // end of namespace surface
} // end of namespace pism
