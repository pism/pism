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

#include "Anomaly.hh"
#include "pism/util/IceGrid.hh"
#include "pism/util/io/io_helpers.hh"
#include "pism/coupler/util/options.hh"

namespace pism {
namespace surface {

Anomaly::Anomaly(IceGrid::ConstPtr g, std::shared_ptr<SurfaceModel> in)
  : SurfaceModel(g, in) {

  ForcingOptions opt(*m_grid->ctx(), "surface.anomaly");

  {
    unsigned int buffer_size = m_config->get_number("input.forcing.buffer_size");

    File file(m_grid->com, opt.filename, PISM_NETCDF3, PISM_READONLY);

    m_ice_surface_temp_anomaly = std::make_shared<IceModelVec2T>(m_grid,
                                                             file,
                                                             "ice_surface_temp_anomaly",
                                                             "", // no standard name
                                                             buffer_size,
                                                             opt.periodic,
                                                             LINEAR);

    m_climatic_mass_balance_anomaly = std::make_shared<IceModelVec2T>(m_grid,
                                                                  file,
                                                                  "climatic_mass_balance_anomaly",
                                                                  "", // no standard name
                                                                  buffer_size,
                                                                  opt.periodic);
  }

  m_ice_surface_temp_anomaly->set_attrs("climate_forcing",
                                        "anomaly of the temperature of the ice at the ice surface"
                                        " but below firn processes",
                                        "Kelvin", "Kelvin", "", 0);
  m_climatic_mass_balance_anomaly->set_attrs("climate_forcing",
                                             "anomaly of the surface mass balance (accumulation/ablation) rate",
                                             "kg m-2 s-1", "kg m-2 year-1", "", 0);

  m_mass_flux = allocate_mass_flux(g);
  m_temperature = allocate_temperature(g);

  m_accumulation = allocate_accumulation(g);
  m_melt         = allocate_melt(g);
  m_runoff       = allocate_runoff(g);
}

void Anomaly::init_impl(const Geometry &geometry) {

  if (m_input_model) {
    m_input_model->init(geometry);
  }

  m_log->message(2,
                 "* Initializing the '-surface ...,anomaly' modifier...\n");

  ForcingOptions opt(*m_grid->ctx(), "surface.anomaly");

  m_log->message(2,
                 "    reading anomalies from %s ...\n", opt.filename.c_str());

  m_ice_surface_temp_anomaly->init(opt.filename, opt.periodic);
  m_climatic_mass_balance_anomaly->init(opt.filename, opt.periodic);
}

void Anomaly::update_impl(const Geometry &geometry, double t, double dt) {
  m_input_model->update(geometry, t, dt);

  m_climatic_mass_balance_anomaly->update(t, dt);
  m_ice_surface_temp_anomaly->update(t, dt);

  m_climatic_mass_balance_anomaly->average(t, dt);
  m_ice_surface_temp_anomaly->average(t, dt);

  m_input_model->mass_flux().add(1.0, *m_climatic_mass_balance_anomaly,
                                 *m_mass_flux);
  m_input_model->temperature().add(1.0, *m_ice_surface_temp_anomaly,
                                   *m_temperature);

  dummy_accumulation(*m_mass_flux, *m_accumulation);
  dummy_melt(*m_mass_flux, *m_melt);
  dummy_runoff(*m_mass_flux, *m_runoff);
}

const array::Scalar &Anomaly::mass_flux_impl() const {
  return *m_mass_flux;
}

const array::Scalar &Anomaly::temperature_impl() const {
  return *m_temperature;
}

const array::Scalar &Anomaly::accumulation_impl() const {
  return *m_accumulation;
}

const array::Scalar &Anomaly::melt_impl() const {
  return *m_melt;
}

const array::Scalar &Anomaly::runoff_impl() const {
  return *m_runoff;
}

} // end of namespace surface
} // end of namespace pism
