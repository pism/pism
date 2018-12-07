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
    unsigned int buffer_size = m_config->get_double("climate_forcing.buffer_size");
    unsigned int evaluations_per_year = m_config->get_double("climate_forcing.evaluations_per_year");
    bool periodic = opt.period > 0;

    PIO file(m_grid->com, "netcdf3", opt.filename, PISM_READONLY);


    m_ice_surface_temp_anomaly = IceModelVec2T::ForcingField(m_grid,
                                                             file,
                                                             "ice_surface_temp_anomaly",
                                                             "", // no standard name
                                                             buffer_size,
                                                             evaluations_per_year,
                                                             periodic);

    m_climatic_mass_balance_anomaly = IceModelVec2T::ForcingField(m_grid,
                                                                  file,
                                                                  "climatic_mass_balance_anomaly",
                                                                  "", // no standard name
                                                                  buffer_size,
                                                                  evaluations_per_year,
                                                                  periodic);
  }

  m_ice_surface_temp_anomaly->set_attrs("climate_forcing",
                                        "anomaly of the temperature of the ice at the ice surface"
                                        " but below firn processes",
                                        "Kelvin", "");
  m_climatic_mass_balance_anomaly->set_attrs("climate_forcing",
                                             "anomaly of the surface mass balance (accumulation/ablation) rate",
                                             "kg m-2 s-1", "");
  m_climatic_mass_balance_anomaly->metadata().set_string("glaciological_units", "kg m-2 year-1");

  m_mass_flux = allocate_mass_flux(g);
  m_temperature = allocate_temperature(g);
}

Anomaly::~Anomaly() {
  // empty
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

  m_ice_surface_temp_anomaly->init(opt.filename, opt.period, opt.reference_time);
  m_climatic_mass_balance_anomaly->init(opt.filename, opt.period, opt.reference_time);
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
}

const IceModelVec2S &Anomaly::mass_flux_impl() const {
  return *m_mass_flux;
}

const IceModelVec2S &Anomaly::temperature_impl() const {
  return *m_temperature;
}

} // end of namespace surface
} // end of namespace pism
