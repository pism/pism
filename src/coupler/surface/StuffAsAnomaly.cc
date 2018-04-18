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

#include "StuffAsAnomaly.hh"
#include "pism/util/IceGrid.hh"
#include "pism/util/Time.hh"
#include "pism/util/pism_utilities.hh"
#include "pism/util/MaxTimestep.hh"

namespace pism {
namespace surface {

StuffAsAnomaly::StuffAsAnomaly(IceGrid::ConstPtr g, std::shared_ptr<SurfaceModel> input)
    : SurfaceModel(g, input) {

  m_mass_flux.create(m_grid, "climatic_mass_balance", WITHOUT_GHOSTS);
  m_mass_flux.set_attrs("climate_state",
                      "surface mass balance (accumulation/ablation) rate",
                      "kg m-2 s-1",
                      "land_ice_surface_specific_mass_balance_flux");
  m_mass_flux.metadata().set_string("glaciological_units", "kg m-2 year-1");

  m_temp.create(m_grid, "ice_surface_temp", WITHOUT_GHOSTS);
  m_temp.set_attrs("climate_state", "ice temperature at the ice surface",
                 "K", "");

  // create special variables
  m_mass_flux_0.create(m_grid, "mass_flux_0", WITHOUT_GHOSTS);
  m_mass_flux_0.set_attrs("internal", "surface mass flux at the beginning of a run",
                        "kg m-2 s-1", "land_ice_surface_specific_mass_balance_flux");

  m_mass_flux_input.create(m_grid, "climatic_mass_balance", WITHOUT_GHOSTS);
  m_mass_flux_input.set_attrs("model_state", "surface mass flux to apply anomalies to",
                            "kg m-2 s-1", "land_ice_surface_specific_mass_balance_flux");

  m_temp_0.create(m_grid, "ice_surface_temp_0", WITHOUT_GHOSTS);
  m_temp_0.set_attrs("internal", "ice-surface temperature and the beginning of a run", "K",
                   "");

  m_temp_input.create(m_grid, "ice_surface_temp", WITHOUT_GHOSTS);
  m_temp_input.set_attrs("model_state", "ice-surface temperature to apply anomalies to",
                       "K", "");
}

StuffAsAnomaly::~StuffAsAnomaly() {
  // empty
}

void StuffAsAnomaly::init_impl(const Geometry &geometry) {
  if (m_input_model != NULL) {
    m_input_model->init(geometry);
  }

  InputOptions opts = process_input_options(m_grid->com, m_config);

  m_log->message(2,
             "* Initializing the 'turn_into_anomaly' modifier\n"
             "  (it applies climate data as anomalies relative to 'ice_surface_temp' and 'climatic_mass_balance'\n"
             "  read from '%s'.\n", opts.filename.c_str());

  if (opts.type == INIT_BOOTSTRAP) {
    m_mass_flux_input.regrid(opts.filename, CRITICAL);
    m_temp_input.regrid(opts.filename, CRITICAL);
  } else {
    m_mass_flux_input.read(opts.filename, opts.record);
    m_temp_input.read(opts.filename, opts.record);
  }
}

MaxTimestep StuffAsAnomaly::max_timestep_impl(double t) const {
  (void) t;
  return MaxTimestep("surface turn_into_anomaly");
}

void StuffAsAnomaly::update_impl(const Geometry &geometry, double t, double dt) {

  if (m_input_model != NULL) {
    m_input_model->update(geometry, t, dt);
    m_temp.copy_from(m_input_model->temperature());
    m_mass_flux.copy_from(m_input_model->mass_flux());

    // if we are at the beginning of the run...
    if (t < m_grid->ctx()->time()->start() + 1) {
      // this is goofy, but time-steps are usually longer than 1 second, so it should work
      m_temp_0.copy_from(m_temp);
      m_mass_flux_0.copy_from(m_mass_flux);
    }
  }

  IceModelVec::AccessList list{&m_mass_flux, &m_mass_flux_0, &m_mass_flux_input,
      &m_temp, &m_temp_0, &m_temp_input};

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    m_mass_flux(i, j) = m_mass_flux(i, j) - m_mass_flux_0(i, j) + m_mass_flux_input(i, j);
    m_temp(i, j)      = m_temp(i, j) - m_temp_0(i, j) + m_temp_input(i, j);
  }
}

const IceModelVec2S &StuffAsAnomaly::mass_flux_impl() const {
  return m_mass_flux;
}

const IceModelVec2S &StuffAsAnomaly::temperature_impl() const {
  return m_temp;
}

} // end of namespace surface
} // end of namespace pism
