// Copyright (C) 2011, 2012, 2013, 2014, 2015, 2016 PISM Authors
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

#include <gsl/gsl_math.h>

#include "PSStuffAsAnomaly.hh"
#include "base/util/IceGrid.hh"
#include "base/util/PISMTime.hh"
#include "base/util/pism_utilities.hh"

namespace pism {
namespace surface {

StuffAsAnomaly::StuffAsAnomaly(IceGrid::ConstPtr g, SurfaceModel *input)
    : SurfaceModifier(g, input) {

  m_mass_flux.create(m_grid, "climatic_mass_balance", WITHOUT_GHOSTS);
  m_mass_flux.set_attrs("climate_state",
                      "surface mass balance (accumulation/ablation) rate",
                      "kg m-2 s-1",
                      "land_ice_surface_specific_mass_balance_flux");
  m_mass_flux.metadata().set_string("glaciological_units", "kg m-2 year-1");
  m_mass_flux.write_in_glaciological_units = true;

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

void StuffAsAnomaly::init_impl() {
  m_t = m_dt = GSL_NAN;  // every re-init restarts the clock

  if (m_input_model != NULL) {
    m_input_model->init();
  }

  InputOptions opts = process_input_options(m_grid->com);

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

void StuffAsAnomaly::update_impl(double my_t, double my_dt) {
  if ((fabs(my_t - m_t) < 1e-12) &&
      (fabs(my_dt - m_dt) < 1e-12)) {
    return;
  }

  m_t  = my_t;
  m_dt = my_dt;

  if (m_input_model != NULL) {
    m_input_model->update(m_t, m_dt);
    m_input_model->ice_surface_temperature(m_temp);
    m_input_model->ice_surface_mass_flux(m_mass_flux);

    // if we are at the beginning of the run...
    if (m_t < m_grid->ctx()->time()->start() + 1) { // this is goofy, but time-steps are
                                      // usually longer than 1 second, so it
                                      // should work
      m_temp_0.copy_from(m_temp);
      m_mass_flux_0.copy_from(m_mass_flux);
    }
  }

  IceModelVec::AccessList list;
  list.add(m_mass_flux);
  list.add(m_mass_flux_0);
  list.add(m_mass_flux_input);

  list.add(m_temp);
  list.add(m_temp_0);
  list.add(m_temp_input);

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    m_mass_flux(i, j) = m_mass_flux(i, j) - m_mass_flux_0(i, j) + m_mass_flux_input(i, j);
    m_temp(i, j)      = m_temp(i, j) - m_temp_0(i, j) + m_temp_input(i, j);
  }
}

void StuffAsAnomaly::ice_surface_mass_flux_impl(IceModelVec2S &result) const {
  result.copy_from(m_mass_flux);
}

void StuffAsAnomaly::ice_surface_temperature_impl(IceModelVec2S &result) const {
  result.copy_from(m_temp);
}

} // end of namespace surface
} // end of namespace pism
