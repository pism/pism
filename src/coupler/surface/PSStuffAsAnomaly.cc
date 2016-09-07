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

  mass_flux.create(m_grid, "climatic_mass_balance", WITHOUT_GHOSTS);
  mass_flux.set_attrs("climate_state",
                      "surface mass balance (accumulation/ablation) rate",
                      "kg m-2 s-1",
                      "land_ice_surface_specific_mass_balance_flux");
  mass_flux.metadata().set_string("glaciological_units", "kg m-2 year-1");
  mass_flux.write_in_glaciological_units = true;

  temp.create(m_grid, "ice_surface_temp", WITHOUT_GHOSTS);
  temp.set_attrs("climate_state", "ice temperature at the ice surface",
                 "K", "");

  // create special variables
  mass_flux_0.create(m_grid, "mass_flux_0", WITHOUT_GHOSTS);
  mass_flux_0.set_attrs("internal", "surface mass flux at the beginning of a run",
                        "kg m-2 s-1", "land_ice_surface_specific_mass_balance_flux");

  mass_flux_input.create(m_grid, "climatic_mass_balance", WITHOUT_GHOSTS);
  mass_flux_input.set_attrs("model_state", "surface mass flux to apply anomalies to",
                            "kg m-2 s-1", "land_ice_surface_specific_mass_balance_flux");

  temp_0.create(m_grid, "ice_surface_temp_0", WITHOUT_GHOSTS);
  temp_0.set_attrs("internal", "ice-surface temperature and the beginning of a run", "K",
                   "");

  temp_input.create(m_grid, "ice_surface_temp", WITHOUT_GHOSTS);
  temp_input.set_attrs("model_state", "ice-surface temperature to apply anomalies to",
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
    mass_flux_input.regrid(opts.filename, CRITICAL);
    temp_input.regrid(opts.filename, CRITICAL);
  } else {
    mass_flux_input.read(opts.filename, opts.record);
    temp_input.read(opts.filename, opts.record);
  }
}

MaxTimestep StuffAsAnomaly::max_timestep_impl(double t) {
  (void) t;
  return MaxTimestep();
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
    m_input_model->ice_surface_temperature(temp);
    m_input_model->ice_surface_mass_flux(mass_flux);

    // if we are at the beginning of the run...
    if (m_t < m_grid->ctx()->time()->start() + 1) { // this is goofy, but time-steps are
                                      // usually longer than 1 second, so it
                                      // should work
      temp_0.copy_from(temp);
      mass_flux_0.copy_from(mass_flux);
    }
  }

  IceModelVec::AccessList list;
  list.add(mass_flux);
  list.add(mass_flux_0);
  list.add(mass_flux_input);

  list.add(temp);
  list.add(temp_0);
  list.add(temp_input);

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    mass_flux(i, j) = mass_flux(i, j) - mass_flux_0(i, j) + mass_flux_input(i, j);
    temp(i, j)      = temp(i, j) - temp_0(i, j) + temp_input(i, j);
  }
}

void StuffAsAnomaly::ice_surface_mass_flux_impl(IceModelVec2S &result) {
  result.copy_from(mass_flux);
}

void StuffAsAnomaly::ice_surface_temperature_impl(IceModelVec2S &result) {
  result.copy_from(temp);
}

void StuffAsAnomaly::add_vars_to_output_impl(const std::string &keyword, std::set<std::string> &result) {
  if (m_input_model != NULL) {
    m_input_model->add_vars_to_output(keyword, result);
  }

  result.insert("ice_surface_temp");
  result.insert("climatic_mass_balance");
}

void StuffAsAnomaly::define_variables_impl(const std::set<std::string> &vars_input,
                                                  const PIO &nc, IO_Type nctype) {
  std::set<std::string> vars = vars_input;

  if (set_contains(vars, "ice_surface_temp")) {
    temp.define(nc, nctype);
  }

  if (set_contains(vars, "climatic_mass_balance")) {
    mass_flux.define(nc, nctype);
  }

  // ensure that no one overwrites these two
  vars.erase("ice_surface_temp");
  vars.erase("climatic_mass_balance");

  if (m_input_model != NULL) {
    m_input_model->define_variables(vars, nc, nctype);
  }
}

void StuffAsAnomaly::write_variables_impl(const std::set<std::string> &vars_input, const PIO &nc) {
  std::set<std::string> vars = vars_input;

  if (set_contains(vars, "ice_surface_temp")) {
    temp.write(nc);
  }

  if (set_contains(vars, "climatic_mass_balance")) {
    mass_flux.write(nc);
  }

  // ensure that no one overwrites these two
  vars.erase("ice_surface_temp");
  vars.erase("climatic_mass_balance");

  if (m_input_model != NULL) {
    m_input_model->write_variables(vars, nc);
  }
}


} // end of namespace surface
} // end of namespace pism
