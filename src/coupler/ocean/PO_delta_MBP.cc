/* Copyright (C) 2013, 2014, 2015, 2016 PISM Authors
 *
 * This file is part of PISM.
 *
 * PISM is free software; you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation; either version 3 of the License, or (at your option) any later
 * version.
 *
 * PISM is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License
 * along with PISM; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 */

#include <gsl/gsl_math.h>

#include "PO_delta_MBP.hh"
#include "base/util/PISMConfigInterface.hh"
#include "base/util/io/io_helpers.hh"
#include "base/util/pism_utilities.hh"

namespace pism {
namespace ocean {

Delta_MBP::Delta_MBP(IceGrid::ConstPtr g, OceanModel* in)
  : PScalarForcing<OceanModel,OceanModifier>(g, in),
    m_shelfbmassflux(m_sys, "effective_shelf_base_mass_flux"),
    m_shelfbtemp(m_sys, "effective_shelf_base_temperature")
{

  m_option_prefix = "-ocean_delta_MBP";
  m_offset_name   = "delta_MBP";

  m_offset = new Timeseries(*m_grid, m_offset_name, m_config->get_string("time_dimension_name"));

  m_offset->metadata().set_string("units", "1");
  m_offset->metadata().set_string("long_name", "melange back pressure fraction");
  m_offset->dimension_metadata().set_string("units", m_grid->ctx()->time()->units_string());

  m_shelfbmassflux.set_string("pism_intent", "climate_state");
  m_shelfbmassflux.set_string("long_name",
                            "ice mass flux from ice shelf base (positive flux is loss from ice shelf)");
  m_shelfbmassflux.set_string("units", "kg m-2 s-1");
  m_shelfbmassflux.set_string("glaciological_units", "kg m-2 year-1");

  m_shelfbtemp.set_string("pism_intent", "climate_state");
  m_shelfbtemp.set_string("long_name",
                        "absolute temperature at ice shelf base");
  m_shelfbtemp.set_string("units", "Kelvin");
}

Delta_MBP::~Delta_MBP() {
  // empty
}

void Delta_MBP::init_impl() {
  m_t = m_dt = GSL_NAN;  // every re-init restarts the clock

  m_input_model->init();

  m_log->message(2, "* Initializing melange back pressure fraction forcing...\n");

  init_internal();
}

MaxTimestep Delta_MBP::max_timestep_impl(double t) {
  (void) t;
  return MaxTimestep();
}

void Delta_MBP::melange_back_pressure_fraction_impl(IceModelVec2S &result) {
  m_input_model->melange_back_pressure_fraction(result);

  offset_data(result);
}

void Delta_MBP::add_vars_to_output_impl(const std::string &keyword, std::set<std::string> &result) {
  m_input_model->add_vars_to_output(keyword, result);

  result.insert(m_shelfbtemp.get_name());
  result.insert(m_shelfbmassflux.get_name());
}

void Delta_MBP::define_variables_impl(const std::set<std::string> &vars_input, const PIO &nc,
                                              IO_Type nctype) {
  std::set<std::string> vars = vars_input;
  std::string order = m_config->get_string("output_variable_order");

  if (set_contains(vars, m_shelfbtemp)) {
    io::define_spatial_variable(m_shelfbtemp, *m_grid, nc, nctype, order, true);
    vars.erase(m_shelfbtemp.get_name());
  }

  if (set_contains(vars, m_shelfbmassflux)) {
    io::define_spatial_variable(m_shelfbmassflux, *m_grid, nc, nctype, order, true);
    vars.erase(m_shelfbmassflux.get_name());
  }

  m_input_model->define_variables(vars, nc, nctype);
}

void Delta_MBP::write_variables_impl(const std::set<std::string> &vars_input, const PIO &nc) {
  std::set<std::string> vars = vars_input;
  IceModelVec2S tmp;

  if (set_contains(vars, m_shelfbtemp)) {
    if (!tmp.was_created()) {
      tmp.create(m_grid, "tmp", WITHOUT_GHOSTS);
    }

    tmp.metadata() = m_shelfbtemp;
    shelf_base_temperature(tmp);
    tmp.write(nc);
    vars.erase(m_shelfbtemp.get_name());
  }

  if (set_contains(vars, m_shelfbmassflux)) {
    if (!tmp.was_created()) {
      tmp.create(m_grid, "tmp", WITHOUT_GHOSTS);
    }

    tmp.metadata() = m_shelfbmassflux;
    tmp.write_in_glaciological_units = true;
    shelf_base_mass_flux(tmp);
    tmp.write(nc);
    vars.erase(m_shelfbmassflux.get_name());
  }

  m_input_model->write_variables(vars, nc);
}


} // end of namespace ocean
} // end of namespace pism
