/* Copyright (C) 2013, 2014 PISM Authors
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

#include "PO_delta_MBP.hh"
#include "PISMConfig.hh"

namespace pism {

PO_delta_MBP::PO_delta_MBP(const IceGrid &g, OceanModel* in)
  : PScalarForcing<OceanModel,POModifier>(g, in),
    shelfbmassflux(g.config.get_unit_system(), "shelfbmassflux", m_grid),
    shelfbtemp(g.config.get_unit_system(), "shelfbtemp", m_grid)
{

  option_prefix = "-ocean_delta_MBP";
  offset_name   = "delta_MBP";

  offset = new Timeseries(&m_grid, offset_name, m_config.get_string("time_dimension_name"));

  offset->get_metadata().set_units("1");
  offset->get_metadata().set_string("long_name", "melange back pressure fraction");
  offset->get_dimension_metadata().set_units(m_grid.time->units_string());

  shelfbmassflux.set_string("pism_intent", "climate_state");
  shelfbmassflux.set_string("long_name",
                            "ice mass flux from ice shelf base (positive flux is loss from ice shelf)");
  shelfbmassflux.set_units("kg m-2 s-1");
  shelfbmassflux.set_glaciological_units("kg m-2 year-1");

  shelfbtemp.set_string("pism_intent", "climate_state");
  shelfbtemp.set_string("long_name",
                        "absolute temperature at ice shelf base");
  shelfbtemp.set_units("Kelvin");
}

PO_delta_MBP::~PO_delta_MBP() {
  // empty
}

void PO_delta_MBP::init() {
  m_t = m_dt = GSL_NAN;  // every re-init restarts the clock

  input_model->init();

  verbPrintf(2, m_grid.com, "* Initializing melange back pressure fraction forcing...\n");

  init_internal();
}

void PO_delta_MBP::melange_back_pressure_fraction(IceModelVec2S &result) {
  input_model->melange_back_pressure_fraction(result);

  result.shift((*offset)(m_t + 0.5*m_dt));
}

void PO_delta_MBP::add_vars_to_output(const std::string &keyword, std::set<std::string> &result) {
  input_model->add_vars_to_output(keyword, result);

  result.insert("shelfbtemp");
  result.insert("shelfbmassflux");
}

void PO_delta_MBP::define_variables(const std::set<std::string> &vars_input, const PIO &nc,
                                              IO_Type nctype) {
  std::set<std::string> vars = vars_input;

  if (set_contains(vars, "shelfbtemp")) {
    shelfbtemp.define(nc, nctype, true);
    vars.erase("shelfbtemp");
  }

  if (set_contains(vars, "shelfbmassflux")) {
    shelfbmassflux.define(nc, nctype, true);
    vars.erase("shelfbmassflux");
  }

  input_model->define_variables(vars, nc, nctype);
}

void PO_delta_MBP::write_variables(const std::set<std::string> &vars_input, const PIO &nc) {
  std::set<std::string> vars = vars_input;
  IceModelVec2S tmp;

  if (set_contains(vars, "shelfbtemp")) {
    if (!tmp.was_created()) {
      tmp.create(m_grid, "tmp", WITHOUT_GHOSTS);
    }

    tmp.metadata() = shelfbtemp;
    shelf_base_temperature(tmp);
    tmp.write(nc);
    vars.erase("shelfbtemp");
  }

  if (set_contains(vars, "shelfbmassflux")) {
    if (!tmp.was_created()) {
      tmp.create(m_grid, "tmp", WITHOUT_GHOSTS);
    }

    tmp.metadata() = shelfbmassflux;
    tmp.write_in_glaciological_units = true;
    shelf_base_mass_flux(tmp);
    tmp.write(nc);
    vars.erase("shelfbmassflux");
  }

  input_model->write_variables(vars, nc);
}


} // end of namespace pism
