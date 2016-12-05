// Copyright (C) 2004--2016 Jed Brown, Ed Bueler and Constantine Khroulev
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

#include <cstring>
#include <cmath>
#include <sstream>
#include <set>

#include "iceModel.hh"

#include "base/basalstrength/PISMYieldStress.hh"
#include "base/calving/OceanKill.hh"
#include "base/energy/BedThermalUnit.hh"
#include "base/hydrology/PISMHydrology.hh"
#include "base/stressbalance/PISMStressBalance.hh"
#include "base/util/IceGrid.hh"
#include "base/util/PISMConfigInterface.hh"
#include "base/util/error_handling.hh"
#include "base/util/pism_options.hh"
#include "coupler/PISMOcean.hh"
#include "coupler/PISMSurface.hh"
#include "earth/PISMBedDef.hh"
#include "base/util/PISMVars.hh"
#include "base/util/pism_utilities.hh"
#include "base/age/AgeModel.hh"

namespace pism {

//! Read some runtime (command line) options and alter the corresponding parameters or flags as appropriate.
void IceModel::setFromOptions() {

  m_log->message(3,
             "Processing physics-related command-line options...\n");

  set_config_from_options(*m_config);

  // Set global attributes using the config database:
  m_output_global_attributes.set_string("title", m_config->get_string("run_info.title"));
  m_output_global_attributes.set_string("institution", m_config->get_string("run_info.institution"));
  m_output_global_attributes.set_string("command", pism_args_string());

  // warn about some option combinations

  if (m_config->get_double("time_stepping.maximum_time_step") <= 0) {
    throw RuntimeError(PISM_ERROR_LOCATION, "time_stepping.maximum_time_step has to be greater than 0.");
  }

  if (not m_config->get_boolean("geometry.update.enabled") &&
      m_config->get_boolean("time_stepping.skip.enabled")) {
    m_log->message(2,
               "PISM WARNING: Both -skip and -no_mass are set.\n"
               "              -skip only makes sense in runs updating ice geometry.\n");
  }

  if (m_config->get_string("calving.methods").find("thickness_calving") != std::string::npos &&
      not m_config->get_boolean("geometry.part_grid.enabled")) {
    m_log->message(2,
               "PISM WARNING: Calving at certain terminal ice thickness (-calving thickness_calving)\n"
               "              without application of partially filled grid cell scheme (-part_grid)\n"
               "              may lead to (incorrect) non-moving ice shelf front.\n");
  }
}

//! Set the output file size using a command-line option.
std::set<std::string> IceModel::output_size_from_option(const std::string &option,
                                                        const std::string &description,
                                                        const std::string &default_value) {

  options::Keyword o_size(option, description, "none,small,medium,big,big_2d",
                          default_value);

  return set_output_size(o_size);
}

//! Assembles a list of variables corresponding to an output file size.
std::set<std::string> IceModel::set_output_size(const std::string &keyword) {
  std::set<std::string> result;

  if (keyword == "none") {
    return result;
  }

  // Add all the model-state variables:
  for (auto v : m_grid->variables().keys()) {
    const SpatialVariableMetadata &m = m_grid->variables().get(v)->metadata();

    std::string intent = m.get_string("pism_intent");

    if (intent == "model_state" or
        intent == "mapping"     or
        intent == "climate_steady") {
      result.insert(v);
    }
  }

  // add cumulative quantities to ensure continuity after restarting
  result.insert("climatic_mass_balance_cumulative");
  result.insert("grounded_basal_flux_cumulative");
  result.insert("floating_basal_flux_cumulative");
  result.insert("nonneg_flux_cumulative");
  result.insert("discharge_flux_cumulative");

  std::string variables;
  if (keyword == "medium") {
    // add all the variables listed in the config file ("medium" size):
    variables = m_config->get_string("output.sizes.medium");
  } else if (keyword == "big_2d") {
    // add all the variables listed in the config file (under "medium" and "big_2d" sizes):
    variables = m_config->get_string("output.sizes.medium");
    variables += "," + m_config->get_string("output.sizes.big_2d");
  } else if (keyword == "big") {
    // add all the variables listed in the config file ("big" size):
    variables = m_config->get_string("output.sizes.medium");
    variables += "," + m_config->get_string("output.sizes.big_2d");
    variables += "," + m_config->get_string("output.sizes.big");
  }

  for (auto name : split(variables, ',')) {
    if (not name.empty()) {
      result.insert(name);
    }
  }

  return result;
}


//! Returns the output size as a keyword, for options "-o_size", "-save_size", "-backup_size", etc.
std::string IceModel::get_output_size(const std::string &option) {
  return options::Keyword(option, "no description", "none,small,medium,big,big_2d", "no default");
}


} // end of namespace pism
