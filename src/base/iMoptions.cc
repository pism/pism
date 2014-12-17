// Copyright (C) 2004--2014 Jed Brown, Ed Bueler and Constantine Khroulev
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
#include "PISMBedDef.hh"
#include "bedrockThermalUnit.hh"
#include "PISMYieldStress.hh"
#include "PISMOceanKill.hh"
#include "PISMHydrology.hh"
#include "PISMStressBalance.hh"
#include "PISMOcean.hh"
#include "PISMSurface.hh"
#include "pism_options.hh"
#include "IceGrid.hh"

#include "error_handling.hh"

namespace pism {

//! Read some runtime (command line) options and alter the corresponding parameters or flags as appropriate.
void  IceModel::setFromOptions() {
  bool flag;

  verbPrintf(3, grid.com,
             "Processing physics-related command-line options...\n");

  set_config_from_options(config);

  OptionsInt("-id", "Specifies the sounding row", id, flag);
  OptionsInt("-jd", "Specifies the sounding column", jd, flag);

  // Set global attributes using the config database:
  global_attributes.set_string("title", config.get_string("run_title"));
  global_attributes.set_string("institution", config.get_string("institution"));
  global_attributes.set_string("command", pism_args_string());

  // warn about some option combinations

  if (config.get("maximum_time_step_years") <= 0) {
    throw RuntimeError("maximum_time_step_years has to be greater than 0.");
  }
  
  if (config.get_flag("do_mass_conserve") == false &&
      config.get_flag("do_skip")) {
    verbPrintf(2, grid.com,
               "PISM WARNING: Both -skip and -no_mass are set.\n"
               "              -skip only makes sense in runs updating ice geometry.\n");
  }

  if (config.get_string("calving_methods").find("thickness_calving") != std::string::npos &&
      config.get_flag("part_grid") == false) {
    verbPrintf(2, grid.com,
               "PISM WARNING: Calving at certain terminal ice thickness (-calving thickness_calving)\n"
               "              without application of partially filled grid cell scheme (-part_grid)\n"
               "              may lead to (incorrect) non-moving ice shelf front.\n");
  }


  // implements an option e.g. described in \ref Greve that is the
  // enhancement factor is coupled to the age of the ice with
  // e = 1 (A < 11'000 years), e = 3 otherwise
  if (config.get_flag("e_age_coupling")) {
    verbPrintf(2, grid.com,
               "  setting age-dependent enhancement factor: "
               "e=1 if A<11'000 years, e=3 otherwise\n");

  }
}

//! Set the output file size using a command-line option.
void IceModel::output_size_from_option(const std::string &option,
                                                 const std::string &description,
                                                 const std::string &default_value,
                                                 std::set<std::string> &result) {

  std::set<std::string> choices;
  std::string keyword;
  bool flag;

  choices.insert("none");
  choices.insert("small");
  choices.insert("medium");
  choices.insert("big");
  OptionsList(option, description, choices,
              default_value, keyword, flag);

  set_output_size(keyword, result);
}

//! Assembles a list of variables corresponding to an output file size.
void IceModel::set_output_size(const std::string &keyword,
                                         std::set<std::string> &result) {
  result.clear();

  if (keyword == "none") {
    return;
  }

  // Add all the model-state variables:
  std::set<std::string> vars = grid.variables().keys();

  std::set<std::string>::const_iterator i = vars.begin();
  while (i != vars.end()) {
    IceModelVec *var = grid.variables().get(*i);
    NCSpatialVariable &m = var->metadata();

    std::string intent = m.get_string("pism_intent");
    if ((intent == "model_state") || (intent == "mapping") || (intent == "climate_steady")) {
      result.insert(*i);
    }
    ++i;
  }

  // add cumulative quantities to ensure continuity after restarting
  if (climatic_mass_balance_cumulative.was_created()) {
    result.insert("climatic_mass_balance_cumulative");
  }
  if (grounded_basal_flux_2D_cumulative.was_created()) {
    result.insert("grounded_basal_flux_2D_cumulative");
  }
  if (floating_basal_flux_2D_cumulative.was_created()) {
    result.insert("floating_basal_flux_2D_cumulative");
  }
  if (nonneg_flux_2D_cumulative.was_created()) {
    result.insert("nonneg_flux_2D_cumulative");
  }
  if (discharge_flux_2D_cumulative.was_created()) {
    result.insert("discharge_flux_cumulative");
  }

  if (keyword == "medium") {
    // add all the variables listed in the config file ("medium" size):
    std::string tmp = config.get_string("output_medium");
    std::istringstream keywords(tmp);

    // split the list; note that this also removes any duplicate entries
    while (getline(keywords, tmp, ' ')) {
      if (!tmp.empty()) {                // this ignores multiple spaces separating variable names
       result.insert(tmp);
      }
    }
  } else if (keyword == "big") {
    // add all the variables listed in the config file ("big" size):
    std::string tmp = config.get_string("output_big");
    std::istringstream keywords(tmp);

    // split the list; note that this also removes any duplicate entries
    while (getline(keywords, tmp, ' ')) {
      if (!tmp.empty()) {                // this ignores multiple spaces separating variable names
       result.insert(tmp);
      }
    }

    if (!config.get_flag("do_age")) {
      result.erase("age");
    }
  }

  if (config.get_flag("do_age")) {
    result.insert("age");
  }

  if (ocean_kill_calving != NULL) {
    ocean_kill_calving->add_vars_to_output(keyword, result);
  }

  if (beddef != NULL) {
    beddef->add_vars_to_output(keyword, result);
  }

  if (btu != NULL) {
    btu->add_vars_to_output(keyword, result);
  }

  if (basal_yield_stress_model != NULL) {
    basal_yield_stress_model->add_vars_to_output(keyword, result);
  }

  // Ask the stress balance module to add more variables:
  if (stress_balance != NULL) {
    stress_balance->add_vars_to_output(keyword, result);
  }

  if (subglacial_hydrology != NULL) {
    subglacial_hydrology->add_vars_to_output(keyword, result);
  }

  // Ask ocean and surface models to add more variables to the list:
  if (ocean != NULL) {
    ocean->add_vars_to_output(keyword, result);
  }

  if (surface != NULL) {
    surface->add_vars_to_output(keyword, result);
  }
}


//! Returns the output size as a keyword, for options "-o_size", "-save_size", "-backup_size", etc.
std::string IceModel::get_output_size(const std::string &option) {
  std::set<std::string> choices;
  std::string keyword;
  bool flag;
  choices.insert("none");
  choices.insert("small");
  choices.insert("medium");
  choices.insert("big");
  OptionsList(option, "UNKNOWN", choices, "UNKNOWN", keyword, flag);
  return keyword;
}


} // end of namespace pism
