// Copyright (C) 2004--2017 Jed Brown, Ed Bueler and Constantine Khroulev
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
void IceModel::process_options() {

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

//! Assembles a list of diagnostics corresponding to an output file size.
std::set<std::string> IceModel::output_variables(const std::string &keyword) {

  std::string variables;

  if (keyword == "none" or
      keyword == "small") {
    variables = "";
  } else if (keyword == "medium") {
    variables = m_config->get_string("output.sizes.medium");
  } else if (keyword == "big_2d") {
    variables = (m_config->get_string("output.sizes.medium") + "," +
                 m_config->get_string("output.sizes.big_2d"));
  } else if (keyword == "big") {
    variables = (m_config->get_string("output.sizes.medium") + "," +
                 m_config->get_string("output.sizes.big_2d") + "," +
                 m_config->get_string("output.sizes.big"));
  }

  return set_split(variables, ',');
}

} // end of namespace pism
