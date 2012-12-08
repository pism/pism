// Copyright (C) 2004--2012 Jed Brown, Ed Bueler and Constantine Khroulev
//
// This file is part of PISM.
//
// PISM is free software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the Free Software
// Foundation; either version 2 of the License, or (at your option) any later
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
#include "PISMHydrology.hh"
#include "PISMStressBalance.hh"
#include "PISMOcean.hh"
#include "PISMSurface.hh"
#include "pism_options.hh"
#include "IceGrid.hh"

//! \file iMoptions.cc Reading runtime options and setting configuration parameters.


//! Read runtime (command line) options and alter the corresponding parameters or flags as appropriate.
/*!
A critical principle of this procedure is that it will not alter IceModel parameters and flags
\e unless the user sets an option to do so.  This base class setFromOptions() can be
called by an IceModel-derived class after the it has set its own defaults.

In fact this procedure only reads the majority of the options.  Some are read in 
initFromOptions(), writeFiles(), and setStartRunEndYearsFromOptions(), among other places.

There are no options to directly set \c dx, \c dy, \c dz, \c Lbz, and \c year as the user 
should not directly set these grid parameters.  There are, however, options for directly 
setting \c Mx, \c My, \c Mz, \c Mbz and also \c Lx, \c Ly, \c Lz.

Note that additional options are read by PISM{Atmosphere|Surface|Ocean}Model
instances, including -pdd... others.
 */
PetscErrorCode  IceModel::setFromOptions() {
  PetscErrorCode ierr;
  bool flag;

  ierr = verbPrintf(3, grid.com,
		    "Processing physics-related command-line options...\n"); CHKERRQ(ierr);

  ierr = PetscOptionsBegin(grid.com, "", "Options overriding config flags and parameters", ""); CHKERRQ(ierr);

  ierr = set_config_from_options(grid.com, config); CHKERRQ(ierr);

  ierr = PISMOptionsInt("-id", "Specifies the sounding row", id, flag); CHKERRQ(ierr);
  ierr = PISMOptionsInt("-jd", "Specifies the sounding column", jd, flag); CHKERRQ(ierr);

  bool initfromT, initfromTandOm;
  ierr = PISMOptionsIsSet("-init_from_temp", initfromT); CHKERRQ(ierr);
  ierr = PISMOptionsIsSet("-init_from_temp_and_liqfrac", initfromTandOm); CHKERRQ(ierr);

  ierr = PetscOptionsEnd(); CHKERRQ(ierr);

  // Set global attributes using the config database:
  global_attributes.set_from_config(config);

  // warn about some option combinations

  if (config.get("maximum_time_step_years") <= 0) {
    PetscPrintf(grid.com, "PISM ERROR: maximum_time_step_years has to be greater than 0.\n");
    PISMEnd();
  }
  
  if (config.get_flag("do_mass_conserve") == false &&
      config.get_flag("do_skip")) {
    ierr = verbPrintf(2, grid.com,
      "PISM WARNING: Both -skip and -no_mass are set.\n"
      "              -skip only makes sense in runs updating ice geometry.\n"); CHKERRQ(ierr);
  }

  if (config.get_flag("do_thickness_calving") &&
      config.get_flag("part_grid") == false) {
    ierr = verbPrintf(2, grid.com,
      "PISM WARNING: Calving at certain terminal ice thickness (-calving_at_thickness)\n"
      "              without application of partially filled grid cell scheme (-part_grid)\n"
      "              may lead to (incorrect) non-moving ice shelf front.\n"); CHKERRQ(ierr);
  }


  // implements an option e.g. described in \ref Greve that is the
  // enhancement factor is coupled to the age of the ice with
  // e = 1 (A < 11'000 years), e = 3 otherwise
  if (config.get_flag("do_e_age_coupling")) {
    ierr = verbPrintf(2, grid.com,
		      "  setting age-dependent enhancement factor: "
                      "e=1 if A<11'000 years, e=3 otherwise\n"); CHKERRQ(ierr);

  }


  // old options
  ierr = check_old_option_and_stop(grid.com, "-eta",    "-gradient"); CHKERRQ(ierr);
  ierr = check_old_option_and_stop(grid.com, "-no_eta", "-gradient"); CHKERRQ(ierr);
  ierr = check_old_option_and_stop(grid.com, "-ssa",
				   "-ssa_sliding' or '-ssa_floating_only"); CHKERRQ(ierr);
  ierr = check_old_option_and_stop(grid.com, "-plastic", "-ssa_sliding"); CHKERRQ(ierr);
  ierr = check_old_option_and_stop(grid.com, "-sliding_scale_brutal",
                                   "-brutal_sliding' and '-brutal_sliding_scale"); CHKERRQ(ierr);

  return 0;
}

//! Assembles a list of variables corresponding to an output file size.
PetscErrorCode IceModel::set_output_size(string option,
					 string description,
					 string default_value,
					 set<string> &result) {
  PetscErrorCode ierr;
  set<string> choices;
  string keyword;
  bool flag;
  map<string, NCSpatialVariable> list;

  result.clear();

  if (keyword == "none") {
    return 0;
  }

  choices.insert("none");
  choices.insert("small");
  choices.insert("medium");
  choices.insert("big");
  ierr = PISMOptionsList(grid.com, option,
			 description, choices,
			 default_value, keyword, flag); CHKERRQ(ierr);

  // Add all the model-state variables:
  set<string> vars = variables.keys();

  set<string>::iterator i = vars.begin();
  while (i != vars.end()) {
    IceModelVec *var = variables.get(*i);

    string intent = var->string_attr("pism_intent");
    if ((intent == "model_state") || (intent == "mapping") || (intent == "climate_steady")) {
      result.insert(*i);
    }
    i++;
  }

  if (keyword == "medium") {
    // add all the variables listed in the config file ("medium" size):
    string tmp = config.get_string("output_medium");
    istringstream keywords(tmp);

    // split the list; note that this also removes any duplicate entries
    while (getline(keywords, tmp, ' ')) {
      if (!tmp.empty())                // this ignores multiple spaces separating variable names
       result.insert(tmp);
    }
  } else if (keyword == "big") {
    // add all the variables listed in the config file ("big" size):
    string tmp = config.get_string("output_big");
    istringstream keywords(tmp);

    // split the list; note that this also removes any duplicate entries
    while (getline(keywords, tmp, ' ')) {
      if (!tmp.empty())                // this ignores multiple spaces separating variable names
       result.insert(tmp);
    }

    if (!config.get_flag("do_age"))
      result.erase("age");
  }

  if (config.get_flag("do_age"))
    result.insert("age");

  if (config.get_flag("ocean_kill"))
    result.insert("ocean_kill_mask");

  if (beddef != NULL)
    beddef->add_vars_to_output(keyword, list);

  if (btu != NULL)
    btu->add_vars_to_output(keyword, list);

  if (basal_yield_stress != NULL)
    basal_yield_stress->add_vars_to_output(keyword, list);

  // Ask the stress balance module to add more variables:
  if (stress_balance != NULL)
    stress_balance->add_vars_to_output(keyword, list);

  if (subglacial_hydrology != NULL)
    subglacial_hydrology->add_vars_to_output(keyword, list);

  // Ask ocean and surface models to add more variables to the list:
  if (ocean != NULL)
    ocean->add_vars_to_output(keyword, list);

  if (surface != NULL)
    surface->add_vars_to_output(keyword, list);

  map<string,NCSpatialVariable>::iterator j = list.begin();
  while(j != list.end()) {
    result.insert(j->first);
    ++j;
  }

  return 0;
}


//! Returns the output size as a keyword, for options "-o_size", "-save_size", "-backup_size", etc.
string IceModel::get_output_size(string option) {
  set<string> choices;
  string keyword;
  bool flag;
  choices.insert("none");
  choices.insert("small");
  choices.insert("medium");
  choices.insert("big");
  PISMOptionsList(grid.com, option,
		  "UNKNOWN", choices,
		  "UNKNOWN", keyword, flag);
  return keyword;
}

