// Copyright (C) 2011, 2012, 2013 PISM Authors
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

#ifndef _PISM_OPTIONS_H_
#define _PISM_OPTIONS_H_

#include "pism_const.hh"

class NCConfigVariable;

PetscErrorCode verbosityLevelFromOptions();

// handy functions for processing options:
PetscErrorCode PISMOptionsList(MPI_Comm com, std::string opt, std::string text, std::set<std::string> choices,
			       std::string default_value, std::string &result, bool &flag);

PetscErrorCode PISMOptionsString(std::string option, std::string text,
				 std::string &result, bool &flag, bool allow_empty_arg = false);
PetscErrorCode PISMOptionsStringArray(std::string opt, std::string text, std::string default_value,
				      std::vector<std::string>& result, bool &flag);
PetscErrorCode PISMOptionsStringSet(std::string opt, std::string text, std::string default_value,
                                    std::set<std::string>& result, bool &flag);

PetscErrorCode PISMOptionsInt(std::string option, std::string text,
			      PetscInt &result, bool &is_set);
PetscErrorCode PISMOptionsIntArray(std::string option, std::string text,
				   std::vector<PetscInt> &result, bool &is_set);

PetscErrorCode PISMOptionsReal(std::string option, std::string text,
			       PetscReal &result, bool &is_set);
PetscErrorCode PISMOptionsRealArray(std::string option, std::string text,
				    std::vector<PetscReal> &result, bool &is_set);

PetscErrorCode PISMOptionsIsSet(std::string option, bool &result);
PetscErrorCode PISMOptionsIsSet(std::string option, std::string descr, bool &result);

PetscErrorCode PISMOptionsHasArgument(std::string option, bool &result);

PetscErrorCode ignore_option(MPI_Comm com, std::string name);
PetscErrorCode check_old_option_and_stop(MPI_Comm com, std::string old_name, std::string new_name);
PetscErrorCode stop_if_set(MPI_Comm com, std::string name);

// usage message and required options; drivers use these
PetscErrorCode stop_on_version_option();

PetscErrorCode show_usage_and_quit(MPI_Comm com, std::string execname, std::string usage);

PetscErrorCode show_usage_check_req_opts(MPI_Comm com, std::string execname,
                                         std::vector<std::string> required_options,
                                         std::string usage);

// config file initialization:
PetscErrorCode init_config(MPI_Comm com, PetscMPIInt rank,
			   NCConfigVariable &config, NCConfigVariable &overrides,
                           bool process_options = false);

PetscErrorCode set_config_from_options(MPI_Comm com,
                                       NCConfigVariable &config);


#endif /* _PISM_OPTIONS_H_ */
