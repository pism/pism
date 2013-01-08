// Copyright (C) 2011, 2012, 2013 PISM Authors
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

#ifndef _PISM_OPTIONS_H_
#define _PISM_OPTIONS_H_

#include "pism_const.hh"

// use namespace std BUT remove trivial namespace browser from doxygen-erated HTML source browser
/// @cond NAMESPACE_BROWSER
using namespace std;
/// @endcond

class NCConfigVariable;

PetscErrorCode verbosityLevelFromOptions();

// handy functions for processing options:
PetscErrorCode PISMOptionsList(MPI_Comm com, string opt, string text, set<string> choices,
			       string default_value, string &result, bool &flag);

PetscErrorCode PISMOptionsString(string option, string text,
				 string &result, bool &flag, bool allow_empty_arg = false);
PetscErrorCode PISMOptionsStringArray(string opt, string text, string default_value,
				      vector<string>& result, bool &flag);
PetscErrorCode PISMOptionsStringSet(string opt, string text, string default_value,
                                    set<string>& result, bool &flag);

PetscErrorCode PISMOptionsInt(string option, string text,
			      PetscInt &result, bool &is_set);
PetscErrorCode PISMOptionsIntArray(string option, string text,
				   vector<PetscInt> &result, bool &is_set);

PetscErrorCode PISMOptionsReal(string option, string text,
			       PetscReal &result, bool &is_set);
PetscErrorCode PISMOptionsRealArray(string option, string text,
				    vector<PetscReal> &result, bool &is_set);

PetscErrorCode PISMOptionsIsSet(string option, bool &result);
PetscErrorCode PISMOptionsIsSet(string option, string descr, bool &result);

PetscErrorCode ignore_option(MPI_Comm com, const char name[]);
PetscErrorCode check_old_option_and_stop(MPI_Comm com, const char old_name[], const char new_name[]);
PetscErrorCode stop_if_set(MPI_Comm com, const char name[]);
PetscErrorCode parse_range(MPI_Comm com, string str, double *a, double *delta, double *b, string &keyword);
PetscErrorCode parse_times(MPI_Comm com, const NCConfigVariable &config, string str,
                           double run_start, double run_end, vector<double> &result);

// usage message and required options; drivers use these
PetscErrorCode stop_on_version_option();

PetscErrorCode show_usage_and_quit(MPI_Comm com, const char execname[], const char usage[]);

PetscErrorCode show_usage_check_req_opts(MPI_Comm com, const char execname[],
                                         vector<string> required_options,
                                         const char usage[]);

// config file initialization:
PetscErrorCode init_config(MPI_Comm com, PetscMPIInt rank,
			   NCConfigVariable &config, NCConfigVariable &overrides,
                           bool process_options = false);

PetscErrorCode set_config_from_options(MPI_Comm com,
                                       NCConfigVariable &config);


#endif /* _PISM_OPTIONS_H_ */
