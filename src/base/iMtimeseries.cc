// Copyright (C) 2009 Constantine Khroulev
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

#include "iceModel.hh"
#include <sstream>
#include <algorithm>

//! Initializes the code writing scalar time-series.
PetscErrorCode IceModel::init_timeseries() {
  PetscErrorCode ierr;
  PetscTruth ts_file_set = PETSC_FALSE, ts_times_set = PETSC_FALSE, ts_vars_set = PETSC_FALSE;
  char tmp[TEMPORARY_STRING_LENGTH] = "\0";

  ierr = PetscOptionsGetString(PETSC_NULL, "-ts_file", tmp,
			       PETSC_MAX_PATH_LEN, &ts_file_set); CHKERRQ(ierr);
  ts_filename = tmp;

  ierr = PetscOptionsGetString(PETSC_NULL, "-ts_times", tmp,
			       TEMPORARY_STRING_LENGTH, &ts_times_set); CHKERRQ(ierr);

  if (ts_file_set ^ ts_times_set) {
    ierr = PetscPrintf(grid.com,
		       "PISM ERROR: you need to specity both -ts_filename and -ts_times to save"
		       "diagnostic time-seties.\n");
    CHKERRQ(ierr);
    PetscEnd();
  }

  // If neither -ts_filename nor -ts_times is set, we're done.
  if (!ts_file_set && !ts_times_set) {
    save_ts = false;
    return 0;
  }
  
  save_ts = true;

  ierr = parse_times(grid.com, tmp, ts_times);
  if (ierr != 0) {
    ierr = PetscPrintf(grid.com, "PISM ERROR: parsing the -ts_times argument failed.\n"); CHKERRQ(ierr);
    PetscEnd();
  }

  ierr = verbPrintf(2, grid.com, "saving scalar time-series to '%s'; ",
		    ts_filename.c_str()); CHKERRQ(ierr);

  ierr = verbPrintf(2, grid.com, "times requested: %s\n", tmp); CHKERRQ(ierr);

  current_ts = 0;

  ierr = PetscOptionsGetString(PETSC_NULL, "-ts_vars", tmp,
			       TEMPORARY_STRING_LENGTH, &ts_vars_set); CHKERRQ(ierr);
  string var_name;
  if (ts_vars_set) {
    ierr = verbPrintf(2, grid.com, "variables requested: %s\n", tmp); CHKERRQ(ierr);
    istringstream arg(tmp);

    while (getline(arg, var_name, ','))
      ts_vars.insert(var_name);

  } else {
    var_name = config.get_string("ts_variables");
    istringstream arg(var_name);
  
    while (getline(arg, var_name, ' ')) {
      if (!var_name.empty()) // this ignores multiple spaces separating variable names
	ts_vars.insert(var_name);
    }
  }

  // This will move the file aside if it exists already.
  NCTool nc(grid.com, grid.rank);
  ierr = nc.open_for_writing(ts_filename.c_str(), false, false); CHKERRQ(ierr);
  ierr = nc.close(); CHKERRQ(ierr);

  ierr = create_timeseries(); CHKERRQ(ierr);
  
  return 0;
}

//! \brief Creates DiagnosticTimeseries objects used to store and report scalar
//! diagnostic quantities.
PetscErrorCode IceModel::create_timeseries() {
  
  if (find(ts_vars.begin(), ts_vars.end(), "ivol") != ts_vars.end()) {
    DiagnosticTimeseries *ivol = new DiagnosticTimeseries(&grid, "ivol", "t");

    ivol->set_units("m3", "");
    ivol->set_dimension_units("years", "");
    ivol->output_filename = ts_filename;

    ivol->set_attr("long_name", "ice volume");
    ivol->set_attr("valid_min", 0.0);

    timeseries.push_back(ivol);
  }

  if (find(ts_vars.begin(), ts_vars.end(), "iarea") != ts_vars.end()) {
    DiagnosticTimeseries *iarea = new DiagnosticTimeseries(&grid, "iarea", "t");

    iarea->set_units("m2", "");
    iarea->set_dimension_units("years", "");
    iarea->output_filename = ts_filename;

    iarea->set_attr("long_name", "ice area");
    iarea->set_attr("valid_min", 0.0);

    timeseries.push_back(iarea);
  }

  if (find(ts_vars.begin(), ts_vars.end(), "iareag") != ts_vars.end()) {
    DiagnosticTimeseries *iareag = new DiagnosticTimeseries(&grid, "iareag", "t");

    iareag->set_units("m2", "");
    iareag->set_dimension_units("years", "");
    iareag->output_filename = ts_filename;

    iareag->set_attr("long_name", "grounded ice area");
    iareag->set_attr("valid_min", 0.0);

    timeseries.push_back(iareag);
  }

  if (find(ts_vars.begin(), ts_vars.end(), "iareaf") != ts_vars.end()) {
    DiagnosticTimeseries *iareaf = new DiagnosticTimeseries(&grid, "iareaf", "t");

    iareaf->set_units("m2", "");
    iareaf->set_dimension_units("years", "");
    iareaf->output_filename = ts_filename;

    iareaf->set_attr("long_name", "floating ice area");
    iareaf->set_attr("valid_min", 0.0);

    timeseries.push_back(iareaf);
  }

  if (find(ts_vars.begin(), ts_vars.end(), "dt") != ts_vars.end()) {
    DiagnosticTimeseries *delta_t = new DiagnosticTimeseries(&grid, "dt", "t");

    delta_t->set_units("s", "years");
    delta_t->set_dimension_units("years", "");
    delta_t->output_filename = ts_filename;

    delta_t->set_attr("long_name", "mass continuity time-step");
    delta_t->set_attr("valid_min", 0.0);

    timeseries.push_back(delta_t);
  }

  return 0;
}

//! Write time-series.
PetscErrorCode IceModel::write_timeseries() {
  PetscErrorCode ierr;

  // return if no time-series requested
  if (!save_ts) return 0;

  // return if wrote all the records already
  if (current_ts == ts_times.size())
    return 0;

  // return if did not yet reach the time we need to save at
  if (ts_times[current_ts] > grid.year)
    return 0;

  // compute values of requested scalar quantities:
  vector<DiagnosticTimeseries*>::iterator i;
  for (i = timeseries.begin(); i < timeseries.end(); ++i) {
    PetscScalar tmp;

    ierr = compute_by_name((*i)->short_name, tmp); CHKERRQ(ierr);
    
    (*i)->append(grid.year, tmp);
  }

  // Interpolate to put them on requested times:
  while ((ts_times[current_ts] <= grid.year) &&
	 (current_ts < ts_times.size())) {
    
    vector<DiagnosticTimeseries*>::iterator i;
    for (i = timeseries.begin(); i < timeseries.end(); ++i) {
      ierr = (*i)->interp(ts_times[current_ts]); CHKERRQ(ierr);
    }

    current_ts++;
  }

  return 0;
}


//! Initialize the code saving sptially-variable diagnostic quantities.
PetscErrorCode IceModel::init_extras() {
  PetscErrorCode ierr;
  PetscTruth times = PETSC_FALSE, file = PETSC_FALSE;
  char tmp[TEMPORARY_STRING_LENGTH] = "\0";
  current_extra = 0;

  ierr = PetscOptionsGetString(PETSC_NULL, "-extra_file", tmp,
			       PETSC_MAX_PATH_LEN, &file); CHKERRQ(ierr);
  extra_filename = tmp;

  ierr = PetscOptionsGetString(PETSC_NULL, "-extra_times", tmp,
			       TEMPORARY_STRING_LENGTH, &times); CHKERRQ(ierr);

  if (file ^ times) {
    ierr = PetscPrintf(grid.com,
		       "PISM ERROR: you need to specify both -extra_file and -extra_times to save spatial time-series.\n");
    CHKERRQ(ierr);
    PetscEnd();
  }

  if (!file && !times) {
    save_extra = false;
    return 0;
  }

  ierr = parse_times(grid.com, tmp, extra_times);
  if (ierr != 0) {
    ierr = PetscPrintf(grid.com, "PISM ERROR: parsing the -extra_times argument failed.\n"); CHKERRQ(ierr);
    PetscEnd();
  }

  save_extra = true;
  extra_file_is_ready = false;
  split_extra = false;

  PetscTruth split;
  ierr = check_option("-split_extras", split); CHKERRQ(ierr);
  if (split) {
    split_extra = true;
  } else if (!ends_with(extra_filename, ".nc")) {
    ierr = verbPrintf(2, grid.com,
		      "PISM WARNING: spatial time-series file name '%s' does not have the '.nc' suffix!\n",
		      extra_filename.c_str());
    CHKERRQ(ierr);
  }
  
  if (split) {
    ierr = verbPrintf(2, grid.com, "saving spatial time-series to '%s+year.nc'; ",
		      extra_filename.c_str()); CHKERRQ(ierr);
  } else {
    ierr = verbPrintf(2, grid.com, "saving spatial time-series to '%s'; ",
		      extra_filename.c_str()); CHKERRQ(ierr);
  }

  if (extra_times.size() > 500) {
    ierr = verbPrintf(2, grid.com,
		      "PISM WARNING: more than 500 times requested. This might fill your hard-drive!\n");
    CHKERRQ(ierr);
  }

  ierr = verbPrintf(2, grid.com, "times requested: %s\n", tmp); CHKERRQ(ierr);

  PetscTruth save_vars;
  ierr = PetscOptionsGetString(PETSC_NULL, "-extra_vars", tmp,
			       TEMPORARY_STRING_LENGTH, &save_vars); CHKERRQ(ierr);

  string var_name;
  if (save_vars) {
    ierr = verbPrintf(2, grid.com, "variables requested: %s\n", tmp); CHKERRQ(ierr);
    istringstream arg(tmp);

    while (getline(arg, var_name, ','))
      extra_vars.insert(var_name);

  } else {
    ierr = verbPrintf(2, grid.com, "PISM WARNING: -extra_vars was not set. Writing the standard set of output variables...\n"); CHKERRQ(ierr);

    var_name = config.get_string("output_variables");
    istringstream arg(var_name);
  
    while (getline(arg, var_name, ' ')) {
      if (!var_name.empty()) // this ignores multiple spaces separating variable names
	extra_vars.insert(var_name);

    }
  }

  return 0;
}

//! Write spatiall-variable diagnostic quantities.
PetscErrorCode IceModel::write_extras() {
  PetscErrorCode ierr;
  NCTool nc(&grid);
  double saving_after = -1.0e30; // initialize to avoid compiler warning; this
				 // value is never used, because saving_after
				 // is only used if save_now == true, and in
				 // this case saving_after is guaranteed to be
				 // initialized. See the code below.
  char filename[PETSC_MAX_PATH_LEN];

  // determine if the user set the -save_at and -save_to options
  if (!save_extra)
    return 0;

  // do we need to save *now*?
  if ( (grid.year >= extra_times[current_extra]) &&
       (current_extra < extra_times.size()) ) {
    saving_after = extra_times[current_extra];

    while (extra_times[current_extra] <= grid.year)
      current_extra++;
  } else {
    // we don't need to save now, so just return
    return 0;
  }

  if (split_extra) {
    extra_file_is_ready = false;	// each time-series record is written to a separate file
    snprintf(filename, PETSC_MAX_PATH_LEN, "%s-%06.0f.nc",
	     extra_filename.c_str(), grid.year);
  } else {
    strncpy(filename, extra_filename.c_str(), PETSC_MAX_PATH_LEN);
  }

  ierr = verbPrintf(3, grid.com, 
		    "\nsaving spatial time-series to %s at %.5f a, for time-step goal %.5f a\n\n",
		    filename, grid.year,saving_after);
  CHKERRQ(ierr);

  // create line for history in .nc file, including time of write
  string date_str = timestamp();
  char tmp[TEMPORARY_STRING_LENGTH];
  snprintf(tmp, TEMPORARY_STRING_LENGTH,
	   "%s: %s saving time-series record at %10.5f a, for time-step goal %10.5f a\n",
	   date_str.c_str(), executable_short_name.c_str(), grid.year, saving_after);

  if (!extra_file_is_ready) {

    // Prepare the file:
    ierr = nc.open_for_writing(filename, false, true); CHKERRQ(ierr);
    // append == false, check_dims == true
    ierr = nc.write_history(history.c_str()); CHKERRQ(ierr); // append the history
    ierr = nc.write_global_attrs(false, "CF-1.4"); CHKERRQ(ierr);
    ierr = nc.close(); CHKERRQ(ierr);

    ierr = polar_stereographic.write(filename); CHKERRQ(ierr);
    extra_file_is_ready = true;
  }
    
  ierr = nc.open_for_writing(filename, true, true); CHKERRQ(ierr);
  // append == true, check_dims == true
  ierr = nc.append_time(grid.year); CHKERRQ(ierr);
  ierr = nc.write_history(tmp); CHKERRQ(ierr); // append the history
  ierr = nc.close(); CHKERRQ(ierr);

  ierr = write_variables(filename, extra_vars);  CHKERRQ(ierr);

  return 0;

  /* FIXME: should we let couplers write their fields?

  if (atmosPCC != PETSC_NULL) {
    ierr = atmosPCC->writeCouplingFieldsToFile(grid.year,filename); CHKERRQ(ierr);
  } else {
    SETERRQ(1,"PISM ERROR: atmosPCC == PETSC_NULL");
  }

  if (oceanPCC != PETSC_NULL) {
    ierr = oceanPCC->writeCouplingFieldsToFile(grid.year,filename); CHKERRQ(ierr);
  } else {
    SETERRQ(1,"PISM ERROR: oceanPCC == PETSC_NULL");
  }
    
  return 0;
  */
}
