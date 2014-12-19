// Copyright (C) 2009-2014 Constantine Khroulev
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

#include <sstream>
#include <algorithm>

#include "iceModel.hh"
#include "PIO.hh"
#include "PISMStressBalance.hh"
#include "PISMDiagnostic.hh"
#include "PISMTime.hh"
#include "pism_options.hh"

#include "error_handling.hh"

namespace pism {

//! Initializes the code writing scalar time-series.
void IceModel::init_timeseries() {

  options::String ts_file("-ts_file", "Specifies the time-series output file name");
  ts_filename = ts_file;

  options::String times("-ts_times", "Specifies a MATLAB-style range or a list of requested times");

  options::StringSet vars("-ts_vars", "Specifies a comma-separated list of veriables to save",
                          "");

  // default behavior is to move the file aside if it exists already; option allows appending
  options::Bool append("-ts_append", "append scalar time-series");


  IO_Mode mode = PISM_READWRITE;
  if (not append.is_set()) {
    mode = PISM_READWRITE_MOVE;
  }

  if (ts_file.is_set() ^ times.is_set()) {
    throw RuntimeError("you need to specity both -ts_file and -ts_times to save diagnostic time-series.");
  }

  // If neither -ts_filen nor -ts_times is set, we're done.
  if (not ts_file.is_set() && not times.is_set()) {
    save_ts = false;
    return;
  }

  save_ts = true;

  try {
    grid.time->parse_times(times, ts_times);  
  } catch (RuntimeError &e) {
    e.add_context("parsing the -ts_times argument");
    throw;
  }

  if (times.value().size() == 0) {
    throw RuntimeError("no argument for -ts_times option.");
  }

  verbPrintf(2, grid.com, "saving scalar time-series to '%s'; ",
             ts_file.c_str());

  verbPrintf(2, grid.com, "times requested: %s\n", times.c_str());

  current_ts = 0;


  std::string var_name;
  if (vars.is_set()) {
    verbPrintf(2, grid.com, "variables requested: %s\n", vars.print().c_str());
    ts_vars = vars;
  } else {
    std::map<std::string,TSDiagnostic*>::iterator j = ts_diagnostics.begin();
    while (j != ts_diagnostics.end()) {
      ts_vars.insert(j->first);
      ++j;
    }
  }

  PIO nc(grid, "netcdf3");      // Use NetCDF-3 to write time-series.
  nc.open(ts_file, mode);

  if (append == true) {
    double time_max;
    std::string time_name = config.get_string("time_dimension_name");
    bool time_exists = false;

    time_exists = nc.inq_var(time_name);
    if (time_exists == true) {
      nc.inq_dim_limits(time_name, NULL, &time_max);

      while (current_ts < ts_times.size() && ts_times[current_ts] < time_max) {
        current_ts++;
      }

      if (current_ts > 0) {
        verbPrintf(2, grid.com,
                   "skipping times before the last record in %s (at %s)\n",
                   ts_file.c_str(), grid.time->date(time_max).c_str());
      }
    }
  }

  write_metadata(nc, false, false);

  nc.close();


  // set the output file:
  std::map<std::string,TSDiagnostic*>::iterator j = ts_diagnostics.begin();
  while (j != ts_diagnostics.end()) {
    (j->second)->init(ts_file);
    ++j;
  }

  // ignore times before (and including) the beginning of the run:
  while (current_ts < ts_times.size() && ts_times[current_ts] < grid.time->start()) {
    current_ts++;
  }

  if (ts_times.size() == current_ts) {
    save_ts = false;
    return;
  }

  // discard requested times before the beginning of the run
  std::vector<double> tmp(ts_times.size() - current_ts);
  for (unsigned int k = 0; k < tmp.size(); ++k) {
    tmp[k] = ts_times[current_ts + k];
  }

  ts_times = tmp;
  current_ts = 0;
}

//! Write time-series.
void IceModel::write_timeseries() {

  // return if no time-series requested
  if (!save_ts) {
     return;
  }

  // return if wrote all the records already
  if (current_ts == ts_times.size()) {
    return;
  }

  // return if did not yet reach the time we need to save at
  if (ts_times[current_ts] > grid.time->current()) {
    return;
  }

  for (std::set<std::string>::iterator j = ts_vars.begin(); j != ts_vars.end(); ++j) {
    TSDiagnostic *diag = ts_diagnostics[*j];

    if (diag != NULL) {
      diag->update(grid.time->current() - dt, grid.time->current());
    }
  }


  // Interpolate to put them on requested times:
  for (; current_ts < ts_times.size() && ts_times[current_ts] <= grid.time->current(); current_ts++) {

    // the very first time (current_ts == 0) defines the left endpoint of the
    // first time interval; we don't write a report at that time
    if (current_ts == 0) {
      continue;
    }

    for (std::set<std::string>::iterator j = ts_vars.begin(); j != ts_vars.end(); ++j) {
      TSDiagnostic *diag = ts_diagnostics[*j];

      if (diag != NULL) {
        diag->save(ts_times[current_ts - 1], ts_times[current_ts]);
      }
    }
  }
}


//! Initialize the code saving spatially-variable diagnostic quantities.
void IceModel::init_extras() {
  bool extra_times_set, extra_file_set, extra_vars_set;
  std::string times, vars;

  last_extra = 0;               // will be set in write_extras()
  next_extra = 0;


  OptionsString("-extra_file", "Specifies the output file",
                extra_filename, extra_file_set);

  OptionsString("-extra_times", "Specifies times to save at",
                times, extra_times_set);

  OptionsString("-extra_vars", "Specifies a comma-separated list of variables to save",
                vars, extra_vars_set);

  options::Bool split("-extra_split", "Specifies whether to save to separate files");
  options::Bool append("-extra_append", "append spatial diagnostics");

  if (extra_file_set ^ extra_times_set) {
    throw RuntimeError("you need to specify both -extra_file and -extra_times to save spatial time-series.");
  }

  if (!extra_file_set && !extra_times_set) {
    save_extra = false;
    return;
  }

  try {
    grid.time->parse_times(times, extra_times);    
  } catch (RuntimeError &e) {
    e.add_context("parsing the -extra_times argument");
    throw;
  }

  if (extra_times.size() == 0) {
    throw RuntimeError("no argument for -extra_times option.");
  }

  if (append.is_set() && split.is_set()) {
    throw RuntimeError("both -extra_split and -extra_append are set.");
  }

  if (append) {
    PIO nc(grid, grid.config.get_string("output_format"));
    std::string time_name = config.get_string("time_dimension_name");
    bool time_exists;

    nc.open(extra_filename, PISM_READONLY);
    time_exists = nc.inq_var(time_name);

    if (time_exists == true) {
      double time_max;
      nc.inq_dim_limits(time_name, NULL, &time_max);

      while (next_extra + 1 < extra_times.size() && extra_times[next_extra + 1] < time_max) {
        next_extra++;
      }

      if (next_extra > 0) {
        verbPrintf(2, grid.com,
                   "skipping times before the last record in %s (at %s)\n",
                   extra_filename.c_str(), grid.time->date(time_max).c_str());
      }

      // discard requested times before the beginning of the run
      std::vector<double> tmp(extra_times.size() - next_extra);
      for (unsigned int k = 0; k < tmp.size(); ++k) {
        tmp[k] = extra_times[next_extra + k];
      }

      extra_times = tmp;
      next_extra = 0;
    }
    nc.close();
  }

  save_extra = true;
  extra_file_is_ready = false;
  split_extra = false;

  if (split.is_set()) {
    split_extra = true;
    verbPrintf(2, grid.com, "saving spatial time-series to '%s+year.nc'; ",
               extra_filename.c_str());
  } else if (!ends_with(extra_filename, ".nc")) {
    verbPrintf(2, grid.com,
               "PISM WARNING: spatial time-series file name '%s' does not have the '.nc' suffix!\n",
               extra_filename.c_str());
    verbPrintf(2, grid.com, "saving spatial time-series to '%s'; ",
               extra_filename.c_str());
  }

  verbPrintf(2, grid.com, "times requested: %s\n", times.c_str());

  if (extra_times.size() > 500) {
    verbPrintf(2, grid.com,
               "PISM WARNING: more than 500 times requested. This might fill your hard-drive!\n");
  }

  std::string var_name;
  if (extra_vars_set) {
    verbPrintf(2, grid.com, "variables requested: %s\n", vars.c_str());
    std::istringstream arg(vars);

    while (getline(arg, var_name, ',')) {
      extra_vars.insert(var_name);
    }

  } else {
    verbPrintf(2, grid.com, "PISM WARNING: -extra_vars was not set."
               " Writing model_state, mapping and climate_steady variables...\n");

    std::set<std::string> vars_set = grid.variables().keys();

    std::set<std::string>::iterator i = vars_set.begin();
    while (i != vars_set.end()) {
      IceModelVec *var = grid.variables().get(*i);
      NCSpatialVariable &m = var->metadata();

      std::string intent = m.get_string("pism_intent");
      if ((intent == "model_state") ||
          (intent == "mapping") ||
          (intent == "climate_steady")) {
        extra_vars.insert(*i);
      }
      i++;
    }

    std::set<std::string> list;
    if (stress_balance) {
      stress_balance->add_vars_to_output("small", extra_vars);
    }

  } // end of the else clause after "if (extra_vars_set)"

  if (extra_vars.size() == 0) {
    verbPrintf(2, grid.com, 
               "PISM WARNING: no variables list after -extra_vars ... writing empty file ...\n");
  }
}

//! Write spatially-variable diagnostic quantities.
void IceModel::write_extras() {
  double saving_after = -1.0e30; // initialize to avoid compiler warning; this
                                 // value is never used, because saving_after
                                 // is only used if save_now == true, and in
                                 // this case saving_after is guaranteed to be
                                 // initialized. See the code below.
  char filename[PETSC_MAX_PATH_LEN];
  unsigned int current_extra;
  // determine if the user set the -save_at and -save_to options
  if (!save_extra) {
    return;
  }

  // do we need to save *now*?
  if (next_extra < extra_times.size() &&
      (grid.time->current() >= extra_times[next_extra] ||
       fabs(grid.time->current() - extra_times[next_extra]) < 1.0)) {
    // the condition above is "true" if we passed a requested time or got to
    // within 1 second from it

    current_extra = next_extra;

    // update next_extra
    while (next_extra < extra_times.size() &&
           (extra_times[next_extra] <= grid.time->current() ||
            fabs(grid.time->current() - extra_times[next_extra]) < 1.0)) {
      next_extra++;
    }

    saving_after = extra_times[current_extra];
  } else {
    return;
  }

  if (current_extra == 0) {
    // The first time defines the left end-point of the first reporting
    // interval; we don't write a report at this time, but we still need to
    // store cumulative quantities that may be needed to compute rates of
    // change.

    std::set<std::string>::iterator j = extra_vars.begin();
    while(j != extra_vars.end()) {
      Diagnostic *diag = diagnostics[*j];

      if (diag != NULL) {
        diag->update_cumulative();
      }
      ++j;
    }

    // This line re-initializes last_extra (the correct value is not known at
    // the time init_extras() is calles).
    last_extra = grid.time->current();

    return;
  }

  if (saving_after < grid.time->start()) {
    // Suppose a user tells PISM to write data at times 0:1000:10000. Suppose
    // also that PISM writes a backup file at year 2500 and gets stopped.
    //
    // When restarted, PISM will decide that it's time to write data for time
    // 2000, but
    // * that record was written already and
    // * PISM will end up writing at year 2500, producing a file containing one
    //   more record than necessary.
    //
    // This check makes sure that this never happens.
    return;
  }

  if (split_extra) {
    extra_file_is_ready = false;        // each time-series record is written to a separate file
    snprintf(filename, PETSC_MAX_PATH_LEN, "%s-%s.nc",
             extra_filename.c_str(), grid.time->date().c_str());
  } else {
    strncpy(filename, extra_filename.c_str(), PETSC_MAX_PATH_LEN);
  }

  verbPrintf(3, grid.com, 
             "\nsaving spatial time-series to %s at %s\n\n",
             filename, grid.time->date().c_str());

  // find out how much time passed since the beginning of the run
  double wall_clock_hours;
  if (grid.rank() == 0) {
    PetscLogDouble current_time;
    GetTime(&current_time);
    wall_clock_hours = (current_time - start_time) / 3600.0;
  }

  MPI_Bcast(&wall_clock_hours, 1, MPI_DOUBLE, 0, grid.com);

  PIO nc(grid, grid.config.get_string("output_format"));

  if (extra_file_is_ready == false) {
    // default behavior is to move the file aside if it exists already; option allows appending
    options::Bool append("-extra_append", "append -extra_file output");

    IO_Mode mode = PISM_READWRITE;
    if (not append.is_set()) {
      mode = PISM_READWRITE_MOVE;
    }

    // Prepare the file:
    nc.open(filename, mode);
    nc.def_time(config.get_string("time_dimension_name"),
                grid.time->calendar(),
                grid.time->CF_units_string());
    nc.put_att_text(config.get_string("time_dimension_name"),
                    "bounds", "time_bounds");

    write_metadata(nc, true, false);

    extra_file_is_ready = true;
  } else {
    // In this case the extra file should be present.
    nc.open(filename, PISM_READWRITE);
  }

  double      current_time = grid.time->current();
  std::string time_name    = config.get_string("time_dimension_name");

  nc.append_time(time_name, current_time);

  unsigned int time_length = nc.inq_dimlen(time_name);

  size_t time_start = static_cast<size_t>(time_length - 1);

  std::vector<double> data(2);
  data[0] = last_extra;
  data[1] = current_time;
  nc.write_time_bounds(extra_bounds, time_start, data);

  nc.write_timeseries(timestamp, time_start, wall_clock_hours);

  write_variables(nc, extra_vars, PISM_FLOAT);

  nc.close();

  // flush time-series buffers
  flush_timeseries();

  last_extra = current_time;
}

//! Computes the maximum time-step we can take and still hit all the requested years.
/*!
  Sets restrict to 'false' if any time-step is OK.
 */
void IceModel::extras_max_timestep(double my_t, double& my_dt, bool &restrict) {

  if (!save_extra) {
    my_dt = -1;
    restrict = false;
    return;
  }

  if (config.get_flag("extras_force_output_times") == false) {
    my_dt = -1;
    restrict = false;
    return;
  }

  std::vector<double>::iterator j;
  j = upper_bound(extra_times.begin(), extra_times.end(), my_t);

  if (j == extra_times.end()) {
    my_dt = -1;
    restrict = false;
    return;
  }

  my_dt = *j - my_t;
  restrict = true;

  // now make sure that we don't end up taking a time-step of less than 1
  // second long
  if (my_dt < 1) {
    if ((j + 1) != extra_times.end()) {
      my_dt = *(j + 1) - my_t;
      restrict = true;
    } else {
      my_dt = -1;
      restrict = false;
    }
  }
}

//! Computes the maximum time-step we can take and still hit all the requested years.
/*!
  Sets restrict to 'false' if any time-step is OK.
 */
void IceModel::ts_max_timestep(double my_t, double& my_dt, bool &restrict) {

  if (!save_ts) {
    my_dt = -1;
    restrict = false;
    return;
  }

  // make sure that we hit the left endpoint of the first report interval
  if (my_t < ts_times[0]) {
    my_dt = ts_times[0] - my_t;
    restrict = true;
    return;
  }

  bool force_times;
  force_times = config.get_flag("ts_force_output_times");

  if (!force_times) {
    my_dt = -1;
    restrict = false;
    return;
  }

  std::vector<double>::iterator j;
  j = upper_bound(ts_times.begin(), ts_times.end(), my_t);

  if (j == ts_times.end()) {
    my_dt = -1;
    restrict = false;
    return;
  }

  my_dt = *j - my_t;
  restrict = true;
}

//! Flush scalar time-series.
void IceModel::flush_timeseries() {
  // flush all the time-series buffers:
  for (std::set<std::string>::iterator j = ts_vars.begin(); j != ts_vars.end(); ++j) {
    TSDiagnostic *diag = ts_diagnostics[*j];

    if (diag != NULL) {
      diag->flush();
    }
  }
}

} // end of namespace pism
