// Copyright (C) 2004-2010 Jed Brown, Ed Bueler and Constantine Khroulev
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
#include <cstdio>
#include <petscda.h>
#include "iceModel.hh"
#include <algorithm>
#include <sstream>
#include <set>
#include "PISMIO.hh"

//! Save model state in NetCDF format.
/*!
Optionally allows saving of full velocity field.

Calls dumpToFile() and writeMatlabVars() to do the actual work.
 */
PetscErrorCode  IceModel::writeFiles(const char* default_filename) {
  PetscErrorCode ierr;
  string filename = default_filename,
    config_out;
  bool o_set, dump_config;

  grid.profiler->begin(event_output);

  ierr = stampHistoryEnd(); CHKERRQ(ierr);

  ierr = PetscOptionsBegin(grid.com, "", "PISM output options", ""); CHKERRQ(ierr);
  {
    ierr = PISMOptionsString("-o", "Output file name", filename, o_set); CHKERRQ(ierr);
    ierr = PISMOptionsString("-dump_config", "File to write the config to",
			     config_out, dump_config); CHKERRQ(ierr);
  }
  ierr = PetscOptionsEnd(); CHKERRQ(ierr);

  if (!ends_with(filename, ".nc")) {
    ierr = verbPrintf(2, grid.com,
		      "PISM WARNING: output file name does not have the '.nc' suffix!\n");
    CHKERRQ(ierr);
  }

  ierr = verbPrintf(2, grid.com, "Writing model state to file `%s'\n", filename.c_str()); CHKERRQ(ierr);


  ierr = dumpToFile(filename.c_str()); CHKERRQ(ierr);


  // save the config file
  if (dump_config) {
    ierr = config.write(config_out.c_str()); CHKERRQ(ierr);
  }

  grid.profiler->end(event_output);

#ifdef PISM_PROFILE
  bool flag;
  ierr = PISMOptionsIsSet("-prof", flag); CHKERRQ(ierr);
  string prof_output_name;
  if (flag) {
    prof_output_name = pism_filename_add_suffix(filename, "-prof", "");

    ierr = verbPrintf(2, grid.com, "Saving profiling data to '%s'...\n",
		      prof_output_name.c_str());
    CHKERRQ(ierr);

    ierr = grid.profiler->save_report(prof_output_name); CHKERRQ(ierr);
  }
#endif

  return 0;
}


PetscErrorCode IceModel::dumpToFile(const char *filename) {
  PetscErrorCode ierr;
  PISMIO nc(&grid);

  // Prepare the file
  ierr = nc.open_for_writing(filename, false, true); CHKERRQ(ierr);
  // append == false, check_dims == true
  ierr = nc.append_time(grid.year); CHKERRQ(ierr);
  ierr = nc.close(); CHKERRQ(ierr);

  // Write metadata *before* variables:

  ierr = mapping.write(filename); CHKERRQ(ierr);
  ierr = global_attributes.write(filename); CHKERRQ(ierr);

  bool override_used;
  ierr = PISMOptionsIsSet("-config_override", override_used); CHKERRQ(ierr);
  if (override_used) {
    overrides.update_from(config);
    ierr = overrides.write(filename); CHKERRQ(ierr);
  }

  if (stress_balance != NULL) {
    ierr = stress_balance->write_model_state(filename); CHKERRQ(ierr);
  } else {
    SETERRQ(1,"PISM ERROR: stress_balance == NULL");
  }

  if (surface != NULL) {
    ierr = surface->write_model_state(grid.year, dt / secpera, filename); CHKERRQ(ierr);
  } else {
    SETERRQ(1,"PISM ERROR: surface == NULL");
  }

  if (ocean != NULL) {
    ierr = ocean->write_model_state(grid.year, dt / secpera, filename); CHKERRQ(ierr);
  } else {
    SETERRQ(1,"PISM ERROR: ocean == NULL");
  }

  ierr = write_model_state(filename);  CHKERRQ(ierr);

  ierr = write_extra_fields(filename); CHKERRQ(ierr); // chance for derived classes to do more

  return 0;
}

//! \brief Writes variables listed in vars to filename, using nctype to write
//! fields stored in dedicated IceModelVecs.
PetscErrorCode IceModel::write_variables(const char *filename, set<string> vars,
					 nc_type nctype) {
  PetscErrorCode ierr;
  IceModelVec *v;

  // FIXME: boundary models should define their variables before writing, too.

  // Ask boundary models to write their variables:
  if (surface != NULL) {
    ierr = surface->write_fields(vars, grid.year, dt / secpera, filename); CHKERRQ(ierr);
  } else {
    SETERRQ(1,"PISM ERROR: surface == NULL");
  }
  if (ocean != NULL) {
    ierr = ocean->write_fields(vars, grid.year, dt / secpera, filename); CHKERRQ(ierr);
  } else {
    SETERRQ(1,"PISM ERROR: ocean == NULL");
  }

  // Define all the variables:
  {
    PISMIO nc(&grid);
    ierr = nc.open_for_writing(filename, true, false); CHKERRQ(ierr);
    // append == true, check_dims == false

    set<string>::iterator i = vars.begin();
    while (i != vars.end()) {
      v = variables.get(*i);

      if (v != NULL) {
        // It has dedicated storage.
        if (*i == "mask") {
          ierr = v->define(nc, NC_BYTE); CHKERRQ(ierr); // use the default data type
        } else {
          ierr = v->define(nc, nctype); CHKERRQ(ierr);
        }
      } else {
        // It might be a diagnostic quantity
        PISMDiagnostic *diag = diagnostics[*i];

        if (diag != NULL) {
          ierr = diag->define(nc); CHKERRQ(ierr);
        }
      }

      ++i;
    }

    ierr = nc.close(); CHKERRQ(ierr);
  }

  // Write all the variables:
  set<string>::iterator i = vars.begin();
  while (i != vars.end()) {
    v = variables.get(*i);

    if (v == NULL) {
      ++i;
    } else {
      ierr = v->write(filename); CHKERRQ(ierr); // use the default data type

      vars.erase(i++);		// note that it only erases variables that were
                                // found (and saved)
    }
  }

  // All the remaining names in vars must be of diagnostic quantities.
  i = vars.begin();
  while (i != vars.end()) {
    PISMDiagnostic *diag = diagnostics[*i];

    if (diag == NULL)
      ++i;
    else {
      v = NULL;

      ierr = diag->compute(v); CHKERRQ(ierr);

      v->write_in_glaciological_units = true;
      ierr = v->write(filename, NC_FLOAT); CHKERRQ(ierr); // diagnostic quantities are always written in float

      delete v;

      vars.erase(i++);
    }
  }

  // check if we have any variables we didn't write
  if (!vars.empty()) {
    int threshold = 3;
    ierr = verbPrintf(threshold, grid.com,
                      "PISM WARNING: the following variables were *not* written by the PISM core: "); CHKERRQ(ierr);
    for (i = vars.begin(); i != vars.end(); ++i) {
      ierr = verbPrintf(threshold, grid.com, "%s, ", (*i).c_str()); CHKERRQ(ierr);
    }
    ierr = verbPrintf(threshold, grid.com, "\b\b\n"); CHKERRQ(ierr);
  }

  return 0;
}


  PetscErrorCode IceModel::write_model_state(const char* filename) {
    PetscErrorCode ierr;

    bool write_temp_pa;
    ierr = PISMOptionsIsSet("-temp_pa", write_temp_pa); CHKERRQ(ierr);
    if (write_temp_pa || (!config.get_flag("do_cold_ice_methods"))) {
      // write temp_pa = pressure-adjusted temp in Celcius
      ierr = verbPrintf(4, grid.com,
                        "  writing pressure-adjusted ice temperature (deg C) 'temp_pa' ...\n"); CHKERRQ(ierr);
      output_vars.insert("temp_pa");
    }

    bool write_liqfrac;
    ierr = PISMOptionsIsSet("-liqfrac", write_liqfrac); CHKERRQ(ierr);
    if (write_liqfrac || (!config.get_flag("do_cold_ice_methods"))) {
      ierr = verbPrintf(4, grid.com,
                        "  writing liquid water fraction 'liqfrac' ...\n"); CHKERRQ(ierr);
      output_vars.insert("liqfrac");
    }

    bool userWantsCTS;
    ierr = PISMOptionsIsSet("-cts", userWantsCTS); CHKERRQ(ierr);
    if (userWantsCTS) {
      ierr = verbPrintf(4, grid.com,
                        "  writing CTS (= E/Es) scalar field 'cts' ...\n"); CHKERRQ(ierr);
      output_vars.insert("cts");
    }

    ierr = write_variables(filename, output_vars, NC_DOUBLE);

    return 0;
  }


  //! Writes extra fields to the output file \c filename. Does nothing in the base class.
  PetscErrorCode IceModel::write_extra_fields(const char* /*filename*/) {
    // Do nothing.
    return 0;
  }


  //! Read a saved PISM model state in NetCDF format, for complete initialization of an evolution or diagnostic run.
  /*!
    Before this is run, the method IceModel::grid_setup() determines the number of
    grid points (Mx,My,Mz,Mbz) and the dimensions (Lx,Ly,Lz) of the computational
    box from the same input file.
  */
  PetscErrorCode IceModel::initFromFile(const char *filename) {
    PetscErrorCode  ierr;
    NCTool nc(grid.com, grid.rank);

    ierr = verbPrintf(2, grid.com, "initializing from NetCDF file '%s'...\n",
                      filename); CHKERRQ(ierr);

    ierr = nc.open_for_reading(filename); CHKERRQ(ierr);

    // Find the index of the last record in the file:
    int last_record;
    ierr = nc.get_dim_length("t", &last_record); CHKERRQ(ierr);
    last_record -= 1;

    // Read the model state, mapping and climate_steady variables:
    set<string> vars = variables.keys();

    set<string>::iterator i = vars.begin();
    while (i != vars.end()) {
      IceModelVec *var = variables.get(*i++);

      string intent = var->string_attr("pism_intent");
      if ((intent == "model_state") || (intent == "mapping") ||
          (intent == "climate_steady")) {
        ierr = var->read(filename, last_record); CHKERRQ(ierr);
      }
    }

    ierr = verbPrintf(3,grid.com,"Setting enthalpy from temperature...\n"); CHKERRQ(ierr);
    if (config.get_flag("do_cold_ice_methods")) {
      ierr = setEnth3FromT3_ColdIce(); CHKERRQ(ierr);
    }


    if (config.get_flag("do_age")) {
      bool age_exists;
      ierr = nc.find_variable("age", NULL, age_exists); CHKERRQ(ierr);

      if (age_exists) {
        ierr = tau3.read(filename, last_record); CHKERRQ(ierr);
      } else {
        ierr = verbPrintf(2,grid.com,
                          "PISM WARNING: input file '%s' does not have the 'age' variable.\n"
                          "  Setting it to zero...\n",
                          filename); CHKERRQ(ierr);
        ierr = tau3.set(0.0); CHKERRQ(ierr);
      }
    }

    // options for initializing enthalpy
    bool initfromT, initfromTandOm;
    ierr = PISMOptionsIsSet("-init_from_temp", initfromT); CHKERRQ(ierr);
    ierr = PISMOptionsIsSet("-init_from_temp_and_liqfrac", initfromTandOm); CHKERRQ(ierr);
    if (initfromT && initfromTandOm) {
      SETERRQ(1,"PISM ERROR: both options -init_from_temp and -init_from_temp_and_liqfrac seen; must choose one or none");
    }
    if (initfromT) {
      ierr = verbPrintf(2, grid.com,
                        "  option -init_from_temp seen; computing enthalpy from ice temperature and thickness ...\n"); CHKERRQ(ierr);
      ierr = setEnth3FromT3_ColdIce(); CHKERRQ(ierr);
    } else if (initfromTandOm) {
      ierr = verbPrintf(2, grid.com,
                        "  option -init_from_temp_and_liqfrac seen; computing enthalpy from ice temperature, liquid water fraction and thickness ...\n"); CHKERRQ(ierr);
      // use vWork3d as already-allocated space
      ierr = vWork3d.set_name("liqfrac"); CHKERRQ(ierr);
      ierr = vWork3d.set_attrs(
                               "internal", "liqfrac; temporary use during initialization",
                               "", ""); CHKERRQ(ierr);
      ierr = vWork3d.read(filename,last_record); CHKERRQ(ierr);
      ierr = setEnth3FromT3AndLiqfrac3(vWork3d); CHKERRQ(ierr);
    }

    string history;
    ierr = nc.get_att_text(NC_GLOBAL, "history", history); CHKERRQ(ierr);
    global_attributes.prepend_history(history);

    ierr = nc.close(); CHKERRQ(ierr);

    return 0;
  }

  //! Manage regridding based on user options.  Call IceModelVec::regrid() to do each selected variable.
  /*!
    For each variable selected by option <tt>-regrid_vars</tt>, we regrid it onto the current grid from
    the NetCDF file specified by <tt>-regrid_file</tt>.

    The default, if <tt>-regrid_vars</tt> is not given, is to regrid the 3
    dimensional quantities \c tau3, \c Tb3 and either \c T3 or \c Enth3. This is
    consistent with one standard purpose of regridding, which is to stick with
    current geometry through the downscaling procedure. Most of the time the user
    should carefully specify which variables to regrid.
  */
  PetscErrorCode IceModel::regrid() {
    PetscErrorCode ierr;
    string filename, tmp;
    bool regridVarsSet, regrid_file_set;
    PISMIO nc(&grid);

    ierr = PetscOptionsBegin(grid.com, PETSC_NULL,
                             "Options controlling regridding",
                             PETSC_NULL); CHKERRQ(ierr);

    // Get the regridding file name:
    ierr = PISMOptionsString("-regrid_file", "Specifies the file to regrid from",
                             filename, regrid_file_set); CHKERRQ(ierr);

    ierr = PISMOptionsString("-regrid_vars", "Specifies the list of variable to regrid",
                             tmp, regridVarsSet); CHKERRQ(ierr);

    // Done with the options.
    ierr = PetscOptionsEnd(); CHKERRQ(ierr);

    // Return if no regridding is requested:
    if (!regrid_file_set) return 0;

    ierr = verbPrintf(2, grid.com, "regridding from file %s ...\n",filename.c_str()); CHKERRQ(ierr);

    string var_name;
    set<string> vars;
    istringstream list(tmp);
    if (regridVarsSet) {
      // split the list; note that this also removes any duplicate entries
      while (getline(list, var_name, ','))
        vars.insert(var_name);
    } else {
      vars.insert("litho_temp");

      if (config.get_flag("do_age"))
	vars.insert("age");

      if (config.get_flag("do_cold_ice_methods"))
        vars.insert("temp");
      else
        vars.insert("enthalpy");
    }

    // create "local interpolation context" from dimensions, limits, and lengths
    //   extracted from regridFile, and from information about the part of the
    //   grid owned by this processor

    ierr = nc.open_for_reading(filename.c_str());

    grid_info g;
    // Note that after this call g.z_len and g.zb_len are zero if the
    // corresponding dimension does not exist.
    ierr = nc.get_grid_info(g); CHKERRQ(ierr);

    double *zlevs = NULL, *zblevs = NULL; // NULLs correspond to 2D-only regridding
    if ((g.z_len != 0) && (g.zb_len != 0)) {
      ierr = nc.get_vertical_dims(zlevs, zblevs); CHKERRQ(ierr);
    } else {
      ierr = verbPrintf(2, grid.com,
                        "PISM WARNING: at least one of 'z' and 'zb' is absent in '%s'.\n"
                        "              3D regridding is disabled.\n",
                        filename.c_str());
      CHKERRQ(ierr);
    }
    ierr = nc.close(); CHKERRQ(ierr);

    LocalInterpCtx lic(g, zlevs, zblevs, grid); // will be de-allocated at 'return 0' below.

    set<string>::iterator i;
    for (i = vars.begin(); i != vars.end(); ++i) {
      IceModelVec *v = variables.get(*i);

      if (v == NULL) {
        ierr = PetscPrintf(grid.com, "PISM ERROR: unknown variable name: %s\n",
                           (*i).c_str()); CHKERRQ(ierr);
        PetscEnd();
      }

      string pism_intent = v->string_attr("pism_intent");
      if (pism_intent != "model_state") {
        ierr = verbPrintf(2, grid.com, "  WARNING: skipping '%s' (only model_state variables can be regridded)...\n",
                          (*i).c_str()); CHKERRQ(ierr);
        continue;
      }

      if ( ((v->grid_type() == GRID_3D) && lic.no_regrid_ice) ||
           ((v->grid_type() == GRID_3D_BEDROCK) && lic.no_regrid_bedrock) )
        {
          ierr = verbPrintf(2, grid.com, "  WARNING: skipping '%s'...\n",
                            (*i).c_str()); CHKERRQ(ierr);
          continue;
        }

      ierr = v->regrid(filename.c_str(), lic, true); CHKERRQ(ierr);

    }

    // Note that deleting a NULL pointer is safe.
    delete [] zlevs;  delete [] zblevs;
    return 0;
  }

  //! Initializes the snapshot-saving mechanism.
  PetscErrorCode IceModel::init_snapshots() {
    PetscErrorCode ierr;
    bool save_at_set, save_to_set, split;
    string tmp;
    current_snapshot = 0;

    ierr = PetscOptionsBegin(grid.com, "", "Options controlling the snapshot-saving mechanism", ""); CHKERRQ(ierr);
    {
      ierr = PISMOptionsString("-save_file", "Specifies a snapshot filename",
                               snapshots_filename, save_to_set); CHKERRQ(ierr);

      ierr = PISMOptionsString("-save_times", "Gives a list or a MATLAB-style range of times to save snapshots at",
                               tmp, save_at_set); CHKERRQ(ierr);

      ierr = PISMOptionsIsSet("-save_split", "Specifies whether to save snapshots to separate files",
                              split); CHKERRQ(ierr);

      ierr = set_output_size("-save_size", "Sets the 'size' of a snapshot file.",
                             "small", snapshot_vars); CHKERRQ(ierr);
    }
    ierr = PetscOptionsEnd(); CHKERRQ(ierr);

    if (save_to_set ^ save_at_set) {
      ierr = PetscPrintf(grid.com,
                         "PISM ERROR: you need to specify both -save_file and -save_times to save snapshots.\n");
      CHKERRQ(ierr);
      PetscEnd();
    }

    if (!save_to_set && !save_at_set) {
      save_snapshots = false;
      return 0;
    }

    ierr = parse_times(grid.com, tmp, snapshot_times);
    if (ierr != 0) {
      ierr = PetscPrintf(grid.com, "PISM ERROR: parsing the -save_times argument failed.\n"); CHKERRQ(ierr);
      PetscEnd();
    }

    if (snapshot_times.size() == 0) {
      PetscPrintf(grid.com, "PISM ERROR: no argument for -save_times option.\n");
      PetscEnd();
    }

    save_snapshots = true;
    snapshots_file_is_ready = false;
    split_snapshots = false;

    if (split) {
      split_snapshots = true;
    } else if (!ends_with(snapshots_filename, ".nc")) {
      ierr = verbPrintf(2, grid.com,
                        "PISM WARNING: snapshots file name does not have the '.nc' suffix!\n");
      CHKERRQ(ierr);
    }

    if (split) {
      ierr = verbPrintf(2, grid.com, "saving snapshots to '%s+year.nc'; ",
                        snapshots_filename.c_str()); CHKERRQ(ierr);
    } else {
      ierr = verbPrintf(2, grid.com, "saving snapshots to '%s'; ",
                        snapshots_filename.c_str()); CHKERRQ(ierr);
    }

    ierr = verbPrintf(2, grid.com, "times requested: %s\n", tmp.c_str()); CHKERRQ(ierr);

    return 0;
  }

  //! Writes a snapshot of the model state (if necessary)
  PetscErrorCode IceModel::write_snapshot() {
    PetscErrorCode ierr;
    PISMIO nc(&grid);
    double saving_after = -1.0e30; // initialize to avoid compiler warning; this
    // value is never used, because saving_after
    // is only used if save_now == true, and in
    // this case saving_after is guaranteed to be
    // initialized. See the code below.
    char filename[PETSC_MAX_PATH_LEN];

    // determine if the user set the -save_times and -save_file options
    if (!save_snapshots)
      return 0;

    // do we need to save *now*?
    if ( (grid.year >= snapshot_times[current_snapshot]) && (current_snapshot < snapshot_times.size()) ) {
      saving_after = snapshot_times[current_snapshot];

      while ((current_snapshot < snapshot_times.size()) &&
             (snapshot_times[current_snapshot] <= grid.year))
        current_snapshot++;
    } else {
      // we don't need to save now, so just return
      return 0;
    }

    // flush time-series buffers
    ierr = flush_timeseries(); CHKERRQ(ierr);

    if (split_snapshots) {
      snapshots_file_is_ready = false;	// each snapshot is written to a separate file
      snprintf(filename, PETSC_MAX_PATH_LEN, "%s-%06.0f.nc",
               snapshots_filename.c_str(), grid.year);
    } else {
      strncpy(filename, snapshots_filename.c_str(), PETSC_MAX_PATH_LEN);
    }

    ierr = verbPrintf(2, grid.com,
                      "\nsaving snapshot to %s at %.5f a, for time-step goal %.5f a\n\n",
                      filename, grid.year,saving_after);
    CHKERRQ(ierr);

    // create line for history in .nc file, including time of write

    string date_str = pism_timestamp();
    char tmp[TEMPORARY_STRING_LENGTH];
    snprintf(tmp, TEMPORARY_STRING_LENGTH,
             "%s: %s snapshot at %10.5f a, for time-step goal %10.5f a\n",
             date_str.c_str(), executable_short_name.c_str(), grid.year, saving_after);

    if (!snapshots_file_is_ready) {

      // Prepare the snapshots file:
      ierr = nc.open_for_writing(filename, false, true); CHKERRQ(ierr);
      // append == false, check_dims == true
      ierr = nc.close(); CHKERRQ(ierr);

      ierr = global_attributes.write(filename); CHKERRQ(ierr);
      ierr = mapping.write(filename); CHKERRQ(ierr);
      snapshots_file_is_ready = true;
    }

    ierr = nc.open_for_writing(filename, true, true); CHKERRQ(ierr);
    // append == true, check_dims == true
    ierr = nc.append_time(grid.year); CHKERRQ(ierr);
    ierr = nc.write_history(tmp); CHKERRQ(ierr); // append the history
    ierr = nc.close(); CHKERRQ(ierr);

    // Let boundary models write their fields:
    if (surface != NULL) {
      ierr = surface->write_model_state(grid.year, dt / secpera, filename); CHKERRQ(ierr);
    } else {
      SETERRQ(1,"PISM ERROR: surface == NULL");
    }

    if (ocean != NULL) {
      ierr = ocean->write_model_state(grid.year, dt / secpera, filename); CHKERRQ(ierr);
    } else {
      SETERRQ(1,"PISM ERROR: ocean == NULL");
    }

    ierr = write_variables(filename, snapshot_vars, NC_DOUBLE);

    ierr = write_extra_fields(filename); CHKERRQ(ierr);

    return 0;
  }

  //! Initialize the backup (snapshot-on-wallclock-time) mechanism.
  PetscErrorCode IceModel::init_backups() {
    PetscErrorCode ierr;
    bool o_set;

    backup_interval = config.get("backup_interval");

    ierr = PetscOptionsBegin(grid.com, "", "PISM output options", ""); CHKERRQ(ierr);
    {
      ierr = PISMOptionsString("-o", "Output file name", backup_filename, o_set); CHKERRQ(ierr);
      if (!o_set)
        backup_filename = executable_short_name + "_backup.nc";
      else
        backup_filename = pism_filename_add_suffix(backup_filename, "_backup", "");

      ierr = PISMOptionsReal("-backup_interval", "Automatic backup interval, hours",
                             backup_interval, o_set); CHKERRQ(ierr);
    }
    ierr = PetscOptionsEnd(); CHKERRQ(ierr);

    last_backup_time = 0.0;

    return 0;
  }

  //! Write a backup (i.e. an intermediate result of a run).
PetscErrorCode IceModel::write_backup() {
  PetscErrorCode ierr;
  double wall_clock_hours;
  PISMIO nc(&grid);

  if (grid.rank == 0) {
    PetscLogDouble current_time;
    ierr = PetscGetTime(&current_time); CHKERRQ(ierr);
    wall_clock_hours = (current_time - start_time) / 3600.0;
  }

  MPI_Bcast(&wall_clock_hours, 1, MPI_DOUBLE, 0, grid.com);

  if (wall_clock_hours - last_backup_time < backup_interval)
    return 0;

  last_backup_time = wall_clock_hours;

  // create a history string:
  string date_str = pism_timestamp();
  char tmp[TEMPORARY_STRING_LENGTH];
  snprintf(tmp, TEMPORARY_STRING_LENGTH,
           "%s automatic backup at %10.5f a, %3.3f hours after the beginning of the run\n",
           executable_short_name.c_str(), grid.year, wall_clock_hours);

  ierr = verbPrintf(2, grid.com,
                    "  Saving an automatic backup to '%s' (%1.3f hours after the beginning of the run)\n",
                    backup_filename.c_str(), wall_clock_hours); CHKERRQ(ierr);

  stampHistory(tmp);

  // write metadata:
  ierr = nc.open_for_writing(backup_filename.c_str(), false, true); CHKERRQ(ierr);
  // append == false, check_dims == true
  ierr = nc.append_time(grid.year); CHKERRQ(ierr);
  ierr = nc.close(); CHKERRQ(ierr);

  ierr = global_attributes.write(backup_filename.c_str()); CHKERRQ(ierr);
  ierr = mapping.write(backup_filename.c_str()); CHKERRQ(ierr);

  // write the model state (this saves only the fields necessary for restarting
  // and *does not* respect -o_size)
  set<string> vars = variables.keys();

  set<string>::iterator i = vars.begin();
  while (i != vars.end()) {
    IceModelVec *var = variables.get(*i++);

    string intent = var->string_attr("pism_intent");
    if ((intent == "model_state") || (intent == "mapping") ||
        (intent == "climate_steady")) {
      ierr = var->write(backup_filename.c_str()); CHKERRQ(ierr);
    }
  }

  if (stress_balance != NULL) {
    ierr = stress_balance->write_model_state(backup_filename); CHKERRQ(ierr);
  } else {
    SETERRQ(1,"PISM ERROR: stress_balance == NULL");
  }

  if (surface != NULL) {
    ierr = surface->write_model_state(grid.year, dt / secpera, backup_filename); CHKERRQ(ierr);
  } else {
    SETERRQ(1,"PISM ERROR: surface == NULL");
  }

  if (ocean != NULL) {
    ierr = ocean->write_model_state(grid.year, dt / secpera, backup_filename); CHKERRQ(ierr);
  } else {
    SETERRQ(1,"PISM ERROR: ocean == NULL");
  }

  // Also flush time-series:
  ierr = flush_timeseries(); CHKERRQ(ierr);

  return 0;
}
