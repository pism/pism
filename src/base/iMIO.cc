// Copyright (C) 2004-2012 Jed Brown, Ed Bueler and Constantine Khroulev
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
#include <petscdmda.h>
#include "iceModel.hh"
#include <algorithm>
#include <sstream>
#include <set>

#include "PIO.hh"
#include "PISMBedDef.hh"
#include "bedrockThermalUnit.hh"
#include "PISMYieldStress.hh"
#include "PISMStressBalance.hh"
#include "PISMSurface.hh"
#include "PISMOcean.hh"
#include "PISMProf.hh"
#include "pism_options.hh"
#include "IceGrid.hh"
#include "PISMTime.hh"
#include "PISMDiagnostic.hh"

//! Save model state in NetCDF format.
/*!
Optionally allows saving of full velocity field.

Calls dumpToFile() to do the actual work.
 */
PetscErrorCode  IceModel::writeFiles(string default_filename) {
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

  if (get_output_size("-o_size") != "none") {
    ierr = verbPrintf(2, grid.com,
                      "Writing model state to file `%s'\n", filename.c_str()); CHKERRQ(ierr);
    ierr = dumpToFile(filename); CHKERRQ(ierr);
  }

  // save the config file
  if (dump_config) {
    ierr = config.write(config_out); CHKERRQ(ierr);
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

//! \brief Write metadata (global attributes, overrides and mapping parameters) to a file.
PetscErrorCode IceModel::write_metadata(string filename, bool write_mapping) {
  PetscErrorCode ierr;

  if (write_mapping) {
    ierr = mapping.write(filename); CHKERRQ(ierr);
  }

  ierr = global_attributes.write(filename); CHKERRQ(ierr);

  bool override_used;
  ierr = PISMOptionsIsSet("-config_override", override_used); CHKERRQ(ierr);
  if (override_used) {
    overrides.update_from(config);
    ierr = overrides.write(filename); CHKERRQ(ierr);
  }

  return 0;
}


PetscErrorCode IceModel::dumpToFile(string filename) {
  PetscErrorCode ierr;
  PIO nc(grid.com, grid.rank, grid.config.get_string("output_format"));

  // Prepare the file
  string time_name = config.get_string("time_dimension_name");
  ierr = nc.open(filename, PISM_WRITE); CHKERRQ(ierr);
  ierr = nc.def_time(time_name, config.get_string("calendar"), grid.time->CF_units()); CHKERRQ(ierr);
  ierr = nc.append_time(time_name, grid.time->current()); CHKERRQ(ierr);
  ierr = nc.close(); CHKERRQ(ierr);

  // Write metadata *before* variables:

  ierr = write_metadata(filename); CHKERRQ(ierr);

  ierr = write_model_state(filename);  CHKERRQ(ierr);

  return 0;
}

//! \brief Writes variables listed in vars to filename, using nctype to write
//! fields stored in dedicated IceModelVecs.
PetscErrorCode IceModel::write_variables(string filename, set<string> vars,
					 PISM_IO_Type nctype) {
  PetscErrorCode ierr;
  IceModelVec *v;

  grid.profiler->begin(event_output_define);

  // Define all the variables:
  {
    string output_format = grid.config.get_string("output_format");

    // This is a kludge: for some reason defining variables using PnetCDF takes
    // a very, very long time. It uses the NetCDF-3 file format though, so we
    // *override* this setting here and a) define variables using NetCDF-3 and
    // then b) write data using PnetCDF.
    //
    // I suspect that dimension and variable lookup is to blame: PISM does not
    // store dimension/variable IDs and looks them up every time they are
    // needed. A simple test executable (pism_netcdf_test) does not have this
    // issue.
    //
    // Note: variable metadata has nothing to do with this -- increasing header
    // padding does not help.
    if (output_format == "pnetcdf")
      output_format = "netcdf3";

    PIO nc(grid.com, grid.rank, output_format);
    ierr = nc.open(filename, PISM_WRITE, true); CHKERRQ(ierr);

    set<string>::iterator i = vars.begin();
    while (i != vars.end()) {
      v = variables.get(*i);

      if (v != NULL) {
        // It has dedicated storage.
        if (*i == "mask") {
          ierr = v->define(nc, PISM_BYTE); CHKERRQ(ierr); // use the default data type
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

    if (beddef != NULL) {
      ierr = beddef->define_variables(vars, nc, nctype); CHKERRQ(ierr);
    }

    if (btu != NULL) {
      ierr = btu->define_variables(vars, nc, nctype); CHKERRQ(ierr);
    }

    if (basal_yield_stress != NULL) {
      ierr = basal_yield_stress->define_variables(vars, nc, nctype); CHKERRQ(ierr);
    }

    if (stress_balance != NULL) {
      ierr = stress_balance->define_variables(vars, nc, nctype); CHKERRQ(ierr);
    } else {
      SETERRQ(grid.com, 1,"PISM ERROR: stress_balance == NULL");
    }

    if (surface != NULL) {
      ierr = surface->define_variables(vars, nc, nctype); CHKERRQ(ierr);
    } else {
      SETERRQ(grid.com, 1,"PISM ERROR: surface == NULL");
    }
    if (ocean != NULL) {
      ierr = ocean->define_variables(vars, nc, nctype); CHKERRQ(ierr);
    } else {
      SETERRQ(grid.com, 1,"PISM ERROR: ocean == NULL");
    }

    ierr = nc.close(); CHKERRQ(ierr);
  }

  grid.profiler->end(event_output_define);

  // Write all the IceModel variables:
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

  // Write bed-deformation-related variables:
  if (beddef != NULL) {
    ierr = beddef->write_variables(vars, filename); CHKERRQ(ierr);
  }

  // Write PISMBedThermalUnit variables:
  if (btu != NULL) {
    ierr = btu->write_variables(vars, filename); CHKERRQ(ierr);
  }

  if (basal_yield_stress != NULL) {
    ierr = basal_yield_stress->write_variables(vars, filename); CHKERRQ(ierr);
  }

  // Write stress balance-related variables:
  if (stress_balance != NULL) {
    ierr = stress_balance->write_variables(vars, filename); CHKERRQ(ierr);
  } else {
    SETERRQ(grid.com, 1,"PISM ERROR: stress_balance == NULL");
  }

  // Ask boundary models to write their variables:
  if (surface != NULL) {
    ierr = surface->write_variables(vars, filename); CHKERRQ(ierr);
  } else {
    SETERRQ(grid.com, 1,"PISM ERROR: surface == NULL");
  }
  if (ocean != NULL) {
    ierr = ocean->write_variables(vars, filename); CHKERRQ(ierr);
  } else {
    SETERRQ(grid.com, 1,"PISM ERROR: ocean == NULL");
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
      ierr = v->write(filename, PISM_FLOAT); CHKERRQ(ierr); // diagnostic quantities are always written in float

      delete v;

      vars.erase(i++);
    }
  }

  // FIXME: we need a way of figuring out if a sub-model did or did not write
  // something.

  if (!vars.empty()) {
    int threshold = 3;
    ierr = verbPrintf(threshold, grid.com,
                      "PISM WARNING: the following variables were *not* written by PISM core (IceModel): "); CHKERRQ(ierr);
    for (i = vars.begin(); i != vars.end(); ++i) {
      ierr = verbPrintf(threshold, grid.com, "%s ", (*i).c_str()); CHKERRQ(ierr);
    }
    ierr = verbPrintf(threshold, grid.com, "\n"); CHKERRQ(ierr);
  }

  return 0;
}


PetscErrorCode IceModel::write_model_state(string filename) {
  PetscErrorCode ierr;
  string o_size = get_output_size("-o_size");

  // only write out these extra diagnostics if decently big
  if (o_size == "medium" || o_size == "big") {
    bool write_temp_pa, write_liqfrac;
    ierr = PISMOptionsIsSet("-temp_pa", write_temp_pa); CHKERRQ(ierr);
    if (write_temp_pa || (!config.get_flag("do_cold_ice_methods"))) {
      // write temp_pa = pressure-adjusted temp in Celcius
      ierr = verbPrintf(4, grid.com,
                        "  writing pressure-adjusted ice temperature (deg C) 'temp_pa' ...\n");
                        CHKERRQ(ierr);
      output_vars.insert("temp_pa");
    }
    ierr = PISMOptionsIsSet("-liqfrac", write_liqfrac); CHKERRQ(ierr);
    if (write_liqfrac || (!config.get_flag("do_cold_ice_methods"))) {
      ierr = verbPrintf(4, grid.com,
                        "  writing liquid water fraction 'liqfrac' ...\n");
                        CHKERRQ(ierr);
      output_vars.insert("liqfrac");
    }
  }

  // if user wants it, give it to them (ignor -o_size, except "none")
  bool userWantsCTS;
  ierr = PISMOptionsIsSet("-cts", userWantsCTS); CHKERRQ(ierr);
  if (userWantsCTS) {
    ierr = verbPrintf(4, grid.com,
                      "  writing CTS (= E/Es) scalar field 'cts' ...\n"); CHKERRQ(ierr);
    output_vars.insert("cts");
  }

  ierr = write_variables(filename, output_vars, PISM_DOUBLE);

  return 0;
}



  //! Read a saved PISM model state in NetCDF format, for complete initialization of an evolution or diagnostic run.
  /*!
    Before this is run, the method IceModel::grid_setup() determines the number of
    grid points (Mx,My,Mz,Mbz) and the dimensions (Lx,Ly,Lz) of the computational
    box from the same input file.
  */
PetscErrorCode IceModel::initFromFile(string filename) {
  PetscErrorCode  ierr;
  PIO nc(grid.com, grid.rank, grid.config.get_string("output_format"));

  ierr = verbPrintf(2, grid.com, "initializing from NetCDF file '%s'...\n",
                    filename.c_str()); CHKERRQ(ierr);

  // options for initializing enthalpy
  bool initfromT, initfromTandOm;
  string enthalpy_pism_intent;
  ierr = PISMOptionsIsSet("-init_from_temp", initfromT); CHKERRQ(ierr);
  ierr = PISMOptionsIsSet("-init_from_temp_and_liqfrac", initfromTandOm); CHKERRQ(ierr);
  if (initfromT && initfromTandOm) {
    PetscPrintf(grid.com,
                "PISM ERROR: both options -init_from_temp and -init_from_temp_and_liqfrac seen; must choose one or none");
    PISMEnd();
  }

  // if initializing enthalpy from temperature or temperature and liquid water
  // fraction, temporarily set pism_intent of Enth3 to "diagnostic" so that
  // PISM does not even try to read the "enthalpy" variable
  if (initfromT || initfromTandOm) {
    enthalpy_pism_intent = Enth3.string_attr("pism_intent");
    ierr = Enth3.set_attr("pism_intent", "diagnostic"); CHKERRQ(ierr);
  }


  ierr = nc.open(filename, PISM_NOWRITE); CHKERRQ(ierr);

  // check if the input file has Href; set its pism_intent to "diagnostic" and
  // set the field itself to 0 if it is not present
  if (config.get_flag("part_grid")) {
    bool exists;
    ierr = nc.inq_var("Href", exists); CHKERRQ(ierr);

    if (!exists) {
      ierr = verbPrintf(2, grid.com,
        "PISM WARNING: Href for PISM-PIK -part_grid not found in '%s'. Setting it to zero...\n",
                        filename.c_str()); CHKERRQ(ierr);

      ierr = vHref.set_attr("pism_intent", "diagnostic"); CHKERRQ(ierr);
      ierr = vHref.set(0.0); CHKERRQ(ierr);
    }
  }

  // Find the index of the last record in the file:
  unsigned int last_record;
  ierr = nc.inq_nrecords(last_record); CHKERRQ(ierr); 
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

  if (config.get_flag("do_cold_ice_methods")) {
    ierr = verbPrintf(3,grid.com,"  setting enthalpy from temperature...\n");
    CHKERRQ(ierr);
    ierr = compute_enthalpy_cold(T3, Enth3); CHKERRQ(ierr);
  }

  if (config.get_flag("do_age")) {
    bool age_exists;
    ierr = nc.inq_var("age", age_exists); CHKERRQ(ierr);

    if (age_exists) {
      ierr = tau3.read(filename, last_record); CHKERRQ(ierr);
    } else {
      ierr = verbPrintf(2,grid.com,
                        "PISM WARNING: input file '%s' does not have the 'age' variable.\n"
                        "  Setting it to zero...\n",
                        filename.c_str()); CHKERRQ(ierr);
      ierr = tau3.set(0.0); CHKERRQ(ierr);
    }
  }

  if (initfromT || initfromTandOm) {
    IceModelVec3 temperature;

    ierr = temperature.create(grid, "temp", false); CHKERRQ(ierr);
    ierr = temperature.set_attrs("internal", "ice temperature; temporary storage during initialization",
                                 "K", "land_ice_temperature"); CHKERRQ(ierr);
    ierr = temperature.set_attr("valid_min", 0.0); CHKERRQ(ierr);

    ierr = temperature.read(filename, last_record); CHKERRQ(ierr);

    if (initfromT) {
      ierr = verbPrintf(2, grid.com,
                        "  option -init_from_temp seen;"
                        " computing enthalpy from ice temperature and thickness ...\n");
      CHKERRQ(ierr);

      ierr = compute_enthalpy_cold(temperature, Enth3); CHKERRQ(ierr);

    } else {
      ierr = verbPrintf(2, grid.com,
                        "  option -init_from_temp_and_liqfrac seen; computing enthalpy from ice temperature,"
                        " liquid water fraction and thickness ...\n"); CHKERRQ(ierr);

      // use vWork3d as already-allocated space
      ierr = vWork3d.set_name("liqfrac"); CHKERRQ(ierr);
      ierr = vWork3d.set_attrs("internal", "liqfrac; temporary use during initialization",
                               "", ""); CHKERRQ(ierr);
      ierr = vWork3d.read(filename, last_record); CHKERRQ(ierr);

      ierr = compute_enthalpy(temperature, vWork3d, Enth3); CHKERRQ(ierr);
    }

    // restore the original value of pism_intent (currently: either
    // "model_state" or "diagnostic" depending on the "do_cold_ice_methods"
    // config parameter).
    ierr = Enth3.set_attr("pism_intent", enthalpy_pism_intent); CHKERRQ(ierr);
  }

  // re-set Href's pism_intent attribute
  if (config.get_flag("part_grid")) {
    ierr = vHref.set_attr("pism_intent", "model_state"); CHKERRQ(ierr);
  }

  string history;
  ierr = nc.get_att_text("PISM_GLOBAL", "history", history); CHKERRQ(ierr);
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

    This \c dimensions argument can be 2 (regrid 2D variables only), 3 (3D
    only) and 0 (everything).
  */
  PetscErrorCode IceModel::regrid(int dimensions) {
    PetscErrorCode ierr;
    string filename;
    bool regrid_vars_set, regrid_file_set;
    vector<string> vars_vector;

    if (! (dimensions == 0 ||
           dimensions == 2 ||
           dimensions == 3))
      SETERRQ(grid.com, 1, "dimensions can only be 0, 2 or 3");

    ierr = PetscOptionsBegin(grid.com, PETSC_NULL, "Options controlling regridding",
                             PETSC_NULL); CHKERRQ(ierr);
    {
      ierr = PISMOptionsString("-regrid_file", "Specifies the file to regrid from",
                               filename, regrid_file_set); CHKERRQ(ierr);

      ierr = PISMOptionsStringArray("-regrid_vars", "Specifies the list of variables to regrid",
                                    "", vars_vector, regrid_vars_set); CHKERRQ(ierr);
    }
    ierr = PetscOptionsEnd(); CHKERRQ(ierr);

    // Return if no regridding is requested:
    if (!regrid_file_set) return 0;

    if (dimensions != 0) {
      ierr = verbPrintf(2, grid.com, "regridding %dD variables from file %s ...\n",
                        dimensions, filename.c_str()); CHKERRQ(ierr);
    } else {
      ierr = verbPrintf(2, grid.com, "regridding from file %s ...\n",filename.c_str()); CHKERRQ(ierr);
    }

    set<string> vars;
    vector<string>::iterator j;
    for (j = vars_vector.begin(); j != vars_vector.end(); ++j)
      vars.insert(*j);

    if (vars.empty()) {
      // defaults if user gives no regrid_vars list
      vars.insert("litho_temp");

      if (config.get_flag("do_age"))
	vars.insert("age");

      if (config.get_flag("do_cold_ice_methods"))
        vars.insert("temp");
      else
        vars.insert("enthalpy");
    }

    if (dimensions == 0) {
      ierr = regrid_variables(filename, vars, 2); CHKERRQ(ierr); 
      ierr = regrid_variables(filename, vars, 3); CHKERRQ(ierr); 
    } else {
      ierr = regrid_variables(filename, vars, dimensions); CHKERRQ(ierr); 
    }

    return 0;
  }

PetscErrorCode IceModel::regrid_variables(string filename, set<string> vars, int ndims) {
  PetscErrorCode ierr;

  set<string>::iterator i;
  for (i = vars.begin(); i != vars.end(); ++i) {
    IceModelVec *v = variables.get(*i);

    if (v == NULL) continue;

    if (v->get_ndims() != ndims) continue;

    string pism_intent = v->string_attr("pism_intent");
    if (pism_intent != "model_state") {
      ierr = verbPrintf(2, grid.com, "  WARNING: skipping '%s' (only model_state variables can be regridded)...\n",
                        (*i).c_str()); CHKERRQ(ierr);
      continue;
    }

    ierr = v->regrid(filename, true); CHKERRQ(ierr);
  }

  return 0;
}


//! Initializes the snapshot-saving mechanism.
PetscErrorCode IceModel::init_snapshots() {
  PetscErrorCode ierr;
  bool save_times_set, save_file_set, split;
  string tmp;
  current_snapshot = 0;

  ierr = PetscOptionsBegin(grid.com, "", "Options controlling the snapshot-saving mechanism", ""); CHKERRQ(ierr);
  {
    ierr = PISMOptionsString("-save_file", "Specifies a snapshot filename",
                             snapshots_filename, save_file_set); CHKERRQ(ierr);

    ierr = PISMOptionsString("-save_times", "Gives a list or a MATLAB-style range of times to save snapshots at",
                             tmp, save_times_set); CHKERRQ(ierr);

    ierr = PISMOptionsIsSet("-save_split", "Specifies whether to save snapshots to separate files",
                            split); CHKERRQ(ierr);

    ierr = set_output_size("-save_size", "Sets the 'size' of a snapshot file.",
                           "small", snapshot_vars); CHKERRQ(ierr);
  }
  ierr = PetscOptionsEnd(); CHKERRQ(ierr);

  if (save_file_set ^ save_times_set) {
    ierr = PetscPrintf(grid.com,
                       "PISM ERROR: you need to specify both -save_file and -save_times to save snapshots.\n");
    CHKERRQ(ierr);
    PISMEnd();
  }

  if (!save_file_set && !save_times_set) {
    save_snapshots = false;
    return 0;
  }

  ierr = parse_times(grid.com, config, tmp, snapshot_times);
  if (ierr != 0) {
    ierr = PetscPrintf(grid.com, "PISM ERROR: parsing the -save_times argument failed.\n"); CHKERRQ(ierr);
    PISMEnd();
  }

  if (snapshot_times.size() == 0) {
    PetscPrintf(grid.com, "PISM ERROR: no argument for -save_times option.\n");
    PISMEnd();
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
    PIO nc(grid.com, grid.rank, grid.config.get_string("output_format"));
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
    if ( (grid.time->current() >= snapshot_times[current_snapshot]) && (current_snapshot < snapshot_times.size()) ) {
      saving_after = snapshot_times[current_snapshot];

      while ((current_snapshot < snapshot_times.size()) &&
             (snapshot_times[current_snapshot] <= grid.time->current()))
        current_snapshot++;
    } else {
      // we don't need to save now, so just return
      return 0;
    }

    grid.profiler->begin(event_snapshots);

    // flush time-series buffers
    ierr = flush_timeseries(); CHKERRQ(ierr);

    if (split_snapshots) {
      snapshots_file_is_ready = false;	// each snapshot is written to a separate file
      snprintf(filename, PETSC_MAX_PATH_LEN, "%s-%s.nc",
               snapshots_filename.c_str(), grid.time->date().c_str());
    } else {
      strncpy(filename, snapshots_filename.c_str(), PETSC_MAX_PATH_LEN);
    }

    ierr = verbPrintf(2, grid.com,
                      "\nsaving snapshot to %s at %s, for time-step goal %s\n\n",
                      filename, grid.time->date().c_str(),
                      grid.time->date(saving_after).c_str());
    CHKERRQ(ierr);

    // create line for history in .nc file, including time of write

    string date_str = pism_timestamp();
    char tmp[TEMPORARY_STRING_LENGTH];
    snprintf(tmp, TEMPORARY_STRING_LENGTH,
             "%s: %s snapshot at %s, for time-step goal %s\n",
             date_str.c_str(), executable_short_name.c_str(), grid.time->date().c_str(),
             grid.time->date(saving_after).c_str());

    if (!snapshots_file_is_ready) {

      // Prepare the snapshots file:
      ierr = nc.open(filename, PISM_WRITE); CHKERRQ(ierr);
      ierr = nc.def_time(config.get_string("time_dimension_name"),
                         config.get_string("calendar"),
                         grid.time->CF_units()); CHKERRQ(ierr);
      ierr = nc.close(); CHKERRQ(ierr);

      ierr = write_metadata(filename); CHKERRQ(ierr);

      snapshots_file_is_ready = true;
    }

    ierr = nc.open(filename, PISM_WRITE, true); CHKERRQ(ierr); // append==true
    ierr = nc.append_time(config.get_string("time_dimension_name"), grid.time->current()); CHKERRQ(ierr);
    ierr = nc.append_history(tmp); CHKERRQ(ierr); // append the history
    ierr = nc.close(); CHKERRQ(ierr);

    ierr = write_variables(filename, snapshot_vars, PISM_DOUBLE);

    grid.profiler->end(event_snapshots);

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

      ierr = set_output_size("-backup_size", "Sets the 'size' of a backup file.",
                             "small", backup_vars); CHKERRQ(ierr);
    }
    ierr = PetscOptionsEnd(); CHKERRQ(ierr);

    last_backup_time = 0.0;

    return 0;
  }

  //! Write a backup (i.e. an intermediate result of a run).
PetscErrorCode IceModel::write_backup() {
  PetscErrorCode ierr;
  double wall_clock_hours;
  PIO nc(grid.com, grid.rank, grid.config.get_string("output_format"));

  if (grid.rank == 0) {
    PetscLogDouble current_time;
    ierr = PetscGetTime(&current_time); CHKERRQ(ierr);
    wall_clock_hours = (current_time - start_time) / 3600.0;
  }

  MPI_Bcast(&wall_clock_hours, 1, MPI_DOUBLE, 0, grid.com);

  if (wall_clock_hours - last_backup_time < backup_interval)
    return 0;

  grid.profiler->begin(event_backups);

  last_backup_time = wall_clock_hours;

  // create a history string:
  string date_str = pism_timestamp();
  char tmp[TEMPORARY_STRING_LENGTH];
  snprintf(tmp, TEMPORARY_STRING_LENGTH,
           "%s automatic backup at %s, %3.3f hours after the beginning of the run\n",
           executable_short_name.c_str(), grid.time->date().c_str(), wall_clock_hours);

  ierr = verbPrintf(2, grid.com,
                    "  Saving an automatic backup to '%s' (%1.3f hours after the beginning of the run)\n",
                    backup_filename.c_str(), wall_clock_hours); CHKERRQ(ierr);

  stampHistory(tmp);

  // write metadata:
  ierr = nc.open(backup_filename, PISM_WRITE); CHKERRQ(ierr);
  ierr = nc.def_time(config.get_string("time_dimension_name"),
                     config.get_string("calendar"),
                     grid.time->CF_units()); CHKERRQ(ierr);
  ierr = nc.append_time(config.get_string("time_dimension_name"), grid.time->current()); CHKERRQ(ierr);
  ierr = nc.close(); CHKERRQ(ierr);

  // Write metadata *before* variables:
  ierr = write_metadata(backup_filename); CHKERRQ(ierr);

  ierr = write_variables(backup_filename, backup_vars, PISM_DOUBLE); CHKERRQ(ierr);

  // Also flush time-series:
  ierr = flush_timeseries(); CHKERRQ(ierr);

  grid.profiler->end(event_backups);

  return 0;
}

