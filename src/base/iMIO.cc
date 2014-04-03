// Copyright (C) 2004-2014 Jed Brown, Ed Bueler and Constantine Khroulev
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
#include "PISMHydrology.hh"
#include "PISMStressBalance.hh"
#include "PISMSurface.hh"
#include "PISMOcean.hh"
#include "pism_options.hh"
#include "IceGrid.hh"
#include "PISMTime.hh"
#include "PISMDiagnostic.hh"
#include "PISMOceanKill.hh"
#include "PISMFloatKill.hh"
#include "PISMCalvingAtThickness.hh"
#include "PISMEigenCalving.hh"

//! Save model state in NetCDF format.
/*!
Optionally allows saving of full velocity field.

Calls dumpToFile() to do the actual work.
 */
PetscErrorCode  IceModel::writeFiles(std::string default_filename) {
  PetscErrorCode ierr;
  std::string filename = default_filename,
    config_out;
  bool o_set, dump_config;

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
    // Here "false" means "do not append"; creates a new file and moves
    // the old one aside
    ierr = config.write(config_out, false); CHKERRQ(ierr);
  }

  return 0;
}

//! \brief Write metadata (global attributes, overrides and mapping parameters) to a file.
PetscErrorCode IceModel::write_metadata(const PIO &nc, bool write_mapping,
                                        bool write_run_stats) {
  PetscErrorCode ierr;

  if (write_mapping) {
    bool mapping_exists = false;
    ierr = nc.inq_var(mapping.get_name(), mapping_exists); CHKERRQ(ierr);
    if (mapping_exists == false) {
      ierr = nc.redef(); CHKERRQ(ierr);
      ierr = nc.def_var(mapping.get_name(), PISM_DOUBLE, std::vector<std::string>()); CHKERRQ(ierr);
    }
    ierr = nc.write_attributes(mapping, PISM_DOUBLE, false); CHKERRQ(ierr);
  }

  if (write_run_stats) {
    ierr = update_run_stats(); CHKERRQ(ierr);
    bool run_stats_exists = false;
    ierr = nc.inq_var(run_stats.get_name(), run_stats_exists); CHKERRQ(ierr);
    if (run_stats_exists == false) {
      ierr = nc.redef(); CHKERRQ(ierr);
      ierr = nc.def_var(run_stats.get_name(), PISM_DOUBLE, std::vector<std::string>()); CHKERRQ(ierr);
    }
    ierr = nc.write_attributes(run_stats, PISM_DOUBLE, false); CHKERRQ(ierr);
  }

  ierr = nc.write_global_attributes(global_attributes); CHKERRQ(ierr);

  bool override_used;
  ierr = PISMOptionsIsSet("-config_override", override_used); CHKERRQ(ierr);
  if (override_used) {
    overrides.update_from(config);
    ierr = overrides.write(nc); CHKERRQ(ierr);
  }

  return 0;
}


PetscErrorCode IceModel::dumpToFile(std::string filename) {
  PetscErrorCode ierr;
  PIO nc(grid, config.get_string("output_format"));

  // Prepare the file
  std::string time_name = config.get_string("time_dimension_name");
  ierr = nc.open(filename, PISM_READWRITE_MOVE); CHKERRQ(ierr);
  ierr = nc.def_time(time_name, grid.time->calendar(),
                     grid.time->CF_units_string()); CHKERRQ(ierr);
  ierr = nc.append_time(time_name, grid.time->current()); CHKERRQ(ierr);

  // Write metadata *before* variables:
  ierr = write_metadata(nc, true, true); CHKERRQ(ierr);

  ierr = write_model_state(nc);  CHKERRQ(ierr);

  ierr = nc.close(); CHKERRQ(ierr);

  return 0;
}

//! \brief Writes variables listed in vars to filename, using nctype to write
//! fields stored in dedicated IceModelVecs.
PetscErrorCode IceModel::write_variables(const PIO &nc, std::set<std::string> vars,
                                         PISM_IO_Type nctype) {
  PetscErrorCode ierr;
  IceModelVec *v;

  // Define all the variables:
  {
    std::set<std::string>::iterator i;
    for (i = vars.begin(); i != vars.end(); ++i) {
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
    }

    if (beddef != NULL) {
      ierr = beddef->define_variables(vars, nc, nctype); CHKERRQ(ierr);
    }

    if (btu != NULL) {
      ierr = btu->define_variables(vars, nc, nctype); CHKERRQ(ierr);
    }

    if (basal_yield_stress_model != NULL) {
      ierr = basal_yield_stress_model->define_variables(vars, nc, nctype); CHKERRQ(ierr);
    }

    if (stress_balance != NULL) {
      ierr = stress_balance->define_variables(vars, nc, nctype); CHKERRQ(ierr);
    } else {
      SETERRQ(grid.com, 1,"PISM ERROR: stress_balance == NULL");
    }

    if (subglacial_hydrology != NULL) {
      ierr = subglacial_hydrology->define_variables(vars, nc, nctype); CHKERRQ(ierr);
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

    if (ocean_kill_calving != NULL) {
      ierr = ocean_kill_calving->define_variables(vars, nc, nctype); CHKERRQ(ierr);
    }

    if (float_kill_calving != NULL) {
      ierr = float_kill_calving->define_variables(vars, nc, nctype); CHKERRQ(ierr);
    }

    if (thickness_threshold_calving != NULL) {
      ierr = thickness_threshold_calving->define_variables(vars, nc, nctype); CHKERRQ(ierr);
    }

    if (eigen_calving != NULL) {
      ierr = eigen_calving->define_variables(vars, nc, nctype); CHKERRQ(ierr);
    }

  }
  // Write all the IceModel variables:
  std::set<std::string>::iterator i;
  for (i = vars.begin(); i != vars.end();) {
    v = variables.get(*i);

    if (v == NULL) {
      ++i;
    } else {
      ierr = v->write(nc); CHKERRQ(ierr); // use the default data type

      vars.erase(i++);          // note that it only erases variables that were
                                // found (and saved)
    }
  }

  // Write bed-deformation-related variables:
  if (beddef != NULL) {
    ierr = beddef->write_variables(vars, nc); CHKERRQ(ierr);
  }

  // Write PISMBedThermalUnit variables:
  if (btu != NULL) {
    ierr = btu->write_variables(vars, nc); CHKERRQ(ierr);
  }

  if (basal_yield_stress_model != NULL) {
    ierr = basal_yield_stress_model->write_variables(vars, nc); CHKERRQ(ierr);
  }

  // Write stress balance-related variables:
  if (stress_balance != NULL) {
    ierr = stress_balance->write_variables(vars, nc); CHKERRQ(ierr);
  } else {
    SETERRQ(grid.com, 1,"PISM ERROR: stress_balance == NULL");
  }

  if (subglacial_hydrology != NULL) {
    ierr = subglacial_hydrology->write_variables(vars, nc); CHKERRQ(ierr);
  }

  // Ask boundary models to write their variables:
  if (surface != NULL) {
    ierr = surface->write_variables(vars, nc); CHKERRQ(ierr);
  } else {
    SETERRQ(grid.com, 1,"PISM ERROR: surface == NULL");
  }
  if (ocean != NULL) {
    ierr = ocean->write_variables(vars, nc); CHKERRQ(ierr);
  } else {
    SETERRQ(grid.com, 1,"PISM ERROR: ocean == NULL");
  }

  if (ocean_kill_calving != NULL) {
    ierr = ocean_kill_calving->write_variables(vars, nc); CHKERRQ(ierr);
  }

  if (float_kill_calving != NULL) {
    ierr = float_kill_calving->write_variables(vars, nc); CHKERRQ(ierr);
  }

  if (thickness_threshold_calving != NULL) {
    ierr = thickness_threshold_calving->write_variables(vars, nc); CHKERRQ(ierr);
  }

  if (eigen_calving != NULL) {
    ierr = eigen_calving->write_variables(vars, nc); CHKERRQ(ierr);
  }

  // All the remaining names in vars must be of diagnostic quantities.
  for (i = vars.begin(); i != vars.end();) {
    PISMDiagnostic *diag = diagnostics[*i];

    if (diag == NULL)
      ++i;
    else {
      v = NULL;

      ierr = diag->compute(v); CHKERRQ(ierr);

      v->write_in_glaciological_units = true;
      ierr = v->write(nc, PISM_FLOAT); CHKERRQ(ierr); // diagnostic quantities are always written in float

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
      ierr = verbPrintf(threshold, grid.com, "%s ", i->c_str()); CHKERRQ(ierr);
    }
    ierr = verbPrintf(threshold, grid.com, "\n"); CHKERRQ(ierr);
  }

  return 0;
}


PetscErrorCode IceModel::write_model_state(const PIO &nc) {
  PetscErrorCode ierr;
  std::string o_size = get_output_size("-o_size");

#if (PISM_USE_PROJ4==1)
  if (global_attributes.has_attribute("proj4")) {
    output_vars.insert("lon_bnds");
    output_vars.insert("lat_bnds");
    vLatitude.metadata().set_string("bounds", "lat_bnds");
    vLongitude.metadata().set_string("bounds", "lon_bnds");
  }
#elif (PISM_USE_PROJ4==0)
  // do nothing
#else  // PISM_USE_PROJ4 is not set
#error "PISM build system error: PISM_USE_PROJ4 is not set."
#endif

  ierr = write_variables(nc, output_vars, PISM_DOUBLE); CHKERRQ(ierr);

  return 0;
}



  //! Read a saved PISM model state in NetCDF format, for complete initialization of an evolution or diagnostic run.
  /*!
    Before this is run, the method IceModel::grid_setup() determines the number of
    grid points (Mx,My,Mz,Mbz) and the dimensions (Lx,Ly,Lz) of the computational
    box from the same input file.
  */
PetscErrorCode IceModel::initFromFile(std::string filename) {
  PetscErrorCode  ierr;
  PIO nc(grid, "guess_mode");

  ierr = verbPrintf(2, grid.com, "initializing from NetCDF file '%s'...\n",
                    filename.c_str()); CHKERRQ(ierr);

  ierr = nc.open(filename, PISM_READONLY); CHKERRQ(ierr);

  // Find the index of the last record in the file:
  unsigned int last_record;
  ierr = nc.inq_nrecords(last_record); CHKERRQ(ierr); 
  last_record -= 1;

  // Read the model state, mapping and climate_steady variables:
  std::set<std::string> vars = variables.keys();

  std::set<std::string>::iterator i = vars.begin();
  while (i != vars.end()) {
    IceModelVec *var = variables.get(*i++);
    NCSpatialVariable &m = var->metadata();

    std::string intent = m.get_string("pism_intent");
    if ((intent == "model_state") || (intent == "mapping") ||
        (intent == "climate_steady")) {

      // skip "age", "enthalpy", and "Href" for now: we'll take care
      // of them a little later
      if (m.get_string("short_name") == "enthalpy" ||
          m.get_string("short_name") == "age"      ||
          m.get_string("short_name") == "Href")
        continue;

      ierr = var->read(filename, last_record); CHKERRQ(ierr);
    }
  }

  if (config.get_flag("do_cold_ice_methods")) {
    ierr = verbPrintf(3,grid.com,"  setting enthalpy from temperature...\n");
    CHKERRQ(ierr);
    ierr = compute_enthalpy_cold(T3, Enth3); CHKERRQ(ierr);
  }

  // check if the input file has Href; set to 0 if it is not present
  if (config.get_flag("part_grid")) {
    bool href_exists;
    ierr = nc.inq_var("Href", href_exists); CHKERRQ(ierr);

    if (href_exists == true) {
      ierr = vHref.read(filename, last_record); CHKERRQ(ierr);
    } else {
      ierr = verbPrintf(2, grid.com,
        "PISM WARNING: Href for PISM-PIK -part_grid not found in '%s'. Setting it to zero...\n",
                        filename.c_str()); CHKERRQ(ierr);
      ierr = vHref.set(0.0); CHKERRQ(ierr);
    }
  }

  // read the age field if present, otherwise set to zero
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


  // Initialize the enthalpy field by reading from a file or by using
  // temperature and liquid water fraction, or by using temperature
  // and assuming that the ice is cold.
  ierr = init_enthalpy(filename, false, last_record); CHKERRQ(ierr);

  std::string history;
  ierr = nc.get_att_text("PISM_GLOBAL", "history", history); CHKERRQ(ierr);
  global_attributes.set_string("history",
                               history + global_attributes.get_string("history"));

  ierr = nc.close(); CHKERRQ(ierr);

  return 0;
}

//! Manage regridding based on user options.  Call IceModelVec::regrid() to do each selected variable.
/*!
  For each variable selected by option `-regrid_vars`, we regrid it onto the current grid from
  the NetCDF file specified by `-regrid_file`.

  The default, if `-regrid_vars` is not given, is to regrid the 3
  dimensional quantities `tau3`, `Tb3` and either `T3` or `Enth3`. This is
  consistent with one standard purpose of regridding, which is to stick with
  current geometry through the downscaling procedure. Most of the time the user
  should carefully specify which variables to regrid.

  This `dimensions` argument can be 2 (regrid 2D variables only), 3 (3D
  only) and 0 (everything).
*/
PetscErrorCode IceModel::regrid(int dimensions) {
  PetscErrorCode ierr;
  std::string filename;
  bool regrid_vars_set, regrid_file_set;
  std::set<std::string> regrid_vars;

  if (! (dimensions == 0 ||
         dimensions == 2 ||
         dimensions == 3))
    SETERRQ(grid.com, 1, "dimensions can only be 0, 2 or 3");

  ierr = PetscOptionsBegin(grid.com, PETSC_NULL, "Options controlling regridding",
                           PETSC_NULL); CHKERRQ(ierr);
  {
    ierr = PISMOptionsString("-regrid_file", "Specifies the file to regrid from",
                             filename, regrid_file_set); CHKERRQ(ierr);

    ierr = PISMOptionsStringSet("-regrid_vars", "Specifies the list of variables to regrid",
                                "", regrid_vars, regrid_vars_set); CHKERRQ(ierr);
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

  if (regrid_vars.empty()) {
    // defaults if user gives no regrid_vars list
    regrid_vars.insert("litho_temp");

    if (config.get_flag("do_age"))
      regrid_vars.insert("age");

    if (config.get_flag("do_cold_ice_methods"))
      regrid_vars.insert("temp");
    else
      regrid_vars.insert("enthalpy");
  }

  if (dimensions == 0) {
    ierr = regrid_variables(filename, regrid_vars, 2); CHKERRQ(ierr);
    ierr = regrid_variables(filename, regrid_vars, 3); CHKERRQ(ierr);
  } else {
    ierr = regrid_variables(filename, regrid_vars, dimensions); CHKERRQ(ierr);
  }

  return 0;
}

PetscErrorCode IceModel::regrid_variables(std::string filename, std::set<std::string> vars, unsigned int ndims) {
  PetscErrorCode ierr;

  std::set<std::string>::iterator i;
  for (i = vars.begin(); i != vars.end(); ++i) {
    IceModelVec *v = variables.get(*i);

    if (v == NULL) continue;

    NCSpatialVariable &m = v->metadata();

    if (v->get_ndims() != ndims) continue;

    std::string pism_intent = m.get_string("pism_intent");
    if (pism_intent != "model_state") {
      ierr = verbPrintf(2, grid.com, "  WARNING: skipping '%s' (only model_state variables can be regridded)...\n",
                        i->c_str()); CHKERRQ(ierr);
      continue;
    }

    if (*i == "enthalpy") {
      ierr = init_enthalpy(filename, true, 0); CHKERRQ(ierr);
      continue;
    }

    ierr = v->regrid(filename, CRITICAL); CHKERRQ(ierr);

    // Check if the current variable is the same as
    // IceModel::ice_thickess, then check the range of the ice
    // thickness
    if (v == &this->ice_thickness) {
      double thk_min = 0.0, thk_max = 0.0;
      ierr = ice_thickness.range(thk_min, thk_max); CHKERRQ(ierr);

      if (thk_max >= grid.Lz + 1e-6) {
        PetscPrintf(grid.com,
                    "PISM ERROR: Maximum ice thickness (%f meters)\n"
                    "            exceeds the height of the computational domain (%f meters).\n"
                    "            Stopping...\n", thk_max, grid.Lz);
        PISMEnd();
      }
    }

  }

  return 0;
}

/**
 * Initialize enthalpy from a file that does not contain it using "temp" and "liqfrac".
 *
 * @param filename input file name
 *
 * @param do_regrid use regridding if 'true', otherwise assume that the
 *               input file has the same grid
 * @param last_record the record to use when 'do_regrid==false'.
 *
 * @return 0 on success
 */
PetscErrorCode IceModel::init_enthalpy(std::string filename,
                                       bool do_regrid, int last_record) {
  PetscErrorCode ierr;
  bool temp_exists  = false,
    liqfrac_exists  = false,
    enthalpy_exists = false;

  PIO nc(grid, "guess_mode");
  ierr = nc.open(filename, PISM_READONLY); CHKERRQ(ierr);
  ierr = nc.inq_var("enthalpy", enthalpy_exists); CHKERRQ(ierr);
  ierr = nc.inq_var("temp", temp_exists); CHKERRQ(ierr);
  ierr = nc.inq_var("liqfrac", liqfrac_exists); CHKERRQ(ierr);
  ierr = nc.close(); CHKERRQ(ierr);

  if (enthalpy_exists == true) {
    if (do_regrid) {
      ierr = Enth3.regrid(filename, CRITICAL); CHKERRQ(ierr);
    } else {
      ierr = Enth3.read(filename, last_record); CHKERRQ(ierr);
    }
  } else if (temp_exists == true) {
    IceModelVec3 &temp = vWork3d,
      &liqfrac         = Enth3;

    NCSpatialVariable enthalpy_metadata = Enth3.metadata();
    ierr = temp.set_name("temp"); CHKERRQ(ierr);
    ierr = temp.set_attrs("temporary", "ice temperature", "Kelvin",
                          "land_ice_temperature"); CHKERRQ(ierr);

    if (do_regrid) {
      ierr = temp.regrid(filename, CRITICAL); CHKERRQ(ierr);
    } else {
      ierr = temp.read(filename, last_record); CHKERRQ(ierr);
    }

    if (liqfrac_exists == true) {
      ierr = liqfrac.set_name("liqfrac"); CHKERRQ(ierr);
      ierr = liqfrac.set_attrs("temporary", "ice liquid water fraction",
                               "1", ""); CHKERRQ(ierr);

      if (do_regrid) {
        ierr = liqfrac.regrid(filename, CRITICAL); CHKERRQ(ierr);
      } else {
        ierr = liqfrac.read(filename, last_record); CHKERRQ(ierr);
      }

      ierr = verbPrintf(2, grid.com,
                        "* Computing enthalpy using ice temperature,"
                        "  liquid water fraction and thickness...\n"); CHKERRQ(ierr);
      ierr = compute_enthalpy(temp, liqfrac, Enth3); CHKERRQ(ierr);
    } else {
      ierr = verbPrintf(2, grid.com,
                        "* Computing enthalpy using ice temperature and thickness...\n");
      CHKERRQ(ierr);
      ierr = compute_enthalpy_cold(temp, Enth3); CHKERRQ(ierr);
    }

    Enth3.metadata() = enthalpy_metadata;
  } else {
    PetscPrintf(grid.com, "PISM ERROR: neither %s nor %s was found in '%s'.\n",
                "enthalpy", "temperature", filename.c_str());
    PISMEnd();
  }

  return 0;
}

//! Initializes the snapshot-saving mechanism.
PetscErrorCode IceModel::init_snapshots() {
  PetscErrorCode ierr;
  bool save_times_set, save_file_set, split;
  std::string tmp;
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

  ierr = grid.time->parse_times(tmp, snapshot_times);
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

  // flush time-series buffers
  ierr = flush_timeseries(); CHKERRQ(ierr);

  if (split_snapshots) {
    snapshots_file_is_ready = false;    // each snapshot is written to a separate file
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

  PIO nc(grid, grid.config.get_string("output_format"));

  if (snapshots_file_is_ready == false) {
    // Prepare the snapshots file:
    ierr = nc.open(filename, PISM_READWRITE_MOVE); CHKERRQ(ierr);
    ierr = nc.def_time(config.get_string("time_dimension_name"),
                       grid.time->calendar(),
                       grid.time->CF_units_string()); CHKERRQ(ierr);

    ierr = write_metadata(nc, true, true); CHKERRQ(ierr);

    snapshots_file_is_ready = true;
  } else {
    // In this case the snapshot file is should be present.
    ierr = nc.open(filename, PISM_READWRITE); CHKERRQ(ierr);
  }

  unsigned int time_length = 0;

  ierr = nc.append_time(config.get_string("time_dimension_name"), grid.time->current()); CHKERRQ(ierr);
  ierr = nc.inq_dimlen(config.get_string("time_dimension_name"), time_length); CHKERRQ(ierr);

  ierr = write_variables(nc, snapshot_vars, PISM_DOUBLE);

  {
    // find out how much time passed since the beginning of the run
    double wall_clock_hours;
    if (grid.rank == 0) {
      PetscLogDouble current_time;
      ierr = PISMGetTime(&current_time); CHKERRQ(ierr);
      wall_clock_hours = (current_time - start_time) / 3600.0;
    }

    MPI_Bcast(&wall_clock_hours, 1, MPI_DOUBLE, 0, grid.com);

    ierr = nc.write_timeseries(timestamp, static_cast<size_t>(time_length - 1),
                               wall_clock_hours); CHKERRQ(ierr);
  }

  ierr = nc.close(); CHKERRQ(ierr);

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

  if (grid.rank == 0) {
    PetscLogDouble current_time;
    ierr = PISMGetTime(&current_time); CHKERRQ(ierr);
    wall_clock_hours = (current_time - start_time) / 3600.0;
  }

  MPI_Bcast(&wall_clock_hours, 1, MPI_DOUBLE, 0, grid.com);

  if (wall_clock_hours - last_backup_time < backup_interval)
    return 0;

  last_backup_time = wall_clock_hours;

  // create a history string:
  std::string date_str = pism_timestamp();
  char tmp[TEMPORARY_STRING_LENGTH];
  snprintf(tmp, TEMPORARY_STRING_LENGTH,
           "%s automatic backup at %s, %3.3f hours after the beginning of the run\n",
           executable_short_name.c_str(), grid.time->date().c_str(), wall_clock_hours);

  ierr = verbPrintf(2, grid.com,
                    "  Saving an automatic backup to '%s' (%1.3f hours after the beginning of the run)\n",
                    backup_filename.c_str(), wall_clock_hours); CHKERRQ(ierr);

  stampHistory(tmp);

  PIO nc(grid, grid.config.get_string("output_format"));

  // write metadata:
  ierr = nc.open(backup_filename, PISM_READWRITE_MOVE); CHKERRQ(ierr);
  ierr = nc.def_time(config.get_string("time_dimension_name"),
                     grid.time->calendar(),
                     grid.time->CF_units_string()); CHKERRQ(ierr);
  ierr = nc.append_time(config.get_string("time_dimension_name"), grid.time->current()); CHKERRQ(ierr);

  // Write metadata *before* variables:
  ierr = write_metadata(nc, true, true); CHKERRQ(ierr);

  ierr = write_variables(nc, backup_vars, PISM_DOUBLE); CHKERRQ(ierr);

  ierr = nc.close(); CHKERRQ(ierr);

  // Also flush time-series:
  ierr = flush_timeseries(); CHKERRQ(ierr);

  return 0;
}

