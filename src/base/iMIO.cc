// Copyright (C) 2004-2016 Jed Brown, Ed Bueler and Constantine Khroulev
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

#include <algorithm>
#include <sstream>
#include <set>

#include "iceModel.hh"

#include "base/basalstrength/PISMYieldStress.hh"
#include "base/calving/CalvingAtThickness.hh"
#include "base/calving/EigenCalving.hh"
#include "base/calving/vonMisesCalving.hh"
#include "base/calving/FloatKill.hh"
#include "base/calving/OceanKill.hh"
#include "base/energy/BedThermalUnit.hh"
#include "base/hydrology/PISMHydrology.hh"
#include "base/stressbalance/PISMStressBalance.hh"
#include "base/util/IceGrid.hh"
#include "base/util/PISMConfigInterface.hh"
#include "base/util/PISMDiagnostic.hh"
#include "base/util/PISMTime.hh"
#include "base/util/error_handling.hh"
#include "base/util/io/PIO.hh"
#include "base/util/pism_options.hh"
#include "coupler/PISMOcean.hh"
#include "coupler/PISMSurface.hh"
#include "earth/PISMBedDef.hh"
#include "base/util/PISMVars.hh"
#include "base/util/io/io_helpers.hh"
#include "base/util/Profiling.hh"
#include "base/util/pism_utilities.hh"

namespace pism {

//! Save model state in NetCDF format.
/*!
Optionally allows saving of full velocity field.

Calls dumpToFile() to do the actual work.
 */
void  IceModel::writeFiles(const std::string &default_filename) {
  std::string config_out;

  stampHistoryEnd();

  options::String filename("-o", "Output file name", default_filename);

  if (!ends_with(filename, ".nc")) {
    m_log->message(2,
               "PISM WARNING: output file name does not have the '.nc' suffix!\n");
  }

  if (get_output_size("-o_size") != "none") {
    m_log->message(2,
               "Writing model state to file `%s'\n", filename->c_str());
    dumpToFile(filename);
  }
}

//! \brief Write metadata (global attributes, overrides and mapping parameters) to a file.
void IceModel::write_metadata(const PIO &nc, MetadataFlag flag) {

  if (flag & WRITE_MAPPING) {
    if (not nc.inq_var(m_mapping.get_name())) {
      nc.redef();
      nc.def_var(m_mapping.get_name(), PISM_DOUBLE,
                 std::vector<std::string>());
    }
    io::write_attributes(nc, m_mapping, PISM_DOUBLE, false);
  }

  if (flag & WRITE_RUN_STATS) {
    update_run_stats();
    if (not nc.inq_var(m_run_stats.get_name())) {
      nc.redef();
      nc.def_var(m_run_stats.get_name(), PISM_DOUBLE,
                 std::vector<std::string>());
    }
    io::write_attributes(nc, m_run_stats, PISM_DOUBLE, false);
  }

  if (flag & WRITE_GLOBAL_ATTRIBUTES) {
    io::write_global_attributes(nc, m_output_global_attributes);
  }

  // write configuration parameters to the file:
  m_config->write(nc);
}


void IceModel::dumpToFile(const std::string &filename) {
  const Profiling &profiling = m_ctx->profiling();

  PIO nc(m_grid->com, m_config->get_string("output.format"));

  profiling.begin("model state dump");

  // Prepare the file
  std::string time_name = m_config->get_string("time.dimension_name");
  nc.open(filename, PISM_READWRITE_MOVE);

  // Write metadata *before* everything else:
  write_metadata(nc, WRITE_ALL);

  io::define_time(nc, time_name, m_time->calendar(),
                  m_time->CF_units_string(),
                  m_sys);

  io::append_time(nc, time_name, m_time->current());

  write_model_state(nc);

  nc.close();

  profiling.end("model state dump");
}

//! \brief Writes variables listed in vars to filename, using nctype to write
//! fields stored in dedicated IceModelVecs.
void IceModel::write_variables(const PIO &nc, const std::set<std::string> &vars_input,
                               IO_Type nctype) {
  std::set<std::string> vars = vars_input;
  const IceModelVec *v;

  // Define all the variables:
  {
    std::set<std::string>::iterator i;
    for (i = vars.begin(); i != vars.end(); ++i) {

      if (m_grid->variables().is_available(*i)) {
        v = m_grid->variables().get(*i);
        // It has dedicated storage.
        if (*i == "mask") {
          v->define(nc, PISM_BYTE); // use the default data type
        } else {
          v->define(nc, nctype);
        }
      } else {
        // It might be a diagnostic quantity
        Diagnostic::Ptr diag = m_diagnostics[*i];

        if (diag) {
          diag->define(nc);
        }
      }
    }

    if (m_beddef != NULL) {
      m_beddef->define_variables(vars, nc, nctype);
    }

    if (m_btu != NULL) {
      m_btu->define_variables(vars, nc, nctype);
    }

    if (m_basal_yield_stress_model != NULL) {
      m_basal_yield_stress_model->define_variables(vars, nc, nctype);
    }

    if (m_stress_balance != NULL) {
      m_stress_balance->define_variables(vars, nc, nctype);
    } else {
      throw RuntimeError("PISM ERROR: stress_balance == NULL");
    }

    if (m_subglacial_hydrology != NULL) {
      m_subglacial_hydrology->define_variables(vars, nc, nctype);
    }

    if (m_surface != NULL) {
      m_surface->define_variables(vars, nc, nctype);
    } else {
      throw RuntimeError("PISM ERROR: surface == NULL");
    }

    if (m_ocean != NULL) {
      m_ocean->define_variables(vars, nc, nctype);
    } else {
      throw RuntimeError("PISM ERROR: ocean == NULL");
    }

    if (m_ocean_kill_calving != NULL) {
      m_ocean_kill_calving->define_variables(vars, nc, nctype);
    }

    if (m_float_kill_calving != NULL) {
      m_float_kill_calving->define_variables(vars, nc, nctype);
    }

    if (m_thickness_threshold_calving != NULL) {
      m_thickness_threshold_calving->define_variables(vars, nc, nctype);
    }

    if (m_eigen_calving != NULL) {
      m_eigen_calving->define_variables(vars, nc, nctype);
    }

    if (m_vonmises_calving != NULL) {
      m_vonmises_calving->define_variables(vars, nc, nctype);
    }

  }
  // Write all the IceModel variables:

  // Make a copy to avoid modifying the container we're iterating over.
  std::set<std::string> vars_copy = vars;
  std::set<std::string>::iterator i;
  for (i = vars.begin(); i != vars.end(); ++i) {
    if (m_grid->variables().is_available(*i)) {
      m_grid->variables().get(*i)->write(nc);

      // note that it only erases variables that were found (and
      // saved)
      vars_copy.erase(*i);
    }
  }
  vars = vars_copy;

  // Write bed-deformation-related variables:
  if (m_beddef != NULL) {
    m_beddef->write_variables(vars, nc);
  }

  // Write BedThermalUnit variables:
  if (m_btu != NULL) {
    m_btu->write_variables(vars, nc);
  }

  if (m_basal_yield_stress_model != NULL) {
    m_basal_yield_stress_model->write_variables(vars, nc);
  }

  // Write stress balance-related variables:
  if (m_stress_balance != NULL) {
    m_stress_balance->write_variables(vars, nc);
  } else {
    throw RuntimeError("PISM ERROR: stress_balance == NULL");
  }

  if (m_subglacial_hydrology != NULL) {
    m_subglacial_hydrology->write_variables(vars, nc);
  }

  // Ask boundary models to write their variables:
  if (m_surface != NULL) {
    m_surface->write_variables(vars, nc);
  } else {
    throw RuntimeError("PISM ERROR: surface == NULL");
  }
  if (m_ocean != NULL) {
    m_ocean->write_variables(vars, nc);
  } else {
    throw RuntimeError("PISM ERROR: ocean == NULL");
  }

  if (m_ocean_kill_calving != NULL) {
    m_ocean_kill_calving->write_variables(vars, nc);
  }

  if (m_float_kill_calving != NULL) {
    m_float_kill_calving->write_variables(vars, nc);
  }

  if (m_thickness_threshold_calving != NULL) {
    m_thickness_threshold_calving->write_variables(vars, nc);
  }

  if (m_eigen_calving != NULL) {
    m_eigen_calving->write_variables(vars, nc);
  }

  if (m_vonmises_calving != NULL) {
    m_vonmises_calving->write_variables(vars, nc);
  }

  // All the remaining names in vars must be of diagnostic quantities.
  for (i = vars.begin(); i != vars.end();) {
    Diagnostic::Ptr diag = m_diagnostics[*i];

    if (not diag) {
      ++i;
    } else {
      IceModelVec::Ptr v_diagnostic = diag->compute();

      v_diagnostic->write_in_glaciological_units = true;
      v_diagnostic->write(nc);

      vars.erase(i++);
    }
  }

  // FIXME: we need a way of figuring out if a sub-model did or did not write
  // something.

  if (not vars.empty()) {
    int threshold = 3;
    m_log->message(threshold,
               "PISM WARNING: the following variables were *not* written by PISM core (IceModel): ");
    for (i = vars.begin(); i != vars.end(); ++i) {
      m_log->message(threshold, "%s ", i->c_str());
    }
    m_log->message(threshold, "\n");
  }
}


void IceModel::write_model_state(const PIO &nc) {
  std::string o_size = get_output_size("-o_size");

#if (PISM_USE_PROJ4==1)
  if (m_output_global_attributes.has_attribute("proj4")) {
    m_output_vars.insert("lon_bnds");
    m_output_vars.insert("lat_bnds");
    m_latitude.metadata().set_string("bounds", "lat_bnds");
    m_longitude.metadata().set_string("bounds", "lon_bnds");
  }
#elif (PISM_USE_PROJ4==0)
  // do nothing
#else  // PISM_USE_PROJ4 is not set
#error "PISM build system error: PISM_USE_PROJ4 is not set."
#endif

  write_variables(nc, m_output_vars, PISM_DOUBLE);
}


//! Manage regridding based on user options.  Call IceModelVec::regrid() to do each selected variable.
/*!
  For each variable selected by option `-regrid_vars`, we regrid it onto the current grid from
  the NetCDF file specified by `-regrid_file`.

  The default, if `-regrid_vars` is not given, is to regrid the 3
  dimensional quantities `m_ice_age`, `Tb3` and either `m_ice_temperature` or `m_ice_enthalpy`. This is
  consistent with one standard purpose of regridding, which is to stick with
  current geometry through the downscaling procedure. Most of the time the user
  should carefully specify which variables to regrid.

  This `dimensions` argument can be 2 (regrid 2D variables only), 3 (3D
  only) and 0 (everything).
*/
void IceModel::regrid(int dimensions) {

  if (not (dimensions == 0 ||
           dimensions == 2 ||
           dimensions == 3)) {
    throw RuntimeError("dimensions can only be 0 (all), 2 or 3");
  }

  options::String regrid_filename("-regrid_file", "Specifies the file to regrid from");

  options::StringSet regrid_vars("-regrid_vars",
                                 "Specifies the list of variables to regrid",
                                 "");


  // Return if no regridding is requested:
  if (not regrid_filename.is_set()) {
     return;
  }

  if (dimensions != 0) {
    m_log->message(2, "regridding %dD variables from file %s ...\n",
               dimensions, regrid_filename->c_str());
  } else {
    m_log->message(2, "regridding from file %s ...\n",regrid_filename->c_str());
  }

  if (regrid_vars->empty()) {
    // defaults if user gives no regrid_vars list
    regrid_vars->insert("litho_temp");

    if (m_config->get_boolean("age.enabled")) {
      regrid_vars->insert("age");
    }

    if (m_config->get_boolean("energy.temperature_based")) {
      regrid_vars->insert("temp");
    } else {
      regrid_vars->insert("enthalpy");
    }
  }

  {
    PIO regrid_file(m_grid->com, "guess_mode");
    regrid_file.open(regrid_filename, PISM_READONLY);

    if (dimensions == 0) {
      regrid_variables(regrid_file, regrid_vars, 2);
      regrid_variables(regrid_file, regrid_vars, 3);
    } else {
      regrid_variables(regrid_file, regrid_vars, dimensions);
    }
  }
}

void IceModel::regrid_variables(const PIO &regrid_file, const std::set<std::string> &vars,
                                unsigned int ndims) {

  std::set<std::string>::iterator i;
  for (i = vars.begin(); i != vars.end(); ++i) {

    if (not m_grid->variables().is_available(*i)) {
      continue;
    }

    // FIXME: remove const_cast. This is bad.
    IceModelVec *v = const_cast<IceModelVec*>(m_grid->variables().get(*i));
    SpatialVariableMetadata &m = v->metadata();

    if (v->get_ndims() != ndims) {
      continue;
    }

    std::string pism_intent = m.get_string("pism_intent");
    if (pism_intent != "model_state") {
      m_log->message(2, "  WARNING: skipping '%s' (only model_state variables can be regridded)...\n",
                 i->c_str());
      continue;
    }

    if (*i == "enthalpy") {
      init_enthalpy(regrid_file,
                    true,       // regrid
                    0);         // record index (ignored)
      continue;
    }

    v->regrid(regrid_file, CRITICAL);
  }

  // Check the range of the ice thickness.
  {
    double max_thickness = m_ice_thickness.range().max,
      Lz = m_grid->Lz();

    if (max_thickness >= Lz + 1e-6) {
      throw RuntimeError::formatted("Maximum ice thickness (%f meters)\n"
                                    "exceeds the height of the computational domain (%f meters).",
                                    max_thickness, Lz);
    }
  }
}

/**
 * Initialize enthalpy from a file that does not contain it using "temp" and "liqfrac".
 *
 * @param input_file input file
 *
 * @param do_regrid use regridding if 'true', otherwise assume that the
 *                  input file has the same grid
 * @param last_record the record to use when 'do_regrid==false'.
 *
 */
void IceModel::init_enthalpy(const PIO &input_file,
                             bool do_regrid, int last_record) {

  const bool
    enthalpy_exists = input_file.inq_var("enthalpy"),
    temp_exists     = input_file.inq_var("temp"),
    liqfrac_exists  = input_file.inq_var("liqfrac");

  if (enthalpy_exists) {
    if (do_regrid) {
      m_ice_enthalpy.regrid(input_file, CRITICAL);
    } else {
      m_ice_enthalpy.read(input_file, last_record);
    }
  } else if (temp_exists) {
    IceModelVec3
      &temp    = m_work3d,
      &liqfrac = m_ice_enthalpy;

    SpatialVariableMetadata enthalpy_metadata = m_ice_enthalpy.metadata();
    temp.set_name("temp");
    temp.metadata(0).set_name("temp");
    temp.set_attrs("temporary", "ice temperature", "Kelvin",
                   "land_ice_temperature");

    if (do_regrid) {
      temp.regrid(input_file, CRITICAL);
    } else {
      temp.read(input_file, last_record);
    }

    if (liqfrac_exists and not m_config->get_boolean("energy.temperature_based")) {
      liqfrac.set_name("liqfrac");
      liqfrac.metadata(0).set_name("liqfrac");
      liqfrac.set_attrs("temporary", "ice liquid water fraction",
                        "1", "");

      if (do_regrid) {
        liqfrac.regrid(input_file, CRITICAL);
      } else {
        liqfrac.read(input_file, last_record);
      }

      m_log->message(2,
                 "* Computing enthalpy using ice temperature,"
                 "  liquid water fraction and thickness...\n");
      compute_enthalpy(temp, liqfrac, m_ice_enthalpy);
    } else {
      m_log->message(2,
                 "* Computing enthalpy using ice temperature and thickness...\n");
      compute_enthalpy_cold(temp, m_ice_enthalpy);
    }

    m_ice_enthalpy.metadata() = enthalpy_metadata;
  } else {
    throw RuntimeError::formatted("neither enthalpy nor temperature was found in '%s'.\n",
                                  input_file.inq_filename().c_str());
  }
}

//! Initializes the snapshot-saving mechanism.
void IceModel::init_snapshots() {
  m_current_snapshot = 0;

  options::String save_file("-save_file", "Specifies a snapshot filename",
                            m_snapshots_filename);
  m_snapshots_filename = save_file;

  options::String save_times("-save_times",
                             "Gives a list or a MATLAB-style range of times to save snapshots at");

  bool split = options::Bool("-save_split", "Specifies whether to save snapshots to separate files");

  m_snapshot_vars = output_size_from_option("-save_size", "Sets the 'size' of a snapshot file.",
                                          "small");


  if (save_file.is_set() ^ save_times.is_set()) {
    throw RuntimeError("you need to specify both -save_file and -save_times to save snapshots.");
  }

  if (!save_file.is_set() && !save_times.is_set()) {
    m_save_snapshots = false;
    return;
  }

  try {
    m_time->parse_times(save_times, m_snapshot_times);
  } catch (RuntimeError &e) {
    e.add_context("parsing the -save_times argument");
    throw;
  }

  if (m_snapshot_times.size() == 0) {
    throw RuntimeError("no argument for -save_times option.");
  }

  m_save_snapshots = true;
  m_snapshots_file_is_ready = false;
  m_split_snapshots = false;

  if (split) {
    m_split_snapshots = true;
  } else if (!ends_with(m_snapshots_filename, ".nc")) {
    m_log->message(2,
               "PISM WARNING: snapshots file name does not have the '.nc' suffix!\n");
  }

  if (split) {
    m_log->message(2, "saving snapshots to '%s+year.nc'; ",
               m_snapshots_filename.c_str());
  } else {
    m_log->message(2, "saving snapshots to '%s'; ",
               m_snapshots_filename.c_str());
  }

  m_log->message(2, "times requested: %s\n", save_times->c_str());
}

  //! Writes a snapshot of the model state (if necessary)
void IceModel::write_snapshot() {
  double saving_after = -1.0e30; // initialize to avoid compiler warning; this
  // value is never used, because saving_after
  // is only used if save_now == true, and in
  // this case saving_after is guaranteed to be
  // initialized. See the code below.
  char filename[PETSC_MAX_PATH_LEN];

  // determine if the user set the -save_times and -save_file options
  if (!m_save_snapshots) {
    return;
  }

  // do we need to save *now*?
  if ((m_time->current() >= m_snapshot_times[m_current_snapshot]) and
      (m_current_snapshot < m_snapshot_times.size())) {
    saving_after = m_snapshot_times[m_current_snapshot];

    while ((m_current_snapshot < m_snapshot_times.size()) and
           (m_snapshot_times[m_current_snapshot] <= m_time->current())) {
      m_current_snapshot++;
    }
  } else {
    // we don't need to save now, so just return
    return;
  }

  // flush time-series buffers
  flush_timeseries();

  if (m_split_snapshots) {
    m_snapshots_file_is_ready = false;    // each snapshot is written to a separate file
    snprintf(filename, PETSC_MAX_PATH_LEN, "%s_%s.nc",
             m_snapshots_filename.c_str(), m_time->date(saving_after).c_str());
  } else {
    strncpy(filename, m_snapshots_filename.c_str(), PETSC_MAX_PATH_LEN);
  }

  m_log->message(2,
             "\nsaving snapshot to %s at %s, for time-step goal %s\n\n",
             filename, m_time->date().c_str(),
             m_time->date(saving_after).c_str());

  PIO nc(m_grid->com, m_config->get_string("output.format"));

  if (not m_snapshots_file_is_ready) {
    // Prepare the snapshots file:
    nc.open(filename, PISM_READWRITE_MOVE);
    io::define_time(nc, m_config->get_string("time.dimension_name"),
                m_time->calendar(),
                m_time->CF_units_string(),
                m_sys);

    m_snapshots_file_is_ready = true;
  } else {
    // In this case the snapshot file is should be present.
    nc.open(filename, PISM_READWRITE);
  }

  // write metadata to the file *every time* we update it
  write_metadata(nc, WRITE_ALL);

  io::append_time(nc, m_config->get_string("time.dimension_name"), m_time->current());

  write_variables(nc, m_snapshot_vars, PISM_DOUBLE);

  {
    // find out how much time passed since the beginning of the run
    double wall_clock_hours = pism::wall_clock_hours(m_grid->com, m_start_time);

    // Get time length now, i.e. after writing variables. This forces PISM to call PIO::enddef(), so
    // that the length of the time dimension is up to date.
    unsigned int time_length = nc.inq_dimlen(m_config->get_string("time.dimension_name"));

    // make sure that time_start is valid even if time_length is zero
    size_t time_start = 0;
    if (time_length > 0) {
      time_start = static_cast<size_t>(time_length - 1);
    } else {
      time_start = 0;
    }

    io::write_timeseries(nc, m_timestamp, time_start, wall_clock_hours);
  }

  nc.close();
}

//! Initialize the backup (snapshot-on-wallclock-time) mechanism.
void IceModel::init_backups() {

  m_backup_interval = m_config->get_double("output.backup_interval");

  options::String backup_file("-o", "Output file name");
  if (backup_file.is_set()) {
    m_backup_filename = pism_filename_add_suffix(backup_file, "_backup", "");
  } else {
    m_backup_filename = "pism_backup.nc";
  }

  m_backup_interval = options::Real("-backup_interval",
                                  "Automatic backup interval, hours", m_backup_interval);

  m_backup_vars = output_size_from_option("-backup_size", "Sets the 'size' of a backup file.",
                                        "small");

  m_last_backup_time = 0.0;
}

  //! Write a backup (i.e. an intermediate result of a run).
void IceModel::write_backup() {
  double wall_clock_hours = pism::wall_clock_hours(m_grid->com, m_start_time);

  if (wall_clock_hours - m_last_backup_time < m_backup_interval) {
    return;
  }

  const Profiling &profiling = m_ctx->profiling();
  profiling.begin("backup");

  m_last_backup_time = wall_clock_hours;

  // create a history string:
  std::string date_str = pism_timestamp(m_grid->com);
  char tmp[TEMPORARY_STRING_LENGTH];
  snprintf(tmp, TEMPORARY_STRING_LENGTH,
           "PISM automatic backup at %s, %3.3f hours after the beginning of the run\n",
           m_time->date().c_str(), wall_clock_hours);

  PetscLogDouble backup_start_time = GetTime();

  m_log->message(2,
                 "  [%s] Saving an automatic backup to '%s' (%1.3f hours after the beginning of the run)\n",
                 pism_timestamp(m_grid->com).c_str(), m_backup_filename.c_str(), wall_clock_hours);

  stampHistory(tmp);

  PIO nc(m_grid->com, m_config->get_string("output.format"));

  // write metadata:
  nc.open(m_backup_filename, PISM_READWRITE_MOVE);
  io::define_time(nc, m_config->get_string("time.dimension_name"),
              m_time->calendar(),
              m_time->CF_units_string(),
              m_sys);
  io::append_time(nc, m_config->get_string("time.dimension_name"), m_time->current());

  // Write metadata *before* variables:
  write_metadata(nc, WRITE_ALL);

  write_variables(nc, m_backup_vars, PISM_DOUBLE);

  nc.close();

  // Also flush time-series:
  flush_timeseries();

  PetscLogDouble backup_end_time = GetTime();
  m_log->message(2,
                 "  [%s] Done saving an automatic backup in %f seconds (%f minutes).\n",
                 pism_timestamp(m_grid->com).c_str(),
                 backup_end_time - backup_start_time,
                 (backup_end_time - backup_start_time) / 60.0);

  profiling.end("backup");
}


} // end of namespace pism
