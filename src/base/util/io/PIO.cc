// Copyright (C) 2012, 2013, 2014 PISM Authors
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

#include <petscvec.h>

#include "PIO.hh"
#include "IceGrid.hh"
#include "pism_const.hh"
#include "LocalInterpCtx.hh"
#include "NCVariable.hh"
#include "PISMConfig.hh"
#include "PISMTime.hh"
#include "PISMNC3File.hh"
#include "PISMNC4_Quilt.hh"
#include <cassert>

#if (PISM_USE_PARALLEL_NETCDF4==1)
#include "PISMNC4_Par.hh"
#endif

#if (PISM_USE_PNETCDF==1)
#include "PISMPNCFile.hh"
#endif

#if (PISM_USE_HDF5==1)
#include "PISMNC4_HDF5.hh"
#endif

#ifdef PISM_USE_TR1
#include <tr1/memory>
using std::tr1::shared_ptr;
#else
#include <memory>
using std::shared_ptr;
#endif

#include "error_handling.hh"

namespace pism {

using std::string;
using std::vector;

static NCFile::Ptr create_backend(MPI_Comm com, string mode) {
  if (mode == "netcdf3") {
    return NCFile::Ptr(new NC3File(com));
  } else if (mode.find("quilt") == 0) {
    size_t n = mode.find(":");
    int compression_level = 0;

    if (n != string::npos) {
      mode.replace(0, 6, "");   // 6 is the length of "quilt:"
      char *endptr;
      compression_level = strtol(mode.c_str(), &endptr, 10);
      if ((*endptr != '\0') || (compression_level < 0) || (compression_level > 9)) {
        PetscPrintf(com, "PISM WARNING: invalid compression level %s. Output compression is disabled.\n",
                    mode.c_str());
        compression_level = 0;
      }
    }

    return NCFile::Ptr(new NC4_Quilt(com, compression_level));
  }
#if (PISM_USE_PARALLEL_NETCDF4==1)
  else if (mode == "netcdf4_parallel") {
    return NCFile::Ptr(new NC4_Par(com));
  }
#endif
#if (PISM_USE_PNETCDF==1)
  else if (mode == "pnetcdf") {
    return NCFile::Ptr(new PNCFile(com));
  }
#endif
#if (PISM_USE_HDF5==1)
  else if (mode == "hdf5") {
    return NCFile::Ptr(new NC4_HDF5(com));
  }
#endif
  else {
    return NCFile::Ptr();       // a "NULL" pointer
  }
}

//! \brief The code shared by different PIO constructors.
void PIO::constructor(MPI_Comm c, const string &mode) {
  m_com  = c;
  m_mode = mode;
  m_nc   = create_backend(m_com, m_mode);

  if (mode != "guess_mode" && not m_nc) {
    throw RuntimeError("failed to allocate an I/O backend (class PIO)");
  }

}

PIO::PIO(MPI_Comm c, const string &mode, const UnitSystem &units_system)
  : m_unit_system(units_system) {
  constructor(c, mode);
}

PIO::PIO(IceGrid &grid, const string &mode)
  : m_unit_system(grid.get_unit_system()) {
  constructor(grid.com, mode);
  if (m_nc) {
    set_local_extent(grid.xs, grid.xm, grid.ys, grid.ym);
  }
}

PIO::PIO(const PIO &other)
  : m_unit_system(other.m_unit_system) {
  m_com  = other.m_com;
  m_nc   = other.m_nc;
  m_mode = other.m_mode;
}

PIO::~PIO() {
}

// Chooses the best I/O backend for reading from 'filename'.
void PIO::detect_mode(const string &filename) {
  assert(m_nc == nullptr);

  string format;
  {
    NC3File nc3(m_com);

    // detect_mode is private, so the caller will handle the failure
    // of open and add context
    nc3.open(filename, PISM_READONLY);
    format = nc3.get_format();
    nc3.close();
  }

  vector<string> modes;
  if (format == "netcdf4") {
    modes.push_back("netcdf4_parallel");
    modes.push_back("netcdf3");
  } else {
    modes.push_back("pnetcdf");
    modes.push_back("netcdf3");
  }

  for (unsigned int j = 0; j < modes.size(); ++j) {
    m_nc = create_backend(m_com, modes[j]);

    if (m_nc) {
      m_mode = modes[j];
      verbPrintf(3, m_com,
                 "  - Using the %s backend to read from %s...\n",
                 modes[j].c_str(), filename.c_str());
      break;
    }
  }

  if (not m_nc) {
    throw RuntimeError("failed to allocate an I/O backend (class PIO)");
  }

  m_nc->set_local_extent(m_xs, m_xm, m_ys, m_ym);
}


/**
 * Check if the storage order of a variable in the current file
 * matches the memory storage order used by PISM.
 *
 * @param var_name name of the variable to check
 * @param result set to false if storage orders match, true otherwise
 *
 * @return 0 on success
 */
void PIO::use_mapped_io(string var_name, bool &result) const {

  vector<string> dimnames = inq_vardims(var_name);

  vector<AxisType> storage, memory;
  memory.push_back(X_AXIS);
  memory.push_back(Y_AXIS);

  for (unsigned int j = 0; j < dimnames.size(); ++j) {
    AxisType dimtype = inq_dimtype(dimnames[j]);

    if (j == 0 && dimtype == T_AXIS) {
      // ignore the time dimension, but only if it is the first
      // dimension in the list
      continue;
    }

    if (dimtype == X_AXIS || dimtype == Y_AXIS) {
      storage.push_back(dimtype);
    } else if (dimtype == Z_AXIS) {
      memory.push_back(dimtype); // now memory = {X_AXIS, Y_AXIS, Z_AXIS}
      // assume that this variable has only one Z_AXIS in the file
      storage.push_back(dimtype);
    } else {
      // an UNKNOWN_AXIS or T_AXIS at index != 0 was found, use mapped I/O
      result = true;
      return;
    }
  }

  // we support 2D and 3D in-memory arrays, but not 4D
  assert(memory.size() <= 3);

  if (storage == memory) {
    // same storage order, do not use mapped I/O
    result = false;
  } else {
    // different storage orders, use mapped I/O
    result = true;
  }

  return;
}

bool PIO::check_if_exists(MPI_Comm com, const string &filename) {
  int file_exists = 0, rank = 0;
  MPI_Comm_rank(com, &rank);

  if (rank == 0) {
    // Check if the file exists:
    if (FILE *f = fopen(filename.c_str(), "r")) {
      file_exists = 1;
      fclose(f);
    } else {
      file_exists = 0;
    }
  }
  MPI_Bcast(&file_exists, 1, MPI_INT, 0, com);

  if (file_exists == 1) {
    return true;
  } else { 
    return false;
  }
}


void PIO::open(const string &filename, IO_Mode mode) {
  try {

    if (mode == PISM_READONLY || mode == PISM_READWRITE) {
      if (not m_nc && m_mode == "guess_mode") {
        detect_mode(filename);
      }
    }

    // opening for reading
    if (mode == PISM_READONLY) {

      assert(m_nc != nullptr);
      m_nc->open(filename, mode);

    } else if (mode == PISM_READWRITE_CLOBBER ||
               mode == PISM_READWRITE_MOVE) {

      assert(m_nc != nullptr);

      if (mode == PISM_READWRITE_MOVE) {
        m_nc->move_if_exists(filename);
      } else {
        m_nc->remove_if_exists(filename);
      }

      m_nc->create(filename);

      int old_fill;
      m_nc->set_fill(PISM_NOFILL, old_fill);
    } else {                      // mode == PISM_READWRITE
      assert(m_nc != nullptr);

      m_nc->open(filename, mode);

      int old_fill;
      m_nc->set_fill(PISM_NOFILL, old_fill);
    }
  } catch (RuntimeError &e) {
    e.add_context("opening or creating " + filename);
    throw;
  }
}


void PIO::close() {
  try {
    m_nc->close();
  } catch (RuntimeError &e) {
    e.add_context("closing " + inq_filename());
    throw;
  }
}

void PIO::redef() const {
  try {
    m_nc->redef();
  } catch (RuntimeError &e) {
    e.add_context("switching to define mode; file " + inq_filename());
    throw;
  }
}


void PIO::enddef() const {
  try {
    m_nc->enddef();
  } catch (RuntimeError &e) {
    e.add_context("switching to data mode; file " + inq_filename());
    throw;
  }
}

string PIO::inq_filename() const {
  return m_nc->get_filename();
}


//! \brief Get the number of records. Uses the length of an unlimited dimension.
unsigned int PIO::inq_nrecords() const {
  try {
    string dim;
    m_nc->inq_unlimdim(dim);

    if (dim.empty()) {
      return 1;                 // one record
    } else {
      return this->inq_dimlen(dim);
    }
  } catch (RuntimeError &e) {
    e.add_context("getting the number of records in file " + inq_filename());
  }
  return 0;                     // will never happen
}

//! \brief Get the number of records of a certain variable. Uses the length of
//! an associated "time" dimension.
unsigned int PIO::inq_nrecords(const string &name, const string &std_name) const {
  try {
    bool exists = false, found_by_standard_name = false;
    string name_found;
    inq_var(name, std_name, exists, name_found, found_by_standard_name);

    if (exists == false) {
      return 0;
    }

    vector<string> dims;
    m_nc->inq_vardimid(name_found, dims);

    for (unsigned int j = 0; j < dims.size(); ++j) {
      AxisType dimtype = inq_dimtype(dims[j]);

      if (dimtype == T_AXIS) {
        return this->inq_dimlen(dims[j]);
      }
    }

    return 1;                   // one record
  } catch (RuntimeError &e) {
    e.add_context("getting the number of records of variable '%s' ('%s') in '%s'", name.c_str(), std_name.c_str(), inq_filename().c_str());
    throw;
  }
  return 0;                     // will never happen
}


//! \brief Find a variable using its standard name and/or short name.
/*!
 * Sets "result" to the short name found.
 */
void PIO::inq_var(const string &short_name, const string &std_name, bool &exists,
                  string &result, bool &found_by_standard_name) const {
  try {
    exists = false;

    if (std_name.empty() == false) {

      int nvars;
      m_nc->inq_nvars(nvars);

      for (int j = 0; j < nvars; ++j) {
        string name;
        m_nc->inq_varname(j, name);

        string attribute = get_att_text(name, "standard_name");

        if (attribute.empty()) {
          continue;
        }

        if (attribute == std_name) {
          if (exists == false) {
            exists = true;
            found_by_standard_name = true;
            result = name;
          } else {
            throw RuntimeError::formatted("inconsistency in '%s': variables '%s' and '%s'\n"
                                          "have the same standard_name (%s)",
                                          inq_filename().c_str(), result.c_str(),
                                          name.c_str(), attribute.c_str());
          }
        }

      } // end of the for loop
    } // end of if (std_name.empty() == false)

    if (exists == false) {
      m_nc->inq_varid(short_name, exists);
      if (exists == true) {
        result = short_name;
      } else {
        result.clear();
      }

      found_by_standard_name = false;
    }

  } catch (RuntimeError &e) {
    e.add_context("searching for variable '%s' ('%s') in '%s'", short_name.c_str(), std_name.c_str(), inq_filename().c_str());
    throw;
  }
}

//! \brief Checks if a variable exists.
bool PIO::inq_var(const string &name) const {
  try {
    bool exists = false;
    m_nc->inq_varid(name, exists);
    return exists;
  } catch (RuntimeError &e) {
    e.add_context("searching for variable '%s' in '%s'", name.c_str(), inq_filename().c_str());
    throw;
  }
}

vector<string> PIO::inq_vardims(const string &name) const {
  try {
    vector<string> result;
    m_nc->inq_vardimid(name, result);
    return result;
  } catch (RuntimeError &e) {
    e.add_context("getting dimensions of variable '%s' in '%s'", name.c_str(), inq_filename().c_str());
    throw;
  }
}


//! \brief Checks if a dimension exists.
bool PIO::inq_dim(const string &name) const {
  try {
    bool exists = false;
    m_nc->inq_dimid(name, exists);
    return exists;
  } catch (RuntimeError &e) {
    e.add_context("searching for dimension '%s' in '%s'", name.c_str(), inq_filename().c_str());
    throw;
  }
}

//! \brief Get the length of a dimension.
/*!
 * Sets result to 0 if a dimension does not exist.
 */
unsigned int PIO::inq_dimlen(const string &name) const {
  try {
    bool exists = false;
    m_nc->inq_dimid(name, exists);
    if (exists == true) {
      unsigned int result = 0;
      m_nc->inq_dimlen(name, result);
      return result;
    } else {
      return 0;
    }
  } catch (RuntimeError &e) {
    e.add_context("getting the length of dimension '%s' in '%s'", name.c_str(), inq_filename().c_str());
    throw;
  }
}

//! \brief Get the "type" of a dimension.
/*!
 * The "type" is one of X_AXIS, Y_AXIS, Z_AXIS, T_AXIS.
 */
AxisType PIO::inq_dimtype(const string &name) const {
  try {
    string axis, standard_name, units;
    Unit tmp_units(m_unit_system, "1");
    bool exists;

    m_nc->inq_varid(name, exists);
    
    if (exists == false) {
      throw RuntimeError("coordinate variable " + name + " is missing");
    }

    axis          = get_att_text(name, "axis");
    standard_name = get_att_text(name, "standard_name");
    units         = get_att_text(name, "units");

    // check if it has units compatible with "seconds":

    tmp_units = Unit(m_unit_system, units);

    Unit seconds(m_unit_system, "seconds");
    if (UnitConverter::are_convertible(tmp_units, seconds)) {
      return T_AXIS;
    }

    // check the standard_name attribute:
    if (standard_name == "time") {
      return T_AXIS;
    } else if (standard_name == "projection_x_coordinate") {
      return X_AXIS;
    } else if (standard_name == "projection_y_coordinate") {
      return Y_AXIS;
    }

    // check the axis attribute:
    if (axis == "T" || axis == "t") {
      return T_AXIS;
    } else if (axis == "X" || axis == "x") {
      return X_AXIS;
    } else if (axis == "Y" || axis == "y") {
      return Y_AXIS;
    } else if (axis == "Z" || axis == "z") {
      return Z_AXIS;
    }

    // check the variable name:
    if (name == "x" || name == "X" ||
        name.find("x") == 0 || name.find("X") == 0) {
      return X_AXIS;
    } else if (name == "y" || name == "Y" ||
               name.find("y") == 0 || name.find("Y") == 0) {
      return Y_AXIS;
    } else if (name == "z" || name == "Z" ||
               name.find("z") == 0 || name.find("Z") == 0) {
      return Z_AXIS;
    } else if (name == "t" || name == "T" || name == "time" ||
               name.find("t") == 0 || name.find("T") == 0) {
      return T_AXIS;
    }

    // we have no clue:
    return UNKNOWN_AXIS;
  } catch (RuntimeError &e) {
    e.add_context("getting the type of dimension '%s' in '%s'", name.c_str(), inq_filename().c_str());
    throw;
  }
  return UNKNOWN_AXIS;          // will never happen
}

void PIO::inq_dim_limits(const string &name, double *min, double *max) const {
  try {
    vector<double> data;
    get_dim(name, data);

    double my_min = data[0],
      my_max = data[0];
    for (unsigned int j = 0; j < data.size(); ++j) {
      my_min = PetscMin(data[j], my_min);
      my_max = PetscMax(data[j], my_max);
    }

    if (min != NULL) {
      *min = my_min;
    }

    if (max != NULL) {
      *max = my_max;
    }

  } catch (RuntimeError &e) {
    e.add_context("getting limits of dimension '%s' in '%s'", name.c_str(), inq_filename().c_str());
    throw;
  }
}


//! \brief Sets grid parameters using data read from the file.
PetscErrorCode PIO::inq_grid(const string &var_name, IceGrid *grid, Periodicity periodicity) const {
  try {
    PetscErrorCode ierr;

    assert(grid != NULL);

    // The following call may fail because var_name does not exist. (And this is fatal!)
    grid_info input = inq_grid_info(var_name, periodicity);

    // if we have no vertical grid information, create a fake 2-level vertical grid.
    if (input.z.size() < 2) {
      double Lz = grid->config.get("grid_Lz");
      verbPrintf(3, m_com,
                 "WARNING: Can't determine vertical grid information using '%s' in %s'\n"
                 "         Using 2 levels and Lz of %3.3fm\n",
                 var_name.c_str(), inq_filename().c_str(), Lz);

      input.z.clear();
      input.z.push_back(0);
      input.z.push_back(Lz);
    }

    grid->Mx = input.x_len;
    grid->My = input.y_len;

    grid->periodicity = periodicity;

    grid->x0 = input.x0;
    grid->y0 = input.y0;
    grid->Lx = input.Lx;
    grid->Ly = input.Ly;

    grid->time->set_start(input.time);
    ierr = grid->time->init(); CHKERRQ(ierr); // re-initialize to take the new start time into account

    ierr = grid->compute_horizontal_spacing(); CHKERRQ(ierr);
    ierr = grid->set_vertical_levels(input.z); CHKERRQ(ierr);

    // We're ready to call grid->allocate().
  } catch (RuntimeError &e) {
    e.add_context("initializing computational grid from " + inq_filename());
    throw;
  }
  return 0;
}

/*! Do not use this method to get units of time and time_bounds
  variables: in these two cases we need to handle the reference date
  correctly.
 */
void PIO::inq_units(const string &name, bool &has_units, Unit &units) const {
  try {
    // Get the string:
    string units_string = get_att_text(name, "units");

    // If a variables does not have the units attribute, set the flag and return:
    if (units_string.empty()) {
      has_units = false;
      units = Unit(units.get_system(), "1");
      return;
    }

    // strip trailing spaces
    while (ends_with(units_string, " ")) {
      units_string.resize(units_string.size() - 1);
    }

    // this may fail if the unit string is invalid
    units = Unit(units.get_system(), units_string);

    has_units = true;
  } catch (RuntimeError &e) {
    e.add_context("getting units of variable '%s' in '%s'", name.c_str(), inq_filename().c_str());
    throw;
  }
}


grid_info PIO::inq_grid_info(const string &name, Periodicity p) const {
  try {
    vector<string> dims;
    bool exists, found_by_standard_name;
    string name_found;
    grid_info result;

    // try "name" as the standard_name first, then as the short name:
    inq_var(name, name, exists, name_found, found_by_standard_name);

    if (exists == false) {
      throw RuntimeError("variable " + name + " is missing");
    }

    m_nc->inq_vardimid(name_found, dims);

    // use "global" dimensions (as opposed to dimensions of a patch)
    if (m_mode == "quilt") {
      for (unsigned int i = 0; i < dims.size(); ++i) {
        if (dims[i] == "x_patch") {
          dims[i] = "x";
        }

        if (dims[i] == "y_patch") {
          dims[i] = "y";
        }
      }
    }

    for (unsigned int i = 0; i < dims.size(); ++i) {
      string dimname = dims[i];

      AxisType dimtype = inq_dimtype(dimname);

      switch (dimtype) {
      case X_AXIS:
        {
          m_nc->inq_dimlen(dimname, result.x_len);
          double x_min = 0.0, x_max = 0.0;
          inq_dim_limits(dimname, &x_min, &x_max);
          get_dim(dimname, result.x);
          result.x0 = 0.5 * (x_min + x_max);
          result.Lx = 0.5 * (x_max - x_min);
          if (p & X_PERIODIC) {
            const double dx = result.x[1] - result.x[0];
            result.Lx += 0.5 * dx;
          }
          break;
        }
      case Y_AXIS:
        {
          m_nc->inq_dimlen(dimname, result.y_len);
          double y_min = 0.0, y_max = 0.0;
          inq_dim_limits(dimname, &y_min, &y_max);
          get_dim(dimname, result.y);
          result.y0 = 0.5 * (y_min + y_max);
          result.Ly = 0.5 * (y_max - y_min);
          if (p & Y_PERIODIC) {
            const double dy = result.y[1] - result.y[0];
            result.Ly += 0.5 * dy;
          }
          break;
        }
      case Z_AXIS:
        {
          m_nc->inq_dimlen(dimname, result.z_len);
          inq_dim_limits(dimname, &result.z_min, &result.z_max);
          get_dim(dimname, result.z);
          break;
        }
      case T_AXIS:
        {
          m_nc->inq_dimlen(dimname, result.t_len);
          inq_dim_limits(dimname, NULL, &result.time);
          break;
        }
      default:
        {
          throw RuntimeError::formatted("can't figure out which direction dimension '%s' corresponds to.",
                                        dimname.c_str());
        }
      } // switch
    }   // for loop

    return result;

  } catch (RuntimeError &e) {
    e.add_context("getting grid information using variable '%s' in '%s'", name.c_str(), inq_filename().c_str());
    throw;
  }
}

//! \brief Define a dimension \b and the associated coordinate variable. Set attributes.
void PIO::def_dim(unsigned long int length, const NCVariable &metadata) const {
  string name = metadata.get_name();
  try {
    m_nc->redef();

    m_nc->def_dim(name, length);

    vector<string> dims(1, name);
    m_nc->def_var(name, PISM_DOUBLE, dims);

    write_attributes(metadata, PISM_DOUBLE, false);
  } catch (RuntimeError &e) {
    e.add_context("defining dimension '%s' in '%s'", name.c_str(), inq_filename().c_str());
    throw;
  }
}

//! \brief Define a variable.
void PIO::def_var(const string &name, IO_Type nctype, const vector<string> &dims) const {
  try {
    m_nc->def_var(name, nctype, dims);
  } catch (RuntimeError &e) {
    e.add_context("defining variable '%s' in '%s'", name.c_str(), inq_filename().c_str());
    throw;
  }
}

void PIO::get_1d_var(const string &name, unsigned int s, unsigned int c,
                     vector<double> &result) const {
  vector<unsigned int> start(1), count(1);

  result.resize(c);

  start[0] = s;
  count[0] = c;

  // this call will handle errors
  get_vara_double(name, start, count, &result[0]);
}


void PIO::put_1d_var(const string &name, unsigned int s, unsigned int c,
                     const vector<double> &data) const {
  vector<unsigned int> start(1), count(1);

  start[0] = s;
  count[0] = c;

  // this call will handle errors
  put_vara_double(name, start, count, const_cast<double*>(&data[0]));
}


//! \brief Get dimension data (a coordinate variable).
void PIO::get_dim(const string &name, vector<double> &data) const {
  try {
    unsigned int dim_length = 0;
    m_nc->inq_dimlen(name, dim_length);

    get_1d_var(name, 0, dim_length, data);
  } catch (RuntimeError &e) {
    e.add_context("reading dimension '%s' from '%s'", name.c_str(), inq_filename().c_str());
    throw;
  }
}

//! \brief Write dimension data (a coordinate variable).
void PIO::put_dim(const string &name, const vector<double> &data) const {
  try {
    put_1d_var(name, 0, (unsigned int)data.size(), data);
  } catch (RuntimeError &e) {
    e.add_context("writing dimension '%s' to '%s'", name.c_str(), inq_filename().c_str());
    throw;
  }
}

void PIO::def_time(const string &name, const string &calendar, const string &units) const {
  try {
    std::map<string,string> attrs;

    bool time_exists;
    m_nc->inq_varid(name, time_exists);
    if (time_exists) {
      return;
    }

    // time
    NCVariable time(name, m_unit_system);
    time.set_string("long_name", "time");
    time.set_string("calendar", calendar);
    time.set_units(units);
    time.set_string("axis", "T");

    def_dim(PISM_UNLIMITED, time);
  } catch (RuntimeError &e) {
    e.add_context("defining the time dimension in " + inq_filename());
    throw;
  }
}


//! \brief Append to the time dimension.
void PIO::append_time(const string &name, double value) const {
  try {
    vector<unsigned int> start(1), count(1);

    unsigned int dim_length = 0;
    m_nc->inq_dimlen(name, dim_length);

    start[0] = dim_length;
    count[0] = 1;

    put_vara_double(name, start, count, &value);
  } catch (RuntimeError &e) {
    e.add_context("appending to the time dimension in " + inq_filename());
    throw;
  }
}

//! \brief Append to the history global attribute.
/*!
 * Use put_att_text("PISM_GLOBAL", "history", ...) to overwrite "history".
 */
void PIO::append_history(const string &history) const {
  try {
    string old_history = get_att_text("PISM_GLOBAL", "history");
    put_att_text("PISM_GLOBAL", "history", history + old_history);
  } catch (RuntimeError &e) {
    e.add_context("appending to the history attribute in " + inq_filename());
    throw;
  }
}

//! \brief Write a multiple-valued double attribute.
void PIO::put_att_double(const string &var_name, const string &att_name, IO_Type nctype,
                         const vector<double> &values) const {
  try {
    m_nc->redef();
    m_nc->put_att_double(var_name, att_name, nctype, values);
  } catch (RuntimeError &e) {
    e.add_context("writing double attribute '%s:%s' in '%s'",
                  var_name.c_str(), att_name.c_str(), inq_filename().c_str());
    throw;
  }
}

//! \brief Write a single-valued double attribute.
void PIO::put_att_double(const string &var_name, const string &att_name, IO_Type nctype,
                         double value) const {
  vector<double> tmp(1, value);

  // this call will handle errors
  this->put_att_double(var_name, att_name, nctype, tmp);
}

//! \brief Write a text attribute.
void PIO::put_att_text(const string &var_name, const string &att_name,
                       const string &value) const {
  try {
    m_nc->redef();

    string tmp = value + "\0";    // ensure that the string is null-terminated

    m_nc->put_att_text(var_name, att_name, tmp);
  } catch (RuntimeError &e) {
    e.add_context("writing text attribute '%s:%s' in '%s'",
                  var_name.c_str(), att_name.c_str(), inq_filename().c_str());
    throw;
  }
}

//! \brief Get a double attribute.
vector<double> PIO::get_att_double(const string &var_name, const string &att_name) const {
  try {
    IO_Type att_type;
    m_nc->inq_atttype(var_name, att_name, att_type);

    // Give an understandable error message if a string attribute was found when
    // a number (or a list of numbers) was expected. (We've seen datasets with
    // "valid_min" stored as a string...)
    if (att_type == PISM_CHAR) {
      string tmp = get_att_text(var_name, att_name);

      throw RuntimeError::formatted("attribute %s is a string '%s'; expected a number or a list of numbers",
                                    att_name.c_str(), tmp.c_str());
    } else {
      // In this case att_type might be PISM_NAT (if an attribute does not
      // exist), but get_att_double can handle that.
      vector<double> result;
      m_nc->get_att_double(var_name, att_name, result);
      return result;
    }
  } catch (RuntimeError &e) {
    e.add_context("reading double attribute '%s:%s' from '%s'",
                  var_name.c_str(), att_name.c_str(), inq_filename().c_str());
    throw;
  }
}

//! \brief Get a text attribute.
string PIO::get_att_text(const string &var_name, const string &att_name) const {
  try {
    string result;
    m_nc->get_att_text(var_name, att_name, result);
    return result;
  } catch (RuntimeError &e) {
    e.add_context("reading text attribute '%s:%s' from %s", var_name.c_str(), att_name.c_str(), inq_filename().c_str());
    throw;
  }
}

//! \brief Read a PETSc Vec using the grid "grid".
/*!
 * Assumes that double corresponds to C++ double.
 *
 * Vec result has to be "global" (i.e. without ghosts).
 */
PetscErrorCode PIO::get_vec(IceGrid *grid, const string &var_name,
                            unsigned int z_count, unsigned int t_start, Vec result) const {
  try {
    PetscErrorCode ierr;

    vector<unsigned int> start, count, imap;
    const unsigned int t_count = 1;
    compute_start_and_count(var_name,
                            t_start, t_count,
                            grid->xs, grid->xm,
                            grid->ys, grid->ym,
                            0, z_count,
                            start, count, imap);

    double *a_petsc;
    ierr = VecGetArray(result, &a_petsc);
    PISM_PETSC_CHK(ierr, "VecGetArray");

    bool mapped_io = true;
    use_mapped_io(var_name, mapped_io);
    if (mapped_io == true) {
      get_varm_double(var_name, start, count, imap, (double*)a_petsc);
    } else {
      get_vara_double(var_name, start, count, (double*)a_petsc);
    }

    ierr = VecRestoreArray(result, &a_petsc);
    PISM_PETSC_CHK(ierr, "VecRestoreArray");

  } catch (RuntimeError &e) {
    e.add_context("reading variable '%s' from '%s'", var_name.c_str(), inq_filename().c_str());
    throw;
  }
  return 0;
}

unsigned int PIO::inq_nattrs(const string &var_name) const {
  try {
    int result = 0;
    m_nc->inq_varnatts(var_name, result);
    return result;
  } catch (RuntimeError &e) {
    e.add_context("getting the number of attributes of variable '%s' in '%s'", var_name.c_str(), inq_filename().c_str());
    throw;
  }
}


string PIO::inq_attname(const string &var_name, unsigned int n) const {
  try {
    string result;
    m_nc->inq_attname(var_name, n, result);
    return result;
  } catch (RuntimeError &e) {
    e.add_context("getting the name of an attribute of variable '%s' in '%s'", var_name.c_str(), inq_filename().c_str());
    throw;
  }
}


IO_Type PIO::inq_atttype(const string &var_name, const string &att_name) const {
  try {
    IO_Type result;
    m_nc->inq_atttype(var_name, att_name, result);
    return result;
  } catch (RuntimeError &e) {
    e.add_context("getting the type of an attribute of variable '%s' in '%s'", var_name.c_str(), inq_filename().c_str());
    throw;
  }
}


//! \brief Write a PETSc Vec using the grid "grid".
/*!
 * Assumes that double corresponds to C++ double.
 *
 * Vec input has to be "global" (i.e. without ghosts).
 *
 * This method always writes to the last record in the file.
 */
PetscErrorCode PIO::put_vec(IceGrid *grid, const string &var_name, unsigned int z_count, Vec input) const {
  try {
    PetscErrorCode ierr;

    unsigned int t_length;
    m_nc->inq_dimlen(grid->config.get_string("time_dimension_name"), t_length);

    assert(t_length >= 1);

    vector<unsigned int> start, count, imap;
    const unsigned int t_count = 1;
    compute_start_and_count(var_name,
                            t_length - 1, t_count,
                            grid->xs, grid->xm,
                            grid->ys, grid->ym,
                            0, z_count,
                            start, count, imap);

    double *a_petsc;
    ierr = VecGetArray(input, &a_petsc);
    PISM_PETSC_CHK(ierr, "VecGetArray");

    if (grid->config.get_string("output_variable_order") == "xyz") {
      // Use the faster and safer (avoids a NetCDF bug) call if the aray storage
      // orders in the memory and in NetCDF files are the same.
      put_vara_double(var_name, start, count, (double*)a_petsc);
    } else {
      // Revert to "mapped" I/O otherwise.
      put_varm_double(var_name, start, count, imap, (double*)a_petsc);
    }

    ierr = VecRestoreArray(input, &a_petsc);
    PISM_PETSC_CHK(ierr, "VecRestoreArray");
  } catch (RuntimeError &e) {
    e.add_context("writing variable '%s' to '%s'", var_name.c_str(), inq_filename().c_str());
    throw;
  }
  return 0;
}

//! \brief Get the interpolation context (grid information) for an input file.
/*!
 * Sets lic to NULL if the variable was not found.
 *
 * @note The *caller* is in charge of destroying lic
 */
LocalInterpCtx* PIO::get_interp_context(const string &name,
                                        const IceGrid &grid,
                                        const vector<double> &zlevels) const {
  bool exists = inq_var(name);

  if (exists == false) {
    throw RuntimeError("variable " + name + " is missing in " + inq_filename());
  } else {
    grid_info gi = inq_grid_info(name, grid.periodicity);

    return new LocalInterpCtx(gi, grid, zlevels.front(), zlevels.back());
  }
}

//! \brief Read a PETSc Vec from a file, using bilinear (or trilinear)
//! interpolation to put it on the grid defined by "grid" and zlevels_out.
void PIO::regrid_vec(IceGrid *grid, const string &var_name,
                     const vector<double> &zlevels_out,
                     unsigned int t_start, Vec result) const {
  try {
    const int X = 1, Y = 2, Z = 3; // indices, just for clarity
    vector<unsigned int> start, count, imap;

    shared_ptr<LocalInterpCtx> lic(get_interp_context(var_name, *grid, zlevels_out));
    assert(lic != nullptr);

    const unsigned int t_count = 1;
    compute_start_and_count(var_name,
                            t_start, t_count,
                            lic->start[X], lic->count[X],
                            lic->start[Y], lic->count[Y],
                            lic->start[Z], lic->count[Z],
                            start, count, imap);

    bool mapped_io = true;
    use_mapped_io(var_name, mapped_io);
    if (mapped_io == true) {
      get_varm_double(var_name, start, count, imap, lic->a);
    } else {
      get_vara_double(var_name, start, count, lic->a);
    }

    regrid(grid, zlevels_out, lic.get(), result);

  } catch (RuntimeError &e) {
    e.add_context("reading variable '%s' (using linear interpolation) from '%s'",
                  var_name.c_str(), inq_filename().c_str());
    throw;
  }
}

/** Regrid `var_name` from a file, replacing missing values with `default_value`.
 *
 * @param grid computational grid; used to initialize interpolation
 * @param var_name variable to regrid
 * @param zlevels_out vertical levels of the resulting grid
 * @param t_start time index of the record to regrid
 * @param default_value default value to replace `_FillValue` with
 * @param[out] result resulting interpolated field
 */
void PIO::regrid_vec_fill_missing(IceGrid *grid, const string &var_name,
                                  const vector<double> &zlevels_out,
                                  unsigned int t_start,
                                  double default_value,
                                  Vec result) const {
  try {
    const int X = 1, Y = 2, Z = 3; // indices, just for clarity
    vector<unsigned int> start, count, imap;

    shared_ptr<LocalInterpCtx> lic(get_interp_context(var_name, *grid, zlevels_out));
    assert(lic != nullptr);

    const unsigned int t_count = 1;
    compute_start_and_count(var_name,
                            t_start, t_count,
                            lic->start[X], lic->count[X],
                            lic->start[Y], lic->count[Y],
                            lic->start[Z], lic->count[Z],
                            start, count, imap);

    bool mapped_io = true;
    use_mapped_io(var_name, mapped_io);
    if (mapped_io == true) {
      get_varm_double(var_name, start, count, imap, lic->a);
    } else {
      get_vara_double(var_name, start, count, lic->a);
    }

    // Replace missing values if the _FillValue attribute is present,
    // and if we have missing values to replace.
    {
      vector<double> attribute;
      m_nc->get_att_double(var_name, "_FillValue", attribute);
      if (attribute.size() == 1) {
        const double fill_value = attribute[0],
          epsilon = 1e-12;
        for (unsigned int i = 0; i < lic->a_len; ++i) {
          if (fabs(lic->a[i] - fill_value) < epsilon) {
            lic->a[i] = default_value;
          }
        }
      }
    }

    regrid(grid, zlevels_out, lic.get(), result);
  } catch (RuntimeError &e) {
    e.add_context("reading variable '%s' (using linear interpolation) from '%s'",
                  var_name.c_str(), inq_filename().c_str());
    throw;
  }
}

int PIO::k_below(double z, const vector<double> &zlevels) const {
  double z_min = zlevels.front(), z_max = zlevels.back();
  int mcurr = 0;

  if (z < z_min - 1.0e-6 || z > z_max + 1.0e-6) {
    throw RuntimeError::formatted("PIO::k_below(): z = %5.4f is outside the allowed range: [%f, %f]",
                                  z, z_min, z_max);
  }

  while (zlevels[mcurr+1] < z) {
    mcurr++;
  }

  return mcurr;
}


//! \brief Bi-(or tri-)linear interpolation.
/*!
 * This is the interpolation code itself.
 *
 * Note that its inputs are (essentially)
 * - the definition of the input grid
 * - the definition of the output grid
 * - input array (lic->a)
 * - output array (Vec result)
 *
 * We should be able to switch to using an external interpolation library
 * fairly easily...
 */
PetscErrorCode PIO::regrid(IceGrid *grid, const vector<double> &zlevels_out,
                           LocalInterpCtx *lic, Vec result) const {
  const int Y = 2, Z = 3; // indices, just for clarity
  PetscErrorCode ierr;

  vector<double> &zlevels_in = lic->zlevels;
  unsigned int nlevels = zlevels_out.size();
  double *input_array = lic->a;

  // array sizes for mapping from logical to "flat" indices
  int y_count = lic->count[Y],
    z_count = lic->count[Z];

  // We'll work with the raw storage here so that the array we are filling is
  // indexed the same way as the buffer we are pulling from (input_array)
  double *output_array;
  ierr = VecGetArray(result, &output_array);
  PISM_PETSC_CHK(ierr, "VecGetArray");

  // NOTE: make sure that the traversal order is correct!
  for (Points p(*grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    const int i0 = i - grid->xs, j0 = j - grid->ys;

    for (unsigned int k = 0; k < nlevels; k++) {
      // location (x,y,z) is in target computational domain
      const double
        z = zlevels_out[k];

      // Indices of neighboring points.
      const int
        Im = lic->x_left[i0],
        Ip = lic->x_right[i0],
        Jm = lic->y_left[j0],
        Jp = lic->y_right[j0];

      double a_mm, a_mp, a_pm, a_pp;  // filled differently in 2d and 3d cases

      if (nlevels > 1) {
        // get the index into the source grid, for just below the level z
        const int kc = k_below(z, zlevels_in);

        // We pretend that there are always 8 neighbors (4 in the map plane,
        // 2 vertical levels). And compute the indices into the input_array for
        // those neighbors.
        const int mmm = (Im * y_count + Jm) * z_count + kc;
        const int mmp = (Im * y_count + Jm) * z_count + kc + 1;
        const int mpm = (Im * y_count + Jp) * z_count + kc;
        const int mpp = (Im * y_count + Jp) * z_count + kc + 1;
        const int pmm = (Ip * y_count + Jm) * z_count + kc;
        const int pmp = (Ip * y_count + Jm) * z_count + kc + 1;
        const int ppm = (Ip * y_count + Jp) * z_count + kc;
        const int ppp = (Ip * y_count + Jp) * z_count + kc + 1;

        // We know how to index the neighbors, but we don't yet know where the
        // point lies within this box.  This is represented by kk in [0,1].
        const double zkc = zlevels_in[kc];
        double dz;
        if (kc == z_count - 1) {
          dz = zlevels_in[kc] - zlevels_in[kc-1];
        } else {
          dz = zlevels_in[kc+1] - zlevels_in[kc];
        }
        const double kk = (z - zkc) / dz;

        // linear interpolation in the z-direction
        a_mm = input_array[mmm] * (1.0 - kk) + input_array[mmp] * kk;
        a_mp = input_array[mpm] * (1.0 - kk) + input_array[mpp] * kk;
        a_pm = input_array[pmm] * (1.0 - kk) + input_array[pmp] * kk;
        a_pp = input_array[ppm] * (1.0 - kk) + input_array[ppp] * kk;
      } else {
        // we don't need to interpolate vertically for the 2-D case
        a_mm = input_array[Im * y_count + Jm];
        a_mp = input_array[Im * y_count + Jp];
        a_pm = input_array[Ip * y_count + Jm];
        a_pp = input_array[Ip * y_count + Jp];
      }

      // interpolation coefficient in the y direction
      const double jj = lic->y_alpha[j0];

      // interpolate in y direction
      const double a_m = a_mm * (1.0 - jj) + a_mp * jj;
      const double a_p = a_pm * (1.0 - jj) + a_pp * jj;

      // interpolation coefficient in the x direction
      const double ii = lic->x_alpha[i0];

      int index = (i0 * grid->ym + j0) * nlevels + k;

      // index into the new array and interpolate in x direction
      output_array[index] = a_m * (1.0 - ii) + a_p * ii;
      // done with the point at (x,y,z)
    }
  }

  ierr = VecRestoreArray(result, &output_array);
  PISM_PETSC_CHK(ierr, "VecRestoreArray");

  return 0;
}


void PIO::compute_start_and_count(const string &short_name,
                                  unsigned int t_start, unsigned int t_count,
                                  unsigned int x_start, unsigned int x_count,
                                  unsigned int y_start, unsigned int y_count,
                                  unsigned int z_start, unsigned int z_count,
                                  vector<unsigned int> &start,
                                  vector<unsigned int> &count,
                                  vector<unsigned int> &imap) const {
  vector<string> dims;

  m_nc->inq_vardimid(short_name, dims);
  unsigned int ndims = dims.size();

  // Resize output vectors:
  start.resize(ndims);
  count.resize(ndims);
  imap.resize(ndims);

  // Assemble start, count and imap:
  for (unsigned int j = 0; j < ndims; j++) {
    string dimname = dims[j];

    AxisType dimtype = inq_dimtype(dimname);

    switch (dimtype) {
    case T_AXIS:
      start[j] = t_start;
      count[j] = t_count;
      imap[j]  = x_count * y_count * z_count;
      break;
    case X_AXIS:
      start[j] = x_start;
      count[j] = x_count;
      imap[j]  = y_count * z_count;
      break;
    case Y_AXIS:
      start[j] = y_start;
      count[j] = y_count;
      imap[j]  = z_count;
      break;
    default:
    case Z_AXIS:
      start[j] = z_start;
      count[j] = z_count;
      imap[j]  = 1;
      break;
    }

    // #if (PISM_DEBUG==1)
    //     fprintf(stderr, "[%d] var=%s start[%d]=%d count[%d]=%d imap[%d]=%d\n",
    //             rank, short_name.c_str(),
    //             j, start[j], j, count[j], j, imap[j]);
    // #endif

  }
}

void PIO::get_vara_double(const string &variable_name,
                          const vector<unsigned int> &start,
                          const vector<unsigned int> &count,
                          double *ip) const {
  try {
    m_nc->enddef();
    m_nc->get_vara_double(variable_name, start, count, ip);
  } catch (RuntimeError &e) {
    e.add_context("reading variable '%s' from '%s'", variable_name.c_str(), inq_filename().c_str());
    throw;
  }
}


void PIO::put_vara_double(const string &variable_name,
                          const vector<unsigned int> &start,
                          const vector<unsigned int> &count,
                          double *op) const {
  try {
    m_nc->enddef();
    m_nc->put_vara_double(variable_name, start, count, op);
  } catch (RuntimeError &e) {
    e.add_context("writing variable '%s' to '%s'", variable_name.c_str(), inq_filename().c_str());
    throw;
  }
}

void PIO::get_varm_double(const string &variable_name,
                          const vector<unsigned int> &start,
                          const vector<unsigned int> &count,
                          const vector<unsigned int> &imap, double *ip) const {
  try {
    m_nc->enddef();
    m_nc->get_varm_double(variable_name, start, count, imap, ip);
  } catch (RuntimeError &e) {
    e.add_context("reading variable '%s' from '%s'", variable_name.c_str(), inq_filename().c_str());
    throw;
  }
}


void PIO::put_varm_double(const string &variable_name,
                          const vector<unsigned int> &start,
                          const vector<unsigned int> &count,
                          const vector<unsigned int> &imap, double *op) const {
  try {
    m_nc->enddef();
    m_nc->put_varm_double(variable_name, start, count, imap, op);
  } catch (RuntimeError &e) {
    e.add_context("writing variable '%s' to '%s'", variable_name.c_str(), inq_filename().c_str());
    throw;
  }
}

void PIO::set_local_extent(unsigned int xs, unsigned int xm,
                           unsigned int ys, unsigned int ym) {
  m_nc->set_local_extent(xs, xm, ys, ym);
  m_xs = xs;
  m_xm = xm;
  m_ys = ys;
  m_ym = ym;
}

void PIO::read_attributes(const string &name, NCVariable &variable) const {
  try {
    bool variable_exists = inq_var(name);

    if (variable_exists == false) {
      throw RuntimeError("variable " + name + " is missing");
    }

    variable.clear_all_strings();
    variable.clear_all_doubles();

    unsigned int nattrs = inq_nattrs(name);

    for (unsigned int j = 0; j < nattrs; ++j) {
      string attribute_name = inq_attname(name, j);
      IO_Type nctype = inq_atttype(name, attribute_name);

      if (nctype == PISM_CHAR) {
        string value = get_att_text(name, attribute_name);

        if (attribute_name == "units") {
          variable.set_units(value);
        } else {
          variable.set_string(attribute_name, value);
        }
      } else {
        vector<double> values = get_att_double(name, attribute_name);
        variable.set_doubles(attribute_name, values);
      }
    } // end of for (int j = 0; j < nattrs; ++j)
  } catch (RuntimeError &e) {
    e.add_context("reading attributes of variable '%s' from '%s'", name.c_str(), inq_filename().c_str());
    throw;
  }
}


//! Write variable attributes to a NetCDF file.
/*!

  \li if write_in_glaciological_units == true, "glaciological_units" are
  written under the name "units" plus the valid range is written in
  glaciological units.

  \li if both valid_min and valid_max are set, then valid_range is written
  instead of the valid_min, valid_max pair.
 */
void PIO::write_attributes(const NCVariable &var, IO_Type nctype,
                           bool write_in_glaciological_units) const {
  string var_name = var.get_name();
  try {
    // units, valid_min, valid_max and valid_range need special treatment:
    if (var.has_attribute("units")) {
      string output_units = var.get_string("units");

      if (write_in_glaciological_units) {
        output_units = var.get_string("glaciological_units");
      }

      put_att_text(var_name, "units", output_units);
    }

    vector<double> bounds(2);
    double fill_value = 0.0;

    if (var.has_attribute("_FillValue")) {
      fill_value = var.get_double("_FillValue");
    }

    // We need to save valid_min, valid_max and valid_range in the units
    // matching the ones in the output.
    if (write_in_glaciological_units) {

      assert(UnitConverter::are_convertible(var.get_units(),
                                            var.get_glaciological_units()) == true);
      UnitConverter c(var.get_units(), var.get_glaciological_units());

      bounds[0]  = c(var.get_double("valid_min"));
      bounds[1]  = c(var.get_double("valid_max"));
      fill_value = c(fill_value);
    } else {
      bounds[0] = var.get_double("valid_min");
      bounds[1] = var.get_double("valid_max");
    }

    if (var.has_attribute("_FillValue")) {
      put_att_double(var_name, "_FillValue", nctype, fill_value);
    }

    if (var.has_attribute("valid_min") && var.has_attribute("valid_max")) {
      put_att_double(var_name, "valid_range", nctype, bounds);
    } else if (var.has_attribute("valid_min")) {
      put_att_double(var_name, "valid_min",   nctype, bounds[0]);
    } else if (var.has_attribute("valid_max")) {
      put_att_double(var_name, "valid_max",   nctype, bounds[1]);
    }

    // Write text attributes:
    const NCVariable::StringAttrs &strings = var.get_all_strings();
    NCVariable::StringAttrs::const_iterator i;
    for (i = strings.begin(); i != strings.end(); ++i) {
      string
        name  = i->first,
        value = i->second;

      if (name == "units" || name == "glaciological_units" || value.empty()) {
        continue;
      }

      put_att_text(var_name, name, value);
    }

    // Write double attributes:
    const NCVariable::DoubleAttrs &doubles = var.get_all_doubles();
    NCVariable::DoubleAttrs::const_iterator j;
    for (j = doubles.begin(); j != doubles.end(); ++j) {
      string name  = j->first;
      vector<double> values = j->second;

      if (name == "valid_min" ||
          name == "valid_max" ||
          name == "valid_range" ||
          name == "_FillValue" ||
          values.empty()) {
        continue;
      }

      put_att_double(var_name, name, nctype, values);
    }

  } catch (RuntimeError &e) {
    e.add_context("writing attributes of variable '%s' to '%s'", var_name.c_str(), inq_filename().c_str());
    throw;
  }
}

/** Write global attributes to a file.
 *
 * Same as `PIO::write_attributes(var, PISM_DOUBLE, false)`, but
 * prepends the history string.
 *
 * @param var metadata object containing attributes
 *
 * @return 0 on success
 */
void PIO::write_global_attributes(const NCVariable &var) const {
  try {
    NCVariable tmp = var;

    string old_history = get_att_text("PISM_GLOBAL", "history");

    tmp.set_name("PISM_GLOBAL");
    tmp.set_string("history", tmp.get_string("history") + old_history);

    write_attributes(tmp, PISM_DOUBLE, false);

  } catch (RuntimeError &e) {
    e.add_context("writing global attributes to " + inq_filename());
    throw;
  }
}

//! Read the valid range information from a file.
/*! Reads `valid_min`, `valid_max` and `valid_range` attributes; if \c
    valid_range is found, sets the pair `valid_min` and `valid_max` instead.
 */
void PIO::read_valid_range(const string &name, NCVariable &variable) const {
  try {
    // Never reset valid_min/max if they were set internally
    if (variable.has_attribute("valid_min") ||
        variable.has_attribute("valid_max")) {
      return;
    }

    // Read the units.
    string input_units_string = get_att_text(name, "units");

    UnitSystem sys = variable.get_units().get_system();
    Unit input_units = Unit(sys, input_units_string);

    UnitConverter c(input_units, variable.get_units());

    vector<double> bounds = get_att_double(name, "valid_range");
    if (bounds.size() == 2) {             // valid_range is present
      variable.set_double("valid_min", c(bounds[0]));
      variable.set_double("valid_max", c(bounds[1]));
    } else {                      // valid_range has the wrong length or is missing
      bounds = get_att_double(name, "valid_min");
      if (bounds.size() == 1) {           // valid_min is present
        variable.set_double("valid_min", c(bounds[0]));
      }

      bounds = get_att_double(name, "valid_max");
      if (bounds.size() == 1) {           // valid_max is present
        variable.set_double("valid_max", c(bounds[0]));
      }
    }
  } catch (RuntimeError &e) {
    e.add_context("reading valid range of variable '%s' from '%s'", name.c_str(), inq_filename().c_str());
    throw;
  }
}

//! Read a time-series variable from a NetCDF file to a vector of doubles.
void PIO::read_timeseries(const NCTimeseries &metadata,
                          Time *time, vector<double> &data) const {
  string name = metadata.get_name();

  try {
    bool variable_exists = false;

    // Find the variable:
    string name_found,
      long_name      = metadata.get_string("long_name"),
      standard_name  = metadata.get_string("standard_name"),
      dimension_name = metadata.get_dimension_name();

    bool found_by_standard_name = false;
    inq_var(name, standard_name, variable_exists, name_found, found_by_standard_name);

    if (not variable_exists) {
      throw RuntimeError("variable " + name + " is missing");
    }

    vector<string> dims = inq_vardims(name_found);
    if (dims.size() != 1) {
      throw RuntimeError("a time-series variable has to be one-dimensional");
    }

    unsigned int length = inq_dimlen(dimension_name);
    if (length <= 0) {
      throw RuntimeError("dimension " + dimension_name + " has length zero");
    }

    data.resize(length);          // memory allocation happens here

    get_1d_var(name_found, 0, length, data);

    bool input_has_units = false;
    Unit internal_units = metadata.get_units(),
      input_units(internal_units.get_system(), "1");

    string input_units_string = get_att_text(name_found, "units");

    if (input_units_string.empty() == true) {
      input_has_units = false;
    } else {
      input_units_string = time->CF_units_to_PISM_units(input_units_string);

      input_units = Unit(internal_units.get_system(), input_units_string);
      input_has_units = true;
    }

    if (metadata.has_attribute("units") == true && input_has_units == false) {
      string units_string = internal_units.format();
      verbPrintf(2, m_com,
                 "PISM WARNING: Variable '%s' ('%s') does not have the units attribute.\n"
                 "              Assuming that it is in '%s'.\n",
                 name.c_str(), long_name.c_str(),
                 units_string.c_str());
      input_units = internal_units;
    }

    convert_doubles(&data[0], data.size(), input_units, internal_units);

  } catch (RuntimeError &e) {
    e.add_context("reading time-series variable '%s' from '%s'", name.c_str(), inq_filename().c_str());
    throw;
  }
}

void PIO::write_timeseries(const NCTimeseries &metadata, size_t t_start,
                           double data, IO_Type nctype) const {
  vector<double> vector_data(1, data);

  // this call will handle errors
  write_timeseries(metadata, t_start, vector_data, nctype);
}

/** @brief Write a time-series `data` to a file.
 *
 * Always use glaciological units when saving time-series.
 */
void PIO::write_timeseries(const NCTimeseries &metadata, size_t t_start,
                           vector<double> &data,
                           IO_Type nctype) const {
  string name = metadata.get_name();
  try {
    bool variable_exists = inq_var(name);

    if (variable_exists == false) {
      metadata.define(*this, nctype, true);
    }

    // create a copy of "data":
    vector<double> tmp = data;

    // convert to glaciological units:
    convert_doubles(&tmp[0], tmp.size(),
                    metadata.get_units(),
                    metadata.get_glaciological_units());

    put_1d_var(name,
               static_cast<unsigned int>(t_start),
               static_cast<unsigned int>(tmp.size()), tmp);

  } catch (RuntimeError &e) {
    e.add_context("writing time-series variable '%s' to '%s'", name.c_str(), inq_filename().c_str());
    throw;
  }
}

void PIO::read_time_bounds(const NCTimeBounds &metadata,
                           Time *time,
                           vector<double> &data) const {
  string name = metadata.get_name();

  try {
    Unit internal_units = metadata.get_units();

    // Find the variable:
    bool variable_exists = inq_var(name);

    if (variable_exists == false) {
      throw RuntimeError("variable " + name + " is missing");
    }

    vector<string> dims = inq_vardims(name);

    if (dims.size() != 2) {
      throw RuntimeError("variable " + name + " has to has two dimensions");
    }

    string
      &dimension_name = dims[0],
      &bounds_name    = dims[1];

    // Check that we have 2 vertices (interval end-points) per time record.
    unsigned int length = inq_dimlen(bounds_name);
    if (length != 2) {
      throw RuntimeError("time-bounds variable " + name + " has to have exactly 2 bounds per time record");
    }

    // Get the number of time records.
    length = inq_dimlen(dimension_name);
    if (length <= 0) {
      throw RuntimeError("dimension " + dimension_name + " has length zero");
    }

    data.resize(2*length);                // memory allocation happens here

    vector<unsigned int> start(2), count(2);
    start[0] = 0;
    start[1] = 0;
    count[0] = length;
    count[1] = 2;

    get_vara_double(name, start, count, &data[0]);

    // Find the corresponding 'time' variable. (We get units from the 'time'
    // variable, because according to CF-1.5 section 7.1 a "boundary variable"
    // may not have metadata set.)
    variable_exists = inq_var(dimension_name);

    if (variable_exists == false) {
      throw RuntimeError("time coordinate variable " + dimension_name + " is missing");
    }

    bool input_has_units = false;
    Unit input_units(internal_units.get_system(), "1");

    string input_units_string = get_att_text(dimension_name, "units");
    input_units_string = time->CF_units_to_PISM_units(input_units_string);

    if (input_units_string.empty() == true) {
      input_has_units = false;
    } else {
      input_units = Unit(internal_units.get_system(), input_units_string);
      input_has_units = true;
    }

    if (metadata.has_attribute("units") && input_has_units == false) {
      string units_string = internal_units.format();
      verbPrintf(2, m_com,
                 "PISM WARNING: Variable '%s' does not have the units attribute.\n"
                 "              Assuming that it is in '%s'.\n",
                 dimension_name.c_str(),
                 units_string.c_str());
      input_units = internal_units;
    }

    convert_doubles(&data[0], data.size(), input_units, internal_units);

    // FIXME: check that time intervals described by the time bounds
    // variable are contiguous (without gaps) and stop if they are not.
  } catch (RuntimeError &e) {
    e.add_context("reading time bounds variable '%s' from '%s'", name.c_str(), inq_filename().c_str());
    throw;
  }
}

void PIO::write_time_bounds(const NCTimeBounds &metadata,
                            size_t t_start,
                            vector<double> &data, IO_Type nctype) const {
  string name = metadata.get_name();
  try {
    bool variable_exists = inq_var(name);
    if (variable_exists == false) {
      metadata.define(*this, nctype, true);
    }

    // make a copy of "data"
    vector<double> tmp = data;

    // convert to glaciological units:
    convert_doubles(&tmp[0], tmp.size(),
                    metadata.get_units(), metadata.get_glaciological_units());

    vector<unsigned int> start(2), count(2);
    start[0] = static_cast<unsigned int>(t_start);
    start[1] = 0;
    count[0] = static_cast<unsigned int>(tmp.size()) / 2;
    count[1] = 2;

    put_vara_double(name, start, count, &tmp[0]);

  } catch (RuntimeError &e) {
    e.add_context("writing time-bounds variable '%s' to '%s'", name.c_str(), inq_filename().c_str());
    throw;
  }
}

grid_info::grid_info() {

  t_len = 0;
  time  = 0;

  x_len = 0;
  x0    = 0;
  Lx    = 0;

  y_len = 0;
  y0    = 0;
  Ly    = 0;

  z_len = 0;
  z_min = 0;
  z_max = 0;
}

} // end of namespace pism
