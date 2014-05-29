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
#include <assert.h>

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

static NCFile* create_backend(MPI_Comm com, std::string mode) {
  if (mode == "netcdf3") {
    return new NC3File(com);
  } else if (mode.find("quilt") == 0) {
    size_t n = mode.find(":");
    int compression_level = 0;

    if (n != std::string::npos) {
      mode.replace(0, 6, "");   // 6 is the length of "quilt:"
      char *endptr;
      compression_level = strtol(mode.c_str(), &endptr, 10);
      if ((*endptr != '\0') || (compression_level < 0) || (compression_level > 9)) {
        PetscPrintf(com, "PISM WARNING: invalid compression level %s. Output compression is disabled.\n",
                    mode.c_str());
        compression_level = 0;
      }
    }

    return new NC4_Quilt(com, compression_level);
  }
#if (PISM_USE_PARALLEL_NETCDF4==1)
  else if (mode == "netcdf4_parallel") {
    return new NC4_Par(com);
  }
#endif
#if (PISM_USE_PNETCDF==1)
  else if (mode == "pnetcdf") {
    return new PNCFile(com);
  }
#endif
#if (PISM_USE_HDF5==1)
  else if (mode == "hdf5") {
    return new NC4_HDF5(com);
  }
#endif
  else {
    return NULL;
  }
}

//! \brief The code shared by different PIO constructors.
void PIO::constructor(MPI_Comm c, const std::string &mode) {
  m_com = c;
  shallow_copy = false;
  m_mode = mode;

  nc = create_backend(m_com, mode);

  if (mode != "guess_mode" && nc == NULL) {
    PetscPrintf(m_com,
                "PISM ERROR: output format '%s' is not supported.\n"
                "Please recompile PISM with the appropriate I/O library.\n",
                mode.c_str());
    PISMEnd();
  }
}

PIO::PIO(MPI_Comm c, const std::string &mode, const UnitSystem &units_system)
  : m_unit_system(units_system) {
  constructor(c, mode);
}

PIO::PIO(IceGrid &grid, const std::string &mode)
  : m_unit_system(grid.get_unit_system()) {
  constructor(grid.com, mode);
  if (nc != NULL)
    set_local_extent(grid.xs, grid.xm, grid.ys, grid.ym);
}

PIO::PIO(const PIO &other)
  : m_unit_system(other.m_unit_system) {
  m_com = other.m_com;
  nc = other.nc;
  m_mode = other.m_mode;

  shallow_copy = true;
}

PIO::~PIO() {
  if (shallow_copy == false)
    delete nc;
}

// Chooses the best I/O backend for reading from 'filename'.
PetscErrorCode PIO::detect_mode(const std::string &filename) {
  int stat;

  // Bail if someone made a decision already.
  if (nc != NULL)
    return 1;

  NC3File nc3(m_com);

  stat = nc3.open(filename, PISM_READONLY);
  if (stat != 0) {
    PetscPrintf(m_com, "PISM ERROR: Can't open '%s'. Exiting...\n", filename.c_str());
    PISMEnd();
  }

  std::string format = nc3.get_format();

  stat = nc3.close(); CHKERRQ(stat);

  std::vector<std::string> modes;
  if (format == "netcdf4") {
    modes.push_back("netcdf4_parallel");
    modes.push_back("netcdf3");
  } else {
    modes.push_back("pnetcdf");
    modes.push_back("netcdf3");
  }

  for (unsigned int j = 0; j < modes.size(); ++j) {
    nc = create_backend(m_com, modes[j]);

    if (nc != NULL) {
      m_mode = modes[j];
      stat = verbPrintf(3, m_com,
                        "  - Using the %s backend to read from %s...\n",
                        modes[j].c_str(), filename.c_str()); CHKERRQ(stat);
      break;
    }
  }

  if (nc == NULL) {
    PetscPrintf(m_com, "PISM ERROR: Unable to allocate an I/O backend. This should never happen!\n");
    PISMEnd();
  }

  nc->set_local_extent(m_xs, m_xm, m_ys, m_ym);

  return 0;
}

PetscErrorCode PIO::check_if_exists(const std::string &filename, bool &result) {
  PetscErrorCode ierr;

  int file_exists = 0, rank = 0;
  MPI_Comm_rank(m_com, &rank);

  if (rank == 0) {
    // Check if the file exists:
    if (FILE *f = fopen(filename.c_str(), "r")) {
      file_exists = 1;
      fclose(f);
    } else {
      file_exists = 0;
    }
  }
  ierr = MPI_Bcast(&file_exists, 1, MPI_INT, 0, m_com); CHKERRQ(ierr);

  if (file_exists == 1)
    result = true;
  else
    result = false;

  return 0;
}


PetscErrorCode PIO::open(const std::string &filename, IO_Mode mode) {
  PetscErrorCode stat;

  // opening for reading

  if (mode == PISM_READONLY) {
    if (nc == NULL && m_mode == "guess_mode") {
      stat = detect_mode(filename); CHKERRQ(stat);
    }

    assert(nc != NULL);

    stat = nc->open(filename, mode);
    if (stat != 0) {
      PetscPrintf(m_com, "PISM ERROR: Can't open '%s'. Exiting...\n", filename.c_str());
      PISMEnd();
    }
    return 0;
  }

  // opening for writing

  if (mode == PISM_READWRITE_CLOBBER || mode == PISM_READWRITE_MOVE) {

    assert(nc != NULL);

    if (mode == PISM_READWRITE_MOVE) {
      stat = nc->move_if_exists(filename); CHKERRQ(stat);
    } else {
      stat = nc->remove_if_exists(filename); CHKERRQ(stat);
    }

    stat = nc->create(filename);
    if (stat != 0) {
      PetscPrintf(m_com, "PISM ERROR: Can't create '%s'. Exiting...\n", filename.c_str());
      PISMEnd();
    }

    int old_fill;
    stat = nc->set_fill(PISM_NOFILL, old_fill); CHKERRQ(stat);
  } else {                      // mode == PISM_READWRITE
    if (nc == NULL && m_mode == "guess_mode") {
      stat = detect_mode(filename); CHKERRQ(stat);
    }

    assert(nc != NULL);

    stat = nc->open(filename, mode);

    if (stat != 0) {
      int rank = 0;
      MPI_Comm_rank(m_com, &rank);
      PetscPrintf(m_com, "PISM ERROR: Can't open '%s' (rank = %d). Exiting...\n", filename.c_str(), rank);
      PISMEnd();
    }

    int old_fill;
    stat = nc->set_fill(PISM_NOFILL, old_fill); CHKERRQ(stat);

  }

  return 0;
}


PetscErrorCode PIO::close() {
  PetscErrorCode ierr = nc->close(); CHKERRQ(ierr);
  return 0;
}

PetscErrorCode PIO::redef() const {
  PetscErrorCode ierr = nc->redef(); CHKERRQ(ierr);
  return 0;
}


PetscErrorCode PIO::enddef() const {
  PetscErrorCode ierr = nc->enddef(); CHKERRQ(ierr);
  return 0;
}

std::string PIO::inq_filename() const {
  return nc->get_filename();
}


//! \brief Get the number of records. Uses the length of an unlimited dimension.
PetscErrorCode PIO::inq_nrecords(unsigned int &result) const {
  PetscErrorCode ierr;
  std::string dim;

  ierr = nc->inq_unlimdim(dim); CHKERRQ(ierr);

  if (dim.empty()) {
    result = 1;
  } else {
    ierr = nc->inq_dimlen(dim, result); CHKERRQ(ierr);
  }

  return 0;
}

//! \brief Get the number of records of a certain variable. Uses the length of
//! an associated "time" dimension.
PetscErrorCode PIO::inq_nrecords(const std::string &name, const std::string &std_name, unsigned int &result) const {
  PetscErrorCode ierr;

  bool exists = false, found_by_standard_name = false;
  std::string name_found;

  ierr = this->inq_var(name, std_name, exists, name_found, found_by_standard_name); CHKERRQ(ierr);

  if (exists == false) {
    result = 0;
    return 0;
  }

  std::vector<std::string> dims;
  ierr = nc->inq_vardimid(name_found, dims); CHKERRQ(ierr);

  for (unsigned int j = 0; j < dims.size(); ++j) {
    AxisType dimtype;

    ierr = this->inq_dimtype(dims[j], dimtype); CHKERRQ(ierr);

    if (dimtype == T_AXIS) {
      ierr = nc->inq_dimlen(dims[j], result); CHKERRQ(ierr);
      return 0;
    }
  }

  result = 1;

  return 0;
}


//! \brief Find a variable using its standard name and/or short name.
/*!
 * Sets "result" to the short name found.
 */
PetscErrorCode PIO::inq_var(const std::string &short_name, const std::string &std_name, bool &exists,
                            std::string &result, bool &found_by_standard_name) const {
  PetscErrorCode ierr;

  exists = false;

  if (std_name.empty() == false) {
    int nvars;

    ierr = nc->inq_nvars(nvars); CHKERRQ(ierr);

    for (int j = 0; j < nvars; ++j) {
      std::string name, attribute;
      ierr = nc->inq_varname(j, name); CHKERRQ(ierr);

      ierr = nc->get_att_text(name, "standard_name", attribute); CHKERRQ(ierr);

      if (attribute.empty())
        continue;

      if (attribute == std_name) {
        if (exists == false) {
          exists = true;
          found_by_standard_name = true;
          result = name;
        } else {
          ierr = PetscPrintf(m_com,
                             "PISM ERROR: Inconsistency in the input file %s:\n  "
                             "Variables '%s' and '%s' have the same standard_name ('%s').\n",
                             inq_filename().c_str(), result.c_str(), name.c_str(), attribute.c_str());
          CHKERRQ(ierr);
          PISMEnd();
        }
      }

    } // end of the for loop
  } // end of if (std_name.empty() == false)

  if (exists == false) {
    ierr = nc->inq_varid(short_name, exists); CHKERRQ(ierr);
    if (exists == true)
      result = short_name;
    else
      result.clear();

    found_by_standard_name = false;
  }

  return 0;
}

//! \brief Checks if a variable exists.
PetscErrorCode PIO::inq_var(const std::string &name, bool &exists) const {

  PetscErrorCode ierr = nc->inq_varid(name, exists); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode PIO::inq_vardims(const std::string &name, std::vector<std::string> &result) const {
  PetscErrorCode ierr = nc->inq_vardimid(name, result); CHKERRQ(ierr);
  return 0;
}


//! \brief Checks if a dimension exists.
PetscErrorCode PIO::inq_dim(const std::string &name, bool &exists) const {

  PetscErrorCode ierr = nc->inq_dimid(name, exists); CHKERRQ(ierr);

  return 0;
}

//! \brief Get the length of a dimension.
/*!
 * Sets result to 0 if a dimension does not exist.
 */
PetscErrorCode PIO::inq_dimlen(const std::string &name, unsigned int &result) const {
  bool exists = false;
  PetscErrorCode ierr;

  ierr = nc->inq_dimid(name, exists); CHKERRQ(ierr);

  if (exists == true) {
    ierr = nc->inq_dimlen(name, result); CHKERRQ(ierr);
  } else {
    result = 0;
  }

  return 0;
}

//! \brief Get the "type" of a dimension.
/*!
 * The "type" is one of X_AXIS, Y_AXIS, Z_AXIS, T_AXIS.
 */
PetscErrorCode PIO::inq_dimtype(const std::string &name, AxisType &result) const {
  PetscErrorCode ierr;
  std::string axis, standard_name, units;
  Unit tmp_units(m_unit_system, "1");
  bool exists;

  ierr = nc->inq_varid(name, exists); CHKERRQ(ierr);

  if (exists == false) {
    PetscPrintf(m_com, "ERROR: coordinate variable '%s' is not present!\n", name.c_str());
    PISMEnd();
  }

  ierr = nc->get_att_text(name, "axis", axis); CHKERRQ(ierr);
  ierr = nc->get_att_text(name, "standard_name", standard_name); CHKERRQ(ierr);
  ierr = nc->get_att_text(name, "units", units); CHKERRQ(ierr);

  // check if it has units compatible with "seconds":

  try {
    tmp_units = Unit(m_unit_system, units);
  }
  catch (RuntimeError &e) {
    std::string message = "processing variable " + name + " in " + this->inq_filename();
    e.add_context(message);
    throw;
  }

  Unit seconds(m_unit_system, "seconds");
  if (UnitConverter::are_convertible(tmp_units, seconds)) {
    result = T_AXIS;
    return 0;
  }

  // check the standard_name attribute:
  if (standard_name == "time") {
    result = T_AXIS;
    return 0;
  } else if (standard_name == "projection_x_coordinate") {
    result = X_AXIS;
    return 0;
  } else if (standard_name == "projection_y_coordinate") {
    result = Y_AXIS;
    return 0;
  }

  // check the axis attribute:
  if (axis == "T" || axis == "t") {
    result = T_AXIS;
    return 0;
  } else if (axis == "X" || axis == "x") {
    result = X_AXIS;
    return 0;
  } else if (axis == "Y" || axis == "y") {
    result = Y_AXIS;
    return 0;
  } else if (axis == "Z" || axis == "z") {
    result = Z_AXIS;
    return 0;
  }

  // check the variable name:
  if (name == "x" || name == "X" ||
      name.find("x") == 0 || name.find("X") == 0) {
    result = X_AXIS;
    return 0;
  } else if (name == "y" || name == "Y" ||
             name.find("y") == 0 || name.find("Y") == 0) {
    result = Y_AXIS;
    return 0;
  } else if (name == "z" || name == "Z" ||
             name.find("z") == 0 || name.find("Z") == 0) {
    result = Z_AXIS;
    return 0;
  } else if (name == "t" || name == "T" || name == "time" ||
             name.find("t") == 0 || name.find("T") == 0) {
    result = T_AXIS;
    return 0;
  }

  // we have no clue:
  result = UNKNOWN_AXIS;
  return 0;
}

PetscErrorCode PIO::inq_dim_limits(const std::string &name, double *min, double *max) const {
  PetscErrorCode ierr;
  std::vector<double> data;

  ierr = this->get_dim(name, data); CHKERRQ(ierr);

  double my_min = data[0],
    my_max = data[0];
  for (unsigned int j = 0; j < data.size(); ++j) {
    my_min = PetscMin(data[j], my_min);
    my_max = PetscMax(data[j], my_max);
  }

  if (min != NULL)
    *min = my_min;

  if (max != NULL)
    *max = my_max;

  return 0;
}


//! \brief Sets grid parameters using data read from the file.
PetscErrorCode PIO::inq_grid(const std::string &var_name, IceGrid *grid, Periodicity periodicity) const {
  PetscErrorCode ierr;

  if (grid == NULL)
    SETERRQ(m_com, 1, "grid == NULL");

  grid_info input;

  // The following call may fail because var_name does not exist. (And this is fatal!)
  ierr = this->inq_grid_info(var_name, periodicity, input); CHKERRQ(ierr);

  // if we have no vertical grid information, create a fake 2-level vertical grid.
  if (input.z.size() < 2) {
    double Lz = grid->config.get("grid_Lz");
    ierr = verbPrintf(3, m_com,
                      "WARNING: Can't determine vertical grid information using '%s' in %s'\n"
                      "         Using 2 levels and Lz of %3.3fm\n",
                      var_name.c_str(), this->inq_filename().c_str(), Lz); CHKERRQ(ierr);

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

  return 0;
}

/*! Do not use this method to get units of time and time_bounds
  variables: in these two cases we need to handle the reference date
  correctly.
 */
PetscErrorCode PIO::inq_units(const std::string &name, bool &has_units, Unit &units) const {
  PetscErrorCode ierr;
  std::string units_string;

  // Get the string:
  ierr = nc->get_att_text(name, "units", units_string); CHKERRQ(ierr);

  // If a variables does not have the units attribute, set the flag and return:
  if (units_string.empty()) {
    has_units = false;
    units = Unit(units.get_system(), "1");
    return 0;
  }

  // strip trailing spaces
  while (ends_with(units_string, " "))
    units_string.resize(units_string.size() - 1);
  
  try {
    units = Unit(units.get_system(), units_string);
  }
  catch (RuntimeError &e) {
    std::string message = "processing variable " + name + " in " + this->inq_filename();
    e.add_context(message);
    throw;
  }

  has_units = true;

  return 0;
}


PetscErrorCode PIO::inq_grid_info(const std::string &name, Periodicity p,
                                  grid_info &g) const {
  PetscErrorCode ierr;

  std::vector<std::string> dims;
  bool exists, found_by_standard_name;
  std::string name_found;

  // try "name" as the standard_name first, then as the short name:
  ierr = this->inq_var(name, name, exists, name_found, found_by_standard_name); CHKERRQ(ierr);

  if (exists == false)
    SETERRQ2(m_com, 1, "Could not find variable %s in %s", name.c_str(),
             this->inq_filename().c_str());

  ierr = nc->inq_vardimid(name_found, dims); CHKERRQ(ierr);

  // use "global" dimensions (as opposed to dimensions of a patch)
  if (m_mode == "quilt") {
    for (unsigned int i = 0; i < dims.size(); ++i) {
      if (dims[i] == "x_patch")
        dims[i] = "x";

      if (dims[i] == "y_patch")
        dims[i] = "y";
    }
  }

  for (unsigned int i = 0; i < dims.size(); ++i) {
    std::string dimname = dims[i];

    AxisType dimtype = UNKNOWN_AXIS;
    ierr = this->inq_dimtype(dimname, dimtype); CHKERRQ(ierr);

    switch (dimtype) {
    case X_AXIS:
      {
        ierr = nc->inq_dimlen(dimname, g.x_len); CHKERRQ(ierr);
        double x_min = 0.0, x_max = 0.0;
        ierr = this->inq_dim_limits(dimname, &x_min, &x_max); CHKERRQ(ierr);
        ierr = this->get_dim(dimname, g.x); CHKERRQ(ierr);
        g.x0 = 0.5 * (x_min + x_max);
        g.Lx = 0.5 * (x_max - x_min);
        if (p & X_PERIODIC) {
          const double dx = g.x[1] - g.x[0];
          g.Lx += 0.5 * dx;
        }
        break;
      }
    case Y_AXIS:
      {
        ierr = nc->inq_dimlen(dimname, g.y_len); CHKERRQ(ierr);
        double y_min = 0.0, y_max = 0.0;
        ierr = this->inq_dim_limits(dimname, &y_min, &y_max); CHKERRQ(ierr);
        ierr = this->get_dim(dimname, g.y); CHKERRQ(ierr);
        g.y0 = 0.5 * (y_min + y_max);
        g.Ly = 0.5 * (y_max - y_min);
        if (p & Y_PERIODIC) {
          const double dy = g.y[1] - g.y[0];
          g.Ly += 0.5 * dy;
        }
        break;
      }
    case Z_AXIS:
      {
        ierr = nc->inq_dimlen(dimname, g.z_len); CHKERRQ(ierr);
        ierr = this->inq_dim_limits(dimname, &g.z_min, &g.z_max); CHKERRQ(ierr);
        ierr = this->get_dim(dimname, g.z); CHKERRQ(ierr);
        break;
      }
    case T_AXIS:
      {
        ierr = nc->inq_dimlen(dimname, g.t_len); CHKERRQ(ierr);
        ierr = this->inq_dim_limits(dimname, NULL, &g.time); CHKERRQ(ierr);
        break;
      }
    default:
      {
        PetscPrintf(m_com, "ERROR: Can't figure out which direction dimension '%s' corresponds to.",
                    dimname.c_str());
        PISMEnd();
      }
    } // switch
  }   // for loop

  return 0;
}

//! \brief Define a dimension \b and the associated coordinate variable. Set attributes.
PetscErrorCode PIO::def_dim(unsigned long int length,
                            const NCVariable &metadata) const {
  PetscErrorCode ierr;

  std::string name = metadata.get_name();

  ierr = nc->redef(); CHKERRQ(ierr);

  ierr = nc->def_dim(name, length); CHKERRQ(ierr);

  std::vector<std::string> dims(1, name);
  ierr = nc->def_var(name, PISM_DOUBLE, dims); CHKERRQ(ierr);

  ierr = this->write_attributes(metadata, PISM_DOUBLE, false); CHKERRQ(ierr);

  return 0;
}

//! \brief Define a variable.
PetscErrorCode PIO::def_var(const std::string &name, IO_Type nctype, const std::vector<std::string> &dims) const {
  PetscErrorCode ierr;

  ierr = nc->def_var(name, nctype, dims); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode PIO::get_1d_var(const std::string &name, unsigned int s, unsigned int c,
                               std::vector<double> &result) const {
  PetscErrorCode ierr;
  std::vector<unsigned int> start(1), count(1);

  result.resize(c);

  start[0] = s;
  count[0] = c;

  ierr = nc->enddef(); CHKERRQ(ierr);

  ierr = nc->get_vara_double(name, start, count, &result[0]); CHKERRQ(ierr);

  return 0;
}


PetscErrorCode PIO::put_1d_var(const std::string &name, unsigned int s, unsigned int c,
                               const std::vector<double> &data) const {
  PetscErrorCode ierr;
  std::vector<unsigned int> start(1), count(1);

  start[0] = s;
  count[0] = c;

  ierr = nc->enddef(); CHKERRQ(ierr);

  ierr = nc->put_vara_double(name, start, count,
                             const_cast<double*>(&data[0])); CHKERRQ(ierr);

  return 0;
}


//! \brief Get dimension data (a coordinate variable).
PetscErrorCode PIO::get_dim(const std::string &name, std::vector<double> &data) const {
  PetscErrorCode ierr;

  unsigned int dim_length = 0;
  ierr = nc->inq_dimlen(name, dim_length); CHKERRQ(ierr);

  ierr = nc->enddef(); CHKERRQ(ierr);

  ierr = this->get_1d_var(name, 0, dim_length, data); CHKERRQ(ierr);

  return 0;
}

//! \brief Write dimension data (a coordinate variable).
PetscErrorCode PIO::put_dim(const std::string &name, const std::vector<double> &data) const {
  PetscErrorCode ierr = this->put_1d_var(name, 0,
                                         (unsigned int)data.size(), data); CHKERRQ(ierr);
  return 0;
}

PetscErrorCode PIO::def_time(const std::string &name, const std::string &calendar, const std::string &units) const {
  PetscErrorCode ierr;
  std::map<std::string,std::string> attrs;

  bool time_exists;
  ierr = nc->inq_varid(name, time_exists); CHKERRQ(ierr);
  if (time_exists)
    return 0;

  // t
  NCVariable t(name, m_unit_system);
  t.set_string("long_name", "time");
  t.set_string("calendar", calendar);
  ierr = t.set_units(units); CHKERRQ(ierr);
  t.set_string("axis", "T");

  ierr = this->def_dim(PISM_UNLIMITED, t); CHKERRQ(ierr);

  return 0;
}


//! \brief Append to the time dimension.
PetscErrorCode PIO::append_time(const std::string &name, double value) const {
  PetscErrorCode ierr;

  std::vector<unsigned int> start(1), count(1);
  unsigned int dim_length = 0;

  ierr = nc->inq_dimlen(name, dim_length); CHKERRQ(ierr);

  start[0] = dim_length;
  count[0] = 1;

  ierr = nc->enddef(); CHKERRQ(ierr);

  ierr = nc->put_vara_double(name, start, count, &value); CHKERRQ(ierr);

  return 0;
}

//! \brief Append to the history global attribute.
/*!
 * Use put_att_text("PISM_GLOBAL", "history", ...) to overwrite "history".
 */
PetscErrorCode PIO::append_history(const std::string &history) const {
  PetscErrorCode ierr;
  std::string old_history;

  ierr = nc->redef(); CHKERRQ(ierr);

  ierr = nc->get_att_text("PISM_GLOBAL", "history", old_history); CHKERRQ(ierr);
  ierr = nc->put_att_text("PISM_GLOBAL", "history", history + old_history); CHKERRQ(ierr);

  return 0;
}

//! \brief Write a multiple-valued double attribute.
PetscErrorCode PIO::put_att_double(const std::string &var_name, const std::string &att_name, IO_Type nctype,
                                   const std::vector<double> &values) const {
  PetscErrorCode ierr;

  ierr = nc->redef(); CHKERRQ(ierr);

  ierr = nc->put_att_double(var_name, att_name, nctype, values); CHKERRQ(ierr);

  return 0;
}

//! \brief Write a single-valued double attribute.
PetscErrorCode PIO::put_att_double(const std::string &var_name, const std::string &att_name, IO_Type nctype,
                                   double value) const {
  PetscErrorCode ierr;
  std::vector<double> tmp; tmp.push_back(value);

  ierr = nc->redef(); CHKERRQ(ierr);

  ierr = nc->put_att_double(var_name, att_name, nctype, tmp); CHKERRQ(ierr);

  return 0;
}

//! \brief Write a text attribute.
PetscErrorCode PIO::put_att_text(const std::string &var_name, const std::string &att_name, const std::string &value) const {
  PetscErrorCode ierr;

  ierr = nc->redef(); CHKERRQ(ierr);

  std::string tmp = value + "\0";    // ensure that the string is null-terminated

  ierr = nc->put_att_text(var_name, att_name, tmp); CHKERRQ(ierr);

  return 0;
}

//! \brief Get a double attribute.
PetscErrorCode PIO::get_att_double(const std::string &var_name, const std::string &att_name,
                                   std::vector<double> &result) const {

  PetscErrorCode ierr;
  IO_Type att_type;
  // virtual int inq_atttype(string variable_name, string att_name, IO_Type &result) const = 0;

  ierr = nc->inq_atttype(var_name, att_name, att_type); CHKERRQ(ierr);

  // Give an understandable error message if a string attribute was found when
  // a number (or a list of numbers) was expected. (We've seen datasets with
  // "valid_min" stored as a string...)
  if (att_type == PISM_CHAR) {
    std::string tmp;
    ierr = nc->get_att_text(var_name, att_name, tmp); CHKERRQ(ierr);

    PetscPrintf(m_com,
                "PISM ERROR: attribute %s:%s in %s is a string (\"%s\");"
                " expected a number (or a list of numbers).\n",
                var_name.c_str(), att_name.c_str(), nc->get_filename().c_str(), tmp.c_str());
    PISMEnd();
  } else {
    // In this case att_type might be PISM_NAT (if an attribute does not
    // exist), but get_att_double can handle that.
    ierr = nc->get_att_double(var_name, att_name, result); CHKERRQ(ierr);
  }

  return 0;
}

//! \brief Get a text attribute.
PetscErrorCode PIO::get_att_text(const std::string &var_name, const std::string &att_name, std::string &result) const {

  PetscErrorCode ierr = nc->get_att_text(var_name, att_name, result); CHKERRQ(ierr);

  return 0;
}

//! \brief Read a PETSc Vec using the grid "grid".
/*!
 * Assumes that double corresponds to C++ double.
 *
 * Vec g has to be "global" (i.e. without ghosts).
 */
PetscErrorCode PIO::get_vec(IceGrid *grid, const std::string &var_name,
                            unsigned int z_count, unsigned int t_start, Vec g) const {
  PetscErrorCode ierr;

  std::vector<unsigned int> start, count, imap;
  ierr = compute_start_and_count(var_name,
                                 t_start,
                                 grid->xs, grid->xm,
                                 grid->ys, grid->ym,
                                 0, z_count,
                                 start, count, imap); CHKERRQ(ierr);

  ierr = nc->enddef(); CHKERRQ(ierr);

  double *a_petsc;
  ierr = VecGetArray(g, &a_petsc); CHKERRQ(ierr);

  // We always use "mapped" I/O here, because we don't know where the input
  // file came from.
  ierr = nc->get_varm_double(var_name, start, count, imap, (double*)a_petsc); CHKERRQ(ierr);

  ierr = VecRestoreArray(g, &a_petsc); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode PIO::inq_nattrs(const std::string &var_name, int &result) const {
  PetscErrorCode ierr = nc->inq_varnatts(var_name, result); CHKERRQ(ierr);
  return 0;
}


PetscErrorCode PIO::inq_attname(const std::string &var_name, unsigned int n, std::string &result) const {
  PetscErrorCode ierr = nc->inq_attname(var_name, n, result); CHKERRQ(ierr);
  return 0;
}


PetscErrorCode PIO::inq_atttype(const std::string &var_name, const std::string &att_name, IO_Type &result) const {
  PetscErrorCode ierr = nc->inq_atttype(var_name, att_name, result); CHKERRQ(ierr);
  return 0;
}



//! \brief Write a PETSc Vec using the grid "grid".
/*!
 * Assumes that double corresponds to C++ double.
 *
 * Vec g has to be "global" (i.e. without ghosts).
 *
 * This method always writes to the last record in the file.
 */
PetscErrorCode PIO::put_vec(IceGrid *grid, const std::string &var_name, unsigned int z_count, Vec g) const {
  PetscErrorCode ierr;

  unsigned int t;
  ierr = nc->inq_dimlen(grid->config.get_string("time_dimension_name"), t); CHKERRQ(ierr);

#if (PISM_DEBUG==1)
  if (t < 1)
    fprintf(stderr, "time dimension length (%d) is less than 1!\n", t);
#endif

  std::vector<unsigned int> start, count, imap;
  ierr = compute_start_and_count(var_name,
                                 t - 1,
                                 grid->xs, grid->xm,
                                 grid->ys, grid->ym,
                                 0, z_count,
                                 start, count, imap); CHKERRQ(ierr);

  ierr = nc->enddef(); CHKERRQ(ierr);

  double *a_petsc;
  ierr = VecGetArray(g, &a_petsc); CHKERRQ(ierr);

  if (grid->config.get_string("output_variable_order") == "xyz") {
    // Use the faster and safer (avoids a NetCDF bug) call if the aray storage
    // orders in the memory and in NetCDF files are the same.
    ierr = nc->put_vara_double(var_name, start, count, (double*)a_petsc); CHKERRQ(ierr);
  } else {
    // Revert to "mapped" I/O otherwise.
    ierr = nc->put_varm_double(var_name, start, count, imap, (double*)a_petsc); CHKERRQ(ierr);
  }

  ierr = VecRestoreArray(g, &a_petsc); CHKERRQ(ierr);

  return 0;
}

//! \brief Get the interpolation context (grid information) for an input file.
/*!
 * Sets lic to NULL if the variable was not found.
 *
 * @note The *caller* is in charge of destroying lic
 */
PetscErrorCode PIO::get_interp_context(const std::string &name,
                                       const IceGrid &grid,
                                       const std::vector<double> &zlevels,
                                       LocalInterpCtx* &lic) const {
  PetscErrorCode ierr;
  bool exists = false;

  ierr = this->inq_var(name, exists); CHKERRQ(ierr);

  if (exists == false) {
    lic = NULL;
    SETERRQ1(grid.com, 1, "Variable %s was not found", name.c_str());
  } else {
    grid_info gi;

    ierr = this->inq_grid_info(name, grid.periodicity, gi); CHKERRQ(ierr);

    lic = new LocalInterpCtx(gi, grid, zlevels.front(), zlevels.back());
  }

  return 0;
}

//! \brief Read a PETSc Vec from a file, using bilinear (or trilinear)
//! interpolation to put it on the grid defined by "grid" and zlevels_out.
PetscErrorCode PIO::regrid_vec(IceGrid *grid, const std::string &var_name,
                               const std::vector<double> &zlevels_out,
                               unsigned int t_start, Vec g) const {
  PetscErrorCode ierr;
  const int X = 1, Y = 2, Z = 3; // indices, just for clarity
  std::vector<unsigned int> start, count, imap;

  LocalInterpCtx *lic = NULL;

  ierr = get_interp_context(var_name, *grid, zlevels_out, lic); CHKERRQ(ierr);

  assert(lic != NULL);

  ierr = compute_start_and_count(var_name, t_start,
                                 lic->start[X], lic->count[X],
                                 lic->start[Y], lic->count[Y],
                                 lic->start[Z], lic->count[Z],
                                 start, count, imap); CHKERRQ(ierr);

  ierr = nc->enddef(); CHKERRQ(ierr);

  // We always use "mapped" I/O here, because we don't know where the input
  // file came from.
  ierr = nc->get_varm_double(var_name, start, count, imap, lic->a); CHKERRQ(ierr);

  ierr = regrid(grid, zlevels_out, lic, g); CHKERRQ(ierr);

  delete lic;

  return 0;
}

/** Regrid `var_name` from a file, replacing missing values with `default_value`.
 *
 * @param grid computational grid; used to initialize interpolation
 * @param var_name variable to regrid
 * @param zlevels_out vertical levels of the resulting grid
 * @param t_start time index of the record to regrid
 * @param default_value default value to replace `_FillValue` with
 * @param[out] g resulting interpolated field
 *
 * @return 0 on success
 */
PetscErrorCode PIO::regrid_vec_fill_missing(IceGrid *grid, const std::string &var_name,
                                            const std::vector<double> &zlevels_out,
                                            unsigned int t_start,
                                            double default_value,
                                            Vec g) const {
  PetscErrorCode ierr;
  const int X = 1, Y = 2, Z = 3; // indices, just for clarity
  std::vector<unsigned int> start, count, imap;

  LocalInterpCtx *lic = NULL;

  ierr = get_interp_context(var_name, *grid, zlevels_out, lic); CHKERRQ(ierr);

  assert(lic != NULL);

  ierr = compute_start_and_count(var_name, t_start,
                                 lic->start[X], lic->count[X],
                                 lic->start[Y], lic->count[Y],
                                 lic->start[Z], lic->count[Z],
                                 start, count, imap); CHKERRQ(ierr);

  ierr = nc->enddef(); CHKERRQ(ierr);

  // We always use "mapped" I/O here, because we don't know where the input
  // file came from.
  ierr = nc->get_varm_double(var_name, start, count, imap, lic->a); CHKERRQ(ierr);

  // Replace missing values if the _FillValue attribute is present,
  // and if we have missing values to replace.
  {
    std::vector<double> attribute;
    ierr = nc->get_att_double(var_name, "_FillValue", attribute); CHKERRQ(ierr);
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

  ierr = regrid(grid, zlevels_out, lic, g); CHKERRQ(ierr);

  delete lic;

  return 0;
}

int PIO::k_below(double z, const std::vector<double> &zlevels) const {
  double z_min = zlevels.front(), z_max = zlevels.back();
  int mcurr = 0;

  if (z < z_min - 1.0e-6 || z > z_max + 1.0e-6) {
    PetscPrintf(m_com,
                "PIO::k_below(): z = %5.4f is outside the allowed range.\n", z);
    PISMEnd();
  }

  while (zlevels[mcurr+1] < z)
    mcurr++;

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
 * - output array (Vec g)
 *
 * We should be able to switch to using an external interpolation library
 * fairly easily...
 */
PetscErrorCode PIO::regrid(IceGrid *grid, const std::vector<double> &zlevels_out,
                           LocalInterpCtx *lic, Vec g) const {
  const int Y = 2, Z = 3; // indices, just for clarity
  PetscErrorCode ierr;

  std::vector<double> &zlevels_in = lic->zlevels;
  unsigned int nlevels = zlevels_out.size();
  double *input_array = lic->a;

  // array sizes for mapping from logical to "flat" indices
  int y_count = lic->count[Y],
    z_count = lic->count[Z];

  // We'll work with the raw storage here so that the array we are filling is
  // indexed the same way as the buffer we are pulling from (input_array)
  double *output_array;
  ierr = VecGetArray(g, &output_array); CHKERRQ(ierr);

  for (int i = grid->xs; i < grid->xs + grid->xm; i++) {
    for (int j = grid->ys; j < grid->ys + grid->ym; j++) {

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
  }

  ierr = VecRestoreArray(g, &output_array); CHKERRQ(ierr);

  return 0;
}


PetscErrorCode PIO::compute_start_and_count(const std::string &short_name, unsigned int t_start,
                                            unsigned int x_start, unsigned int x_count,
                                            unsigned int y_start, unsigned int y_count,
                                            unsigned int z_start, unsigned int z_count,
                                            std::vector<unsigned int> &start,
                                            std::vector<unsigned int> &count,
                                            std::vector<unsigned int> &imap) const {
  PetscErrorCode ierr;
  std::vector<std::string> dims;

  ierr = nc->inq_vardimid(short_name, dims); CHKERRQ(ierr);
  unsigned int ndims = dims.size();

  // Resize output vectors:
  start.resize(ndims);
  count.resize(ndims);
  imap.resize(ndims);

  // Assemble start, count and imap:
  for (unsigned int j = 0; j < ndims; j++) {
    std::string dimname = dims[j];

    AxisType dimtype;
    ierr = this->inq_dimtype(dimname, dimtype); CHKERRQ(ierr);

    switch (dimtype) {
    case T_AXIS:
      start[j] = t_start;
      count[j] = 1;             // t_count is always 1
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

  return 0;
}

PetscErrorCode PIO::get_vara_double(const std::string &variable_name,
                                    const std::vector<unsigned int> &start,
                                    const std::vector<unsigned int> &count,
                                    double *ip) const {
  PetscErrorCode ierr;

  ierr = nc->enddef(); CHKERRQ(ierr);

  ierr = nc->get_vara_double(variable_name, start, count, ip); CHKERRQ(ierr);

  return 0;
}


PetscErrorCode PIO::put_vara_double(const std::string &variable_name,
                                    const std::vector<unsigned int> &start,
                                    const std::vector<unsigned int> &count,
                                    double *op) const {
  PetscErrorCode ierr;

  ierr = nc->enddef(); CHKERRQ(ierr);

  ierr = nc->put_vara_double(variable_name, start, count, op); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode PIO::get_varm_double(const std::string &variable_name,
                                    const std::vector<unsigned int> &start,
                                    const std::vector<unsigned int> &count,
                                    const std::vector<unsigned int> &imap, double *ip) const {
  PetscErrorCode ierr;

  ierr = nc->enddef(); CHKERRQ(ierr);

  ierr = nc->get_varm_double(variable_name, start, count, imap, ip); CHKERRQ(ierr);

  return 0;
}


PetscErrorCode PIO::put_varm_double(const std::string &variable_name,
                                    const std::vector<unsigned int> &start,
                                    const std::vector<unsigned int> &count,
                                    const std::vector<unsigned int> &imap, double *op) const {
  PetscErrorCode ierr;

  ierr = nc->enddef(); CHKERRQ(ierr);

  ierr = nc->put_varm_double(variable_name, start, count, imap, op); CHKERRQ(ierr);

  return 0;
}

void PIO::set_local_extent(unsigned int xs, unsigned int xm,
                           unsigned int ys, unsigned int ym) {
  nc->set_local_extent(xs, xm, ys, ym);
  m_xs = xs;
  m_xm = xm;
  m_ys = ys;
  m_ym = ym;
}

PetscErrorCode PIO::read_attributes(const std::string &name, NCVariable &variable) const {
  PetscErrorCode ierr;
  bool variable_exists;
  int nattrs;

  ierr = this->inq_var(name, variable_exists); CHKERRQ(ierr);

  if (variable_exists == false) {
    ierr = PetscPrintf(m_com,
                       "PISM ERROR: variable %s was not found in %s.\n"
                       "            Exiting...\n",
                       name.c_str(),
                       this->inq_filename().c_str()); CHKERRQ(ierr);
    PISMEnd();
  }

  variable.clear_all_strings();
  variable.clear_all_doubles();

  ierr = this->inq_nattrs(name, nattrs); CHKERRQ(ierr);

  for (int j = 0; j < nattrs; ++j) {
    std::string attribute_name;
    IO_Type nctype;
    ierr = this->inq_attname(name, j, attribute_name); CHKERRQ(ierr);
    ierr = this->inq_atttype(name, attribute_name, nctype); CHKERRQ(ierr);

    if (nctype == PISM_CHAR) {
      std::string value;
      ierr = this->get_att_text(name, attribute_name, value); CHKERRQ(ierr);

      if (attribute_name == "units") {
        ierr = variable.set_units(value); CHKERRQ(ierr);
      } else {
        variable.set_string(attribute_name, value);
      }
    } else {
      std::vector<double> values;

      ierr = this->get_att_double(name, attribute_name, values); CHKERRQ(ierr);
      variable.set_doubles(attribute_name, values);
    }
  } // end of for (int j = 0; j < nattrs; ++j)

  return 0;
}


//! Write variable attributes to a NetCDF file.
/*!

  \li if write_in_glaciological_units == true, "glaciological_units" are
  written under the name "units" plus the valid range is written in
  glaciological units.

  \li if both valid_min and valid_max are set, then valid_range is written
  instead of the valid_min, valid_max pair.
 */
PetscErrorCode PIO::write_attributes(const NCVariable &var, IO_Type nctype,
                                     bool write_in_glaciological_units) const {
  int ierr;

  std::string var_name = var.get_name();

  ierr = this->redef(); CHKERRQ(ierr);

  // units, valid_min, valid_max and valid_range need special treatment:
  if (var.has_attribute("units")) {
    std::string output_units = var.get_string("units");

    if (write_in_glaciological_units)
      output_units = var.get_string("glaciological_units");

    ierr = this->put_att_text(var_name, "units", output_units); CHKERRQ(ierr);
  }

  std::vector<double> bounds(2);
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
    ierr = this->put_att_double(var_name, "_FillValue", nctype, fill_value); CHKERRQ(ierr);
  }

  if (var.has_attribute("valid_min") && var.has_attribute("valid_max")) {
    ierr = this->put_att_double(var_name, "valid_range", nctype, bounds);
  } else if (var.has_attribute("valid_min")) {
    ierr = this->put_att_double(var_name, "valid_min",   nctype, bounds[0]);
  } else if (var.has_attribute("valid_max")) {
    ierr = this->put_att_double(var_name, "valid_max",   nctype, bounds[1]);
  }

  CHKERRQ(ierr);

  // Write text attributes:
  const NCVariable::StringAttrs &strings = var.get_all_strings();
  NCVariable::StringAttrs::const_iterator i;
  for (i = strings.begin(); i != strings.end(); ++i) {
    std::string
      name  = i->first,
      value = i->second;

    if (name == "units" || name == "glaciological_units" || value.empty())
      continue;

    ierr = this->put_att_text(var_name, name, value); CHKERRQ(ierr);
  }

  // Write double attributes:
  const NCVariable::DoubleAttrs &doubles = var.get_all_doubles();
  NCVariable::DoubleAttrs::const_iterator j;
  for (j = doubles.begin(); j != doubles.end(); ++j) {
    std::string name  = j->first;
    std::vector<double> values = j->second;

    if (name == "valid_min" ||
        name == "valid_max" ||
        name == "valid_range" ||
        name == "_FillValue" ||
        values.empty())
      continue;

    ierr = this->put_att_double(var_name, name, nctype, values); CHKERRQ(ierr);
  }

  return 0;
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
PetscErrorCode PIO::write_global_attributes(const NCVariable &var) const {
  PetscErrorCode ierr;
  std::string old_history;
  NCVariable tmp = var;

  ierr = this->get_att_text("PISM_GLOBAL", "history", old_history); CHKERRQ(ierr);

  tmp.set_name("PISM_GLOBAL");
  tmp.set_string("history", tmp.get_string("history") + old_history);

  ierr = this->write_attributes(tmp, PISM_DOUBLE, false); CHKERRQ(ierr);

  return 0;
}

//! Read the valid range information from a file.
/*! Reads `valid_min`, `valid_max` and `valid_range` attributes; if \c
    valid_range is found, sets the pair `valid_min` and `valid_max` instead.
 */
PetscErrorCode PIO::read_valid_range(const std::string &name, NCVariable &variable) const {
  std::string input_units_string;
  std::vector<double> bounds;
  int ierr;

  // Never reset valid_min/max if any of them was set internally.
  if (variable.has_attribute("valid_min") ||
      variable.has_attribute("valid_max"))
    return 0;

  // Read the units: The following code ignores the units in the input file if
  // a) they are absent :-) b) they are invalid c) they are not compatible with
  // internal units.
  ierr = this->get_att_text(name, "units", input_units_string); CHKERRQ(ierr);

  UnitSystem sys = variable.get_units().get_system();
  Unit input_units = Unit(sys, "1");
  try {
    input_units = Unit(sys, input_units_string);
  }
  catch (...) {
    // ignore the failure above (WHY? -- CK)
    input_units = variable.get_units();
  }

  shared_ptr<UnitConverter> c(new UnitConverter());
  if (UnitConverter::are_convertible(input_units, variable.get_units()) == true) {
    c.reset(new UnitConverter(input_units, variable.get_units()));
  }

  ierr = this->get_att_double(name, "valid_range", bounds); CHKERRQ(ierr);
  if (bounds.size() == 2) {             // valid_range is present
    variable.set_double("valid_min", (*c)(bounds[0]));
    variable.set_double("valid_max", (*c)(bounds[1]));
  } else {                      // valid_range has the wrong length or is missing
    ierr = this->get_att_double(name, "valid_min", bounds); CHKERRQ(ierr);
    if (bounds.size() == 1) {           // valid_min is present
      variable.set_double("valid_min", (*c)(bounds[0]));
    }

    ierr = this->get_att_double(name, "valid_max", bounds); CHKERRQ(ierr);
    if (bounds.size() == 1) {           // valid_max is present
      variable.set_double("valid_max", (*c)(bounds[0]));
    }
  }

  return 0;
}

//! Read a time-series variable from a NetCDF file to a vector of doubles.
PetscErrorCode PIO::read_timeseries(const NCTimeseries &metadata,
                                    Time *time,
                                    std::vector<double> &data) const {
  PetscErrorCode ierr;
  bool variable_exists;

  // Find the variable:
  std::string name_found,
    name           = metadata.get_name(),
    long_name      = metadata.get_string("long_name"),
    standard_name  = metadata.get_string("standard_name"),
    dimension_name = metadata.get_dimension_name();

  bool found_by_standard_name;
  ierr = this->inq_var(name, standard_name,
                       variable_exists, name_found, found_by_standard_name); CHKERRQ(ierr);

  if (!variable_exists) {
    ierr = PetscPrintf(m_com,
                       "PISM ERROR: Can't find '%s' (%s) in '%s'.\n",
                       name.c_str(),
                       standard_name.c_str(), this->inq_filename().c_str());
    CHKERRQ(ierr);
    PISMEnd();
  }

  std::vector<std::string> dims;
  ierr = this->inq_vardims(name_found, dims); CHKERRQ(ierr);

  if (dims.size() != 1) {
    ierr = PetscPrintf(m_com,
                       "PISM ERROR: Variable '%s' in '%s' depends on %d dimensions,\n"
                       "            but a time-series variable can only depend on 1 dimension.\n",
                       name.c_str(), this->inq_filename().c_str(), dims.size()); CHKERRQ(ierr);
    PISMEnd();
  }

  unsigned int length;
  ierr = this->inq_dimlen(dimension_name, length); CHKERRQ(ierr);

  if (length <= 0) {
    ierr = PetscPrintf(m_com,
                       "PISM ERROR: Dimension %s has zero length!\n",
                       dimension_name.c_str()); CHKERRQ(ierr);
    PISMEnd();
  }

  data.resize(length);          // memory allocation happens here

  ierr = this->enddef(); CHKERRQ(ierr);

  ierr = this->get_1d_var(name_found, 0, length, data); CHKERRQ(ierr);

  bool input_has_units;
  std::string input_units_string;
  Unit internal_units = metadata.get_units(),
    input_units(internal_units.get_system(), "1");

  ierr = this->get_att_text(name_found, "units", input_units_string); CHKERRQ(ierr);

  if (input_units_string.empty() == true) {
    input_has_units = false;
  } else {
    input_units_string = time->CF_units_to_PISM_units(input_units_string);

    try {
      input_units = Unit(internal_units.get_system(), input_units_string);
    }
    catch (RuntimeError &e) {
      std::string message = "processing variable " + name + " in " + this->inq_filename();
      e.add_context(message);
      throw;
    }
    input_has_units = true;
  }

  if (metadata.has_attribute("units") == true && input_has_units == false) {
    std::string units_string = internal_units.format();
    ierr = verbPrintf(2, m_com,
                      "PISM WARNING: Variable '%s' ('%s') does not have the units attribute.\n"
                      "              Assuming that it is in '%s'.\n",
                      name.c_str(), long_name.c_str(),
                      units_string.c_str()); CHKERRQ(ierr);
    input_units = internal_units;
  }

  ierr = convert_doubles(&data[0], data.size(),
                         input_units, internal_units); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode PIO::write_timeseries(const NCTimeseries &metadata, size_t t_start,
                                     double data, IO_Type nctype) const {
  PetscErrorCode ierr;
  std::vector<double> vector_data(1, data);

  ierr = this->write_timeseries(metadata, t_start, vector_data, nctype); CHKERRQ(ierr);

  return 0;
}

/** @brief Write a time-series `data` to a file.
 *
 * Always use glaciological units when saving time-series.
 */
PetscErrorCode PIO::write_timeseries(const NCTimeseries &metadata, size_t t_start,
                                     std::vector<double> &data,
                                     IO_Type nctype) const {
  PetscErrorCode ierr;
  bool variable_exists = false;

  ierr = this->inq_var(metadata.get_name(), variable_exists); CHKERRQ(ierr);
  if (variable_exists == false) {
    ierr = metadata.define(*this, nctype, true); CHKERRQ(ierr);
  }

  ierr = this->enddef(); CHKERRQ(ierr);

  // convert to glaciological units:
  ierr = convert_doubles(&data[0], data.size(),
                         metadata.get_units(),
                         metadata.get_glaciological_units()); CHKERRQ(ierr);

  ierr = this->put_1d_var(metadata.get_name(),
                          static_cast<unsigned int>(t_start),
                          static_cast<unsigned int>(data.size()), data); CHKERRQ(ierr);

  // restore internal units:
  ierr = convert_doubles(&data[0], data.size(),
                         metadata.get_glaciological_units(),
                         metadata.get_units()); CHKERRQ(ierr);
  return 0;
}

PetscErrorCode PIO::read_time_bounds(const NCTimeBounds &metadata,
                                     Time *time,
                                     std::vector<double> &data) const {
  PetscErrorCode ierr;
  bool variable_exists = false;

  std::string
    name     = metadata.get_name(),
    filename = this->inq_filename();

  Unit internal_units = metadata.get_units();

  // Find the variable:
  ierr = this->inq_var(name, variable_exists); CHKERRQ(ierr);

  if (variable_exists == false) {
    ierr = PetscPrintf(m_com,
                      "PISM ERROR: Can't find '%s' in '%s'.\n",
                       name.c_str(), filename.c_str()); CHKERRQ(ierr);
    PISMEnd();
  }

  std::vector<std::string> dims;
  ierr = this->inq_vardims(name, dims); CHKERRQ(ierr);

  if (dims.size() != 2) {
    ierr = PetscPrintf(m_com,
                       "PISM ERROR: Variable '%s' in '%s' depends on %d dimensions,\n"
                       "            but a time-bounds variable can only depend on 2 dimension.\n",
                       name.c_str(),
                       filename.c_str(), dims.size()); CHKERRQ(ierr);
    PISMEnd();
  }

  std::string
    &dimension_name = dims[0],
    &bounds_name    = dims[1];

  unsigned int length = 0;

  // Check that we have 2 vertices (interval end-points) per time record.
  ierr = this->inq_dimlen(bounds_name, length); CHKERRQ(ierr);
  if (length != 2) {
    PetscPrintf(m_com,
                "PISM ERROR: A time-bounds variable has to have exactly 2 bounds per time record.\n"
                "            Please check that the dimension corresponding to 'number of vertices' goes\n"
                "            last in the 'ncdump -h %s' output.\n",
                filename.c_str());
    PISMEnd();
  }

  // Get the number of time records.
  ierr = this->inq_dimlen(dimension_name, length); CHKERRQ(ierr);
  if (length <= 0) {
    ierr = PetscPrintf(m_com,
                       "PISM ERROR: Dimension %s has zero length!\n",
                       dimension_name.c_str()); CHKERRQ(ierr);
    PISMEnd();
  }

  data.resize(2*length);                // memory allocation happens here

  ierr = this->enddef(); CHKERRQ(ierr);

  std::vector<unsigned int> start(2), count(2);
  start[0] = 0;
  start[1] = 0;
  count[0] = length;
  count[1] = 2;

  ierr = this->get_vara_double(name, start, count, &data[0]); CHKERRQ(ierr);

  // Find the corresponding 'time' variable. (We get units from the 'time'
  // variable, because according to CF-1.5 section 7.1 a "boundary variable"
  // may not have metadata set.)
  ierr = this->inq_var(dimension_name, variable_exists); CHKERRQ(ierr);

  if (variable_exists == false) {
    PetscPrintf(m_com, "PISM ERROR: Can't find '%s' in %s.\n",
                dimension_name.c_str(), filename.c_str());
    PISMEnd();
  }

  bool input_has_units = false;
  std::string input_units_string;
  Unit input_units(internal_units.get_system(), "1");

  ierr = this->get_att_text(dimension_name, "units", input_units_string); CHKERRQ(ierr);
  input_units_string = time->CF_units_to_PISM_units(input_units_string);

  if (input_units_string.empty() == true) {
    input_has_units = false;
  } else {
    try {
      input_units = Unit(internal_units.get_system(), input_units_string);
    }
    catch (RuntimeError &e) {
      std::string message = "processing variable " + name + " in " + this->inq_filename();
      e.add_context(message);
      throw;
    }
    input_has_units = true;
  }

  if (metadata.has_attribute("units") && input_has_units == false) {
    std::string units_string = internal_units.format();
    ierr = verbPrintf(2, m_com,
                      "PISM WARNING: Variable '%s' does not have the units attribute.\n"
                      "              Assuming that it is in '%s'.\n",
                      dimension_name.c_str(),
                      units_string.c_str()); CHKERRQ(ierr);
    input_units = internal_units;
  }

  ierr = convert_doubles(&data[0], data.size(), input_units, internal_units); CHKERRQ(ierr);

  // FIXME: check that time intervals described by the time bounds
  // variable are contiguous (without gaps) and stop if they are not.

  return 0;
}

PetscErrorCode PIO::write_time_bounds(const NCTimeBounds &metadata,
                                      size_t t_start,
                                      std::vector<double> &data, IO_Type nctype) const {
  PetscErrorCode ierr;
  bool variable_exists = false;

  ierr = this->inq_var(metadata.get_name(), variable_exists); CHKERRQ(ierr);
  if (variable_exists == false) {
    ierr = metadata.define(*this, nctype, true); CHKERRQ(ierr);
  }

  // convert to glaciological units:
  ierr = convert_doubles(&data[0], data.size(),
                         metadata.get_units(), metadata.get_glaciological_units()); CHKERRQ(ierr);

  ierr = this->enddef(); CHKERRQ(ierr);

  std::vector<unsigned int> start(2), count(2);
  start[0] = static_cast<unsigned int>(t_start);
  start[1] = 0;
  count[0] = static_cast<unsigned int>(data.size()) / 2;
  count[1] = 2;

  ierr = this->put_vara_double(metadata.get_name(), start, count, &data[0]); CHKERRQ(ierr);

  // restore internal units:
  ierr = convert_doubles(&data[0], data.size(),
                         metadata.get_glaciological_units(), metadata.get_units()); CHKERRQ(ierr);

  return 0;
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
