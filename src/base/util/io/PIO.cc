// Copyright (C) 2012, 2013 PISM Authors
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

#include <petscvec.h>

#include "PIO.hh"
#include "IceGrid.hh"
#include "pism_const.hh"
#include "LocalInterpCtx.hh"
#include "NCVariable.hh"
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

static PISMNCFile* create_backend(MPI_Comm com, int rank, string mode) {
  if (mode == "netcdf3") {
    return new PISMNC3File(com, rank);
  } else if (mode.find("quilt") == 0) {
    size_t n = mode.find(":");
    int compression_level = 0;

    if (n != string::npos) {
      mode.replace(0, 6, "");
      char *endptr;
      compression_level = strtol(mode.c_str(), &endptr, 10);
      if ((*endptr != '\0') || (compression_level < 0) || (compression_level > 9)) {
        PetscPrintf(com, "PISM WARNING: invalid compression level %s. Output compression is disabled.\n",
                    mode.c_str());
        compression_level = 0;
      }
    }

    return new PISMNC4_Quilt(com, rank, compression_level);
  }
#if (PISM_USE_PARALLEL_NETCDF4==1)
  else if (mode == "netcdf4_parallel") {
    return new PISMNC4_Par(com, rank);
  }
#endif
#if (PISM_USE_PNETCDF==1)
  else if (mode == "pnetcdf") {
    return new PISMPNCFile(com, rank);
  }
#endif
#if (PISM_USE_HDF5==1)
  else if (mode == "hdf5") {
    return new PISMNC4_HDF5(com, rank);
  }
#endif
  else {
    return NULL;
  }
}

//! \brief The code shared by different PIO constructors.
void PIO::constructor(MPI_Comm c, int r, string mode, PISMUnitSystem unit_system) {
  com = c;
  rank = r;
  shallow_copy = false;
  m_mode = mode;
  m_unit_system = unit_system;

  nc = create_backend(com, rank, mode);

  if (mode != "guess_mode" && nc == NULL) {
    PetscPrintf(com,
                "PISM ERROR: output format '%s' is not supported.\n"
                "Please recompile PISM with the appropriate I/O library.\n",
                mode.c_str());
    PISMEnd();
  }
}

PIO::PIO(MPI_Comm c, int r, string mode, PISMUnitSystem units_system) {
  constructor(c, r, mode, units_system);
}

PIO::PIO(IceGrid &grid, string mode) {
  constructor(grid.com, grid.rank, mode, grid.get_unit_system());
  if (nc != NULL)
    set_local_extent(grid.xs, grid.xm, grid.ys, grid.ym);
}

PIO::PIO(const PIO &other) {
  com = other.com;
  rank = other.rank;
  nc = other.nc;
  m_mode = other.m_mode;
  m_unit_system = other.m_unit_system;

  shallow_copy = true;
}

PIO::~PIO() {
  if (shallow_copy == false)
    delete nc;
}

// Chooses the best I/O backend for reading from 'filename'.
PetscErrorCode PIO::detect_mode(string filename) {
  int stat;

  // Bail if someone made a decision already.
  if (nc != NULL)
    return 1;

  PISMNC3File nc3(com, rank);

  stat = nc3.open(filename, PISM_NOWRITE);
  if (stat != 0) {
    PetscPrintf(com, "PISM ERROR: Can't open '%s'. Exiting...\n", filename.c_str());
    PISMEnd();
  }

  string format = nc3.get_format();

  stat = nc3.close(); CHKERRQ(stat);

  vector<string> modes;
  if (format == "netcdf4") {
    modes.push_back("netcdf4_parallel");
    modes.push_back("netcdf3");
  } else {
    modes.push_back("pnetcdf");
    modes.push_back("netcdf3");
  }

  for (unsigned int j = 0; j < modes.size(); ++j) {
    nc = create_backend(com, rank, modes[j]);

    if (nc != NULL) {
      m_mode = modes[j];
      stat = verbPrintf(3, com,
                        "  - Using the %s backend to read from %s...\n",
                        modes[j].c_str(), filename.c_str()); CHKERRQ(stat);
      break;
    }
  }

  if (nc == NULL) {
    PetscPrintf(com, "PISM ERROR: Unable to allocate an I/O backend. This should never happen!\n");
    PISMEnd();
  }

  nc->set_local_extent(m_xs, m_xm, m_ys, m_ym);

  return 0;
}

PetscErrorCode PIO::open(string filename, int mode, bool append) {
  PetscErrorCode stat;

  // opening for reading

  if (!(mode & PISM_WRITE)) {
    if (nc == NULL && m_mode == "guess_mode") {
      stat = detect_mode(filename); CHKERRQ(stat);
    }

    assert(nc != NULL);

    stat = nc->open(filename, mode);
    if (stat != 0) {
      PetscPrintf(com, "PISM ERROR: Can't open '%s'. Exiting...\n", filename.c_str());
      PISMEnd();
    }
    return 0;
  }

  // opening for writing

  if (append == false) {
    assert(nc != NULL);

    stat = nc->move_if_exists(filename); CHKERRQ(stat);

    stat = nc->create(filename);
    if (stat != 0) {
      PetscPrintf(com, "PISM ERROR: Can't create '%s'. Exiting...\n", filename.c_str());
      PISMEnd();
    }

    int old_fill;
    stat = nc->set_fill(PISM_NOFILL, old_fill); CHKERRQ(stat);
  } else {
    if (nc == NULL && m_mode == "guess_mode") {
      stat = detect_mode(filename); CHKERRQ(stat);
    }

    assert(nc != NULL);

    stat = nc->open(filename, mode);

    if (stat != 0) {
      PetscPrintf(com, "PISM ERROR: Can't open '%s' (rank = %d). Exiting...\n", filename.c_str(), rank);
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

string PIO::inq_filename() const {
  return nc->get_filename();
}


//! \brief Get the number of records. Uses the length of an unlimited dimension.
PetscErrorCode PIO::inq_nrecords(unsigned int &result) const {
  PetscErrorCode ierr;
  string dim;

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
PetscErrorCode PIO::inq_nrecords(string name, string std_name, unsigned int &result) const {
  PetscErrorCode ierr;

  bool exists = false, found_by_standard_name = false;
  string name_found;

  ierr = this->inq_var(name, std_name, exists, name_found, found_by_standard_name); CHKERRQ(ierr);

  if (exists == false) {
    result = 0;
    return 0;
  }

  vector<string> dims;
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
PetscErrorCode PIO::inq_var(string short_name, string std_name, bool &exists,
                            string &result, bool &found_by_standard_name) const {
  PetscErrorCode ierr;

  exists = false;

  if (std_name.empty() == false) {
    int nvars;

    ierr = nc->inq_nvars(nvars); CHKERRQ(ierr);

    for (int j = 0; j < nvars; ++j) {
      string name, attribute;
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
	  ierr = PetscPrintf(com,
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
PetscErrorCode PIO::inq_var(string name, bool &exists) const {

  PetscErrorCode ierr = nc->inq_varid(name, exists); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode PIO::inq_vardims(string name, vector<string> &result) const {
  PetscErrorCode ierr = nc->inq_vardimid(name, result); CHKERRQ(ierr);
  return 0;
}


//! \brief Checks if a dimension exists.
PetscErrorCode PIO::inq_dim(string name, bool &exists) const {

  PetscErrorCode ierr = nc->inq_dimid(name, exists); CHKERRQ(ierr);

  return 0;
}

//! \brief Get the length of a dimension.
/*!
 * Sets result to 0 if a dimension does not exist.
 */
PetscErrorCode PIO::inq_dimlen(string name, unsigned int &result) const {
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
PetscErrorCode PIO::inq_dimtype(string name, AxisType &result) const {
  PetscErrorCode ierr;
  string axis, standard_name, units;
  PISMUnit tmp_units;
  bool exists;

  ierr = nc->inq_varid(name, exists); CHKERRQ(ierr);

  if (exists == false) {
    PetscPrintf(com, "ERROR: coordinate variable '%s' is not present!\n", name.c_str());
    PISMEnd();
  }

  ierr = nc->get_att_text(name, "axis", axis); CHKERRQ(ierr);
  ierr = nc->get_att_text(name, "standard_name", standard_name); CHKERRQ(ierr);
  ierr = nc->get_att_text(name, "units", units); CHKERRQ(ierr);

  // check if it has units compatible with "seconds":

  if (tmp_units.parse(m_unit_system, units) != 0) {
    ierr = PetscPrintf(com, "ERROR: units specification '%s' is unknown or invalid (processing variable '%s').\n",
		       units.c_str(), name.c_str());
    PISMEnd();
  }

  PISMUnit seconds;
  int errcode = seconds.parse(m_unit_system, "seconds");
  assert(errcode == 0);
  if (ut_are_convertible(tmp_units.get(), seconds.get())) {
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

PetscErrorCode PIO::inq_dim_limits(string name, double *min, double *max) const {
  PetscErrorCode ierr;
  vector<double> data;

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
PetscErrorCode PIO::inq_grid(string var_name, IceGrid *grid, Periodicity periodicity) const {
  PetscErrorCode ierr;

  if (grid == NULL)
    SETERRQ(com, 1, "grid == NULL");

  grid_info input(m_unit_system);

  // The following call may fail because var_name does not exist. (And this is fatal!)
  ierr = this->inq_grid_info(var_name, input); CHKERRQ(ierr);

  // if we have no vertical grid information, create a fake 2-level vertical grid.
  if (input.z.size() < 2) {
    double Lz = grid->config.get("grid_Lz");
    ierr = verbPrintf(3, com,
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

  // The grid dimensions Lx/Ly are computed differently depending on the grid's
  // periodicity. For x-periodic grids, e.g., the length Lx is a little longer
  // than the difference between the maximum and minimum x-coordinates.
  PetscReal x_max = input.x_max, x_min = input.x_min;
  if (periodicity & X_PERIODIC) {
    PetscReal dx = (x_max-x_min)/(grid->Mx-1);
    x_max += dx;
  }

  PetscReal y_max = input.y_max, y_min = input.y_min;
  if (periodicity & Y_PERIODIC) {
    PetscReal dy = (y_max-y_min)/(grid->My-1);
    y_max += dy;
  }

  grid->x0 = (x_max + x_min) / 2.0;
  grid->y0 = (y_max + y_min) / 2.0;
  grid->Lx = (x_max - x_min) / 2.0;
  grid->Ly = (y_max - y_min) / 2.0;

  grid->time->set_start(input.time);
  ierr = grid->time->init(); CHKERRQ(ierr); // re-initialize to take the new start time into account

  ierr = grid->compute_horizontal_spacing(); CHKERRQ(ierr);
  ierr = grid->set_vertical_levels(input.z); CHKERRQ(ierr);

  // We're ready to call grid->allocate().

  return 0;
}


PetscErrorCode PIO::inq_units(string name, bool &has_units, PISMUnit &units,
                              bool use_reference_date) const {
  PetscErrorCode ierr;
  string units_string;

  // Get the string:
  ierr = nc->get_att_text(name, "units", units_string); CHKERRQ(ierr);

  // If a variables does not have the units attribute, set the flag and return:
  if (units_string.empty()) {
    has_units = false;
    units.reset();
    return 0;
  }

  int n = (int)units_string.find("since");
  if (use_reference_date == false) {
    /*!
      \note This code finds the string "since" in the units_string and
      terminates it on the first 's' of "since", if this sub-string was found.
      This is done to ignore the reference date in the time units string (the
      reference date specification always starts with this word).
    */
    if (n != -1) units_string.resize(n);
  } else if (n == -1) {
    ierr = PetscPrintf(com,
                       "PISM ERROR: units specification '%s' does not contain a reference date (processing variable '%s').\n",
                       units_string.c_str(), name.c_str());
    PISMEnd();
  }

  // strip trailing spaces
  while (ends_with(units_string, " "))
    units_string.resize(units_string.size() - 1);
  
  if (units.parse(m_unit_system, units_string) != 0) {
    ierr = PetscPrintf(com, "PISM ERROR: units specification '%s' is unknown or invalid (processing variable '%s').\n",
		       units_string.c_str(), name.c_str());
    PISMEnd();
  }

  has_units = true;

  return 0;
}


PetscErrorCode PIO::inq_grid_info(string name, grid_info &g) const {
  PetscErrorCode ierr;

  vector<string> dims;
  bool exists, found_by_standard_name;
  string name_found;

  // try "name" as the standard_name first, then as the short name:
  ierr = this->inq_var(name, name, exists, name_found, found_by_standard_name); CHKERRQ(ierr);

  if (exists == false)
    SETERRQ2(com, 1, "Could not find variable %s in %s", name.c_str(),
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
    string dimname = dims[i];

    AxisType dimtype = UNKNOWN_AXIS;
    ierr = this->inq_dimtype(dimname, dimtype); CHKERRQ(ierr);

    switch (dimtype) {
    case X_AXIS:
      {
        ierr = nc->inq_dimlen(dimname, g.x_len); CHKERRQ(ierr);
        ierr = this->inq_dim_limits(dimname, &g.x_min, &g.x_max); CHKERRQ(ierr);
        ierr = this->get_dim(dimname, g.x); CHKERRQ(ierr);
        break;
      }
    case Y_AXIS:
      {
        ierr = nc->inq_dimlen(dimname, g.y_len); CHKERRQ(ierr);
        ierr = this->inq_dim_limits(dimname, &g.y_min, &g.y_max); CHKERRQ(ierr);
        ierr = this->get_dim(dimname, g.y); CHKERRQ(ierr);
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
        PetscPrintf(com, "ERROR: Can't figure out which direction dimension '%s' corresponds to.",
                    dimname.c_str());
        PISMEnd();
      }
    } // switch
  }   // for loop

  return 0;
}

//! \brief Define a dimension \b and the associated coordinate variable. Set string attributes.
PetscErrorCode PIO::def_dim(string name, long int length, map<string,string> attrs) const {
  PetscErrorCode ierr;

  ierr = nc->redef();
  CHKERRQ(ierr);

  ierr = nc->def_dim(name, length); CHKERRQ(ierr);

  vector<string> dims(1);
  dims[0] = name;

  ierr = nc->def_var(name, PISM_DOUBLE, dims); CHKERRQ(ierr);

  map<string,string>::iterator j;
  for (j = attrs.begin(); j != attrs.end(); ++j) {
    const string &att_name = j->first,
      &att_value = j->second;

    ierr = nc->put_att_text(name, att_name, att_value); CHKERRQ(ierr);
  }

  return 0;
}

//! \brief Define a variable.
PetscErrorCode PIO::def_var(string name, PISM_IO_Type nctype, vector<string> dims) const {
  PetscErrorCode ierr;

  ierr = nc->def_var(name, nctype, dims); CHKERRQ(ierr);

  return 0;
}
PetscErrorCode PIO::get_1d_var(string name, unsigned int s, unsigned int c,
                               vector<double> &result) const {
  PetscErrorCode ierr;
  vector<unsigned int> start(1), count(1);

  result.resize(c);

  start[0] = s;
  count[0] = c;

  ierr = nc->enddef(); CHKERRQ(ierr);

  ierr = nc->get_vara_double(name, start, count, &result[0]); CHKERRQ(ierr);

  return 0;
}


PetscErrorCode PIO::put_1d_var(string name, unsigned int s, unsigned int c,
                               const vector<double> &data) const {
  PetscErrorCode ierr;
  vector<unsigned int> start(1), count(1);

  start[0] = s;
  count[0] = c;

  ierr = nc->enddef(); CHKERRQ(ierr);

  ierr = nc->put_vara_double(name, start, count,
                             const_cast<double*>(&data[0])); CHKERRQ(ierr);

  return 0;
}


//! \brief Get dimension data (a coordinate variable).
PetscErrorCode PIO::get_dim(string name, vector<double> &data) const {
  PetscErrorCode ierr;

  unsigned int dim_length = 0;
  ierr = nc->inq_dimlen(name, dim_length); CHKERRQ(ierr);

  ierr = nc->enddef(); CHKERRQ(ierr);

  ierr = this->get_1d_var(name, 0, dim_length, data); CHKERRQ(ierr);

  return 0;
}

//! \brief Write dimension data (a coordinate variable).
PetscErrorCode PIO::put_dim(string name, const vector<double> &data) const {
  PetscErrorCode ierr = this->put_1d_var(name, 0,
                                         (unsigned int)data.size(), data); CHKERRQ(ierr);
  return 0;
}

PetscErrorCode PIO::def_time(string name, string calendar, string units) const {
  PetscErrorCode ierr;
  map<string,string> attrs;

  bool time_exists;
  ierr = nc->inq_varid(name, time_exists); CHKERRQ(ierr);
  if (time_exists)
    return 0;

  // t
  attrs["long_name"] = "time";
  attrs["calendar"]  = calendar;
  attrs["units"]     = units;
  attrs["axis"]      = "T";

  ierr = this->def_dim(name, PISM_UNLIMITED, attrs); CHKERRQ(ierr);

  return 0;
}


//! \brief Append to the time dimension.
PetscErrorCode PIO::append_time(string name, double value) const {
  PetscErrorCode ierr;

  vector<unsigned int> start(1), count(1);
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
PetscErrorCode PIO::append_history(string history) const {
  PetscErrorCode ierr;
  string old_history;

  ierr = nc->redef(); CHKERRQ(ierr);

  ierr = nc->get_att_text("PISM_GLOBAL", "history", old_history); CHKERRQ(ierr);
  ierr = nc->put_att_text("PISM_GLOBAL", "history", history + old_history); CHKERRQ(ierr);

  return 0;
}

//! \brief Write a multiple-valued double attribute.
PetscErrorCode PIO::put_att_double(string var_name, string att_name, PISM_IO_Type nctype,
                                   vector<double> values) const {
  PetscErrorCode ierr;

  ierr = nc->redef(); CHKERRQ(ierr);

  ierr = nc->put_att_double(var_name, att_name, nctype, values); CHKERRQ(ierr);

  return 0;
}

//! \brief Write a single-valued double attribute.
PetscErrorCode PIO::put_att_double(string var_name, string att_name, PISM_IO_Type nctype,
                                   double value) const {
  PetscErrorCode ierr;
  vector<double> tmp; tmp.push_back(value);

  ierr = nc->redef(); CHKERRQ(ierr);

  ierr = nc->put_att_double(var_name, att_name, nctype, tmp); CHKERRQ(ierr);

  return 0;
}

//! \brief Write a text attribute.
PetscErrorCode PIO::put_att_text(string var_name, string att_name, string value) const {
  PetscErrorCode ierr;

  ierr = nc->redef(); CHKERRQ(ierr);

  string tmp = value + "\0";    // ensure that the string is null-terminated

  ierr = nc->put_att_text(var_name, att_name, tmp); CHKERRQ(ierr);

  return 0;
}

//! \brief Get a double attribute.
PetscErrorCode PIO::get_att_double(string var_name, string att_name,
                                   vector<double> &result) const {

  PetscErrorCode ierr;
  PISM_IO_Type att_type;
  // virtual int inq_atttype(string variable_name, string att_name, PISM_IO_Type &result) const = 0;

  ierr = nc->inq_atttype(var_name, att_name, att_type); CHKERRQ(ierr);

  // Give an understandable error message if a string attribute was found when
  // a number (or a list of numbers) was expected. (We've seen datasets with
  // "valid_min" stored as a string...)
  if (att_type == PISM_CHAR) {
    string tmp;
    ierr = nc->get_att_text(var_name, att_name, tmp); CHKERRQ(ierr);

    PetscPrintf(com,
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
PetscErrorCode PIO::get_att_text(string var_name, string att_name, string &result) const {

  PetscErrorCode ierr = nc->get_att_text(var_name, att_name, result); CHKERRQ(ierr);

  return 0;
}

//! \brief Read a PETSc Vec using the grid "grid".
/*!
 * Assumes that PetscScalar corresponds to C++ double.
 *
 * Vec g has to be "global" (i.e. without ghosts).
 */
PetscErrorCode PIO::get_vec(IceGrid *grid, string var_name, unsigned int z_count, int t, Vec g) const {
  PetscErrorCode ierr;

  vector<unsigned int> start, count, imap;
  ierr = compute_start_and_count(var_name,
                                 t,
                                 grid->xs, grid->xm,
                                 grid->ys, grid->ym,
                                 0, z_count,
                                 start, count, imap); CHKERRQ(ierr);

  ierr = nc->enddef(); CHKERRQ(ierr);

  PetscScalar *a_petsc;
  ierr = VecGetArray(g, &a_petsc); CHKERRQ(ierr);

  // We always use "mapped" I/O here, because we don't know where the input
  // file came from.
  ierr = nc->get_varm_double(var_name, start, count, imap, (double*)a_petsc); CHKERRQ(ierr);

  ierr = VecRestoreArray(g, &a_petsc); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode PIO::inq_nattrs(string var_name, int &result) const {
  PetscErrorCode ierr = nc->inq_varnatts(var_name, result); CHKERRQ(ierr);
  return 0;
}


PetscErrorCode PIO::inq_attname(string var_name, unsigned int n, string &result) const {
  PetscErrorCode ierr = nc->inq_attname(var_name, n, result); CHKERRQ(ierr);
  return 0;
}


PetscErrorCode PIO::inq_atttype(string var_name, string att_name, PISM_IO_Type &result) const {
  PetscErrorCode ierr = nc->inq_atttype(var_name, att_name, result); CHKERRQ(ierr);
  return 0;
}



//! \brief Write a PETSc Vec using the grid "grid".
/*!
 * Assumes that PetscScalar corresponds to C++ double.
 *
 * Vec g has to be "global" (i.e. without ghosts).
 *
 * This method always writes to the last record in the file.
 */
PetscErrorCode PIO::put_vec(IceGrid *grid, string var_name, unsigned int z_count, Vec g) const {
  PetscErrorCode ierr;

  unsigned int t;
  ierr = nc->inq_dimlen(grid->config.get_string("time_dimension_name"), t); CHKERRQ(ierr);

#if (PISM_DEBUG==1)
  if (t < 1)
    fprintf(stderr, "time dimension length (%d) is less than 1!\n", t);
#endif

  vector<unsigned int> start, count, imap;
  ierr = compute_start_and_count(var_name,
                                 t - 1,
                                 grid->xs, grid->xm,
                                 grid->ys, grid->ym,
                                 0, z_count,
                                 start, count, imap); CHKERRQ(ierr);

  ierr = nc->enddef(); CHKERRQ(ierr);

  PetscScalar *a_petsc;
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

//! \brief Read a PETSc Vec from a file, using bilinear (or trilinear)
//! interpolation to put it on the grid defined by "grid" and zlevels_out.
PetscErrorCode PIO::regrid_vec(IceGrid *grid, string var_name,
                               const vector<double> &zlevels_out,
                               LocalInterpCtx *lic, Vec g) const {
  PetscErrorCode ierr;
  const int T = 0, X = 1, Y = 2, Z = 3; // indices, just for clarity
  vector<unsigned int> start, count, imap;

  int t_start = lic->start[T],
    x_start = lic->start[X],
    y_start = lic->start[Y],
    x_count = lic->count[X],
    y_count = lic->count[Y],
    z_start = lic->start[Z],
    z_count = lic->count[Z];

  ierr = compute_start_and_count(var_name, t_start,
                                 x_start, x_count,
                                 y_start, y_count,
                                 z_start, z_count,
                                 start, count, imap); CHKERRQ(ierr);

  ierr = nc->enddef(); CHKERRQ(ierr);

  // We always use "mapped" I/O here, because we don't know where the input
  // file came from.
  ierr = nc->get_varm_double(var_name, start, count, imap, lic->a); CHKERRQ(ierr);

  ierr = regrid(grid, zlevels_out, lic, g);

  return 0;
}

int PIO::k_below(double z, const vector<double> &zlevels) const {
  double z_min = zlevels.front(), z_max = zlevels.back();
  PetscInt mcurr = 0;

  if (z < z_min - 1.0e-6 || z > z_max + 1.0e-6) {
    PetscPrintf(com,
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
PetscErrorCode PIO::regrid(IceGrid *grid, const vector<double> &zlevels_out, LocalInterpCtx *lic, Vec g) const {
  const int Y = 2, Z = 3; // indices, just for clarity
  PetscErrorCode ierr;

  vector<double> &zlevels_in = lic->zlevels;
  unsigned int nlevels = (int)zlevels_out.size();
  double *input_array = lic->a;

  // array sizes for mapping from logical to "flat" indices
  int y_count = lic->count[Y],
    z_count = lic->count[Z];

  // We'll work with the raw storage here so that the array we are filling is
  // indexed the same way as the buffer we are pulling from (input_array)
  PetscScalar *output_array;
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


PetscErrorCode PIO::compute_start_and_count(string short_name, int t_start,
                                            int x_start, int x_count,
                                            int y_start, int y_count,
                                            int z_start, int z_count,
                                            vector<unsigned int> &start,
                                            vector<unsigned int> &count,
                                            vector<unsigned int> &imap) const {
  PetscErrorCode ierr;
  vector<string> dims;

  ierr = nc->inq_vardimid(short_name, dims); CHKERRQ(ierr);
  int ndims = (int)dims.size();

  // Resize output vectors:
  start.resize(ndims);
  count.resize(ndims);
  imap.resize(ndims);

  // Assemble start, count and imap:
  for (int j = 0; j < ndims; j++) {
    string dimname = dims[j];

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
    case Z_AXIS:
      start[j] = z_start;
      count[j] = z_count;
      imap[j]  = 1;
      break;
    default:
      {
        SETERRQ(com, 1, "dimtype is not one of T_AXIS, X_AXIS, Y_AXIS, Z_AXIS");
      }
    }

// #if (PISM_DEBUG==1)
//     fprintf(stderr, "[%d] var=%s start[%d]=%d count[%d]=%d imap[%d]=%d\n",
//             rank, short_name.c_str(),
//             j, start[j], j, count[j], j, imap[j]);
// #endif

  }

  return 0;
}

PetscErrorCode PIO::get_vara_double(string variable_name,
                                    vector<unsigned int> start,
                                    vector<unsigned int> count,
                                    double *ip) const {
  PetscErrorCode ierr;

  ierr = nc->enddef(); CHKERRQ(ierr);

  ierr = nc->get_vara_double(variable_name, start, count, ip); CHKERRQ(ierr);

  return 0;
}


PetscErrorCode PIO::put_vara_double(string variable_name,
                                    vector<unsigned int> start,
                                    vector<unsigned int> count,
                                    double *op) const {
  PetscErrorCode ierr;

  ierr = nc->enddef(); CHKERRQ(ierr);

  ierr = nc->put_vara_double(variable_name, start, count, op); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode PIO::get_varm_double(string variable_name,
                                    vector<unsigned int> start,
                                    vector<unsigned int> count,
                                    vector<unsigned int> imap, double *ip) const {
  PetscErrorCode ierr;

  ierr = nc->enddef(); CHKERRQ(ierr);

  ierr = nc->get_varm_double(variable_name, start, count, imap, ip); CHKERRQ(ierr);

  return 0;
}


PetscErrorCode PIO::put_varm_double(string variable_name,
                                    vector<unsigned int> start,
                                    vector<unsigned int> count,
                                    vector<unsigned int> imap, double *op) const {
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
