// Copyright (C) 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019 PISM Authors
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

#include <cassert>
#include <cstdio>
#include <memory>
#include <map>
using std::shared_ptr;

#include <petscvec.h>

#include "File.hh"
#include "pism/util/IceGrid.hh"
#include "pism/util/pism_utilities.hh"
#include "pism/util/VariableMetadata.hh"
#include "pism/util/ConfigInterface.hh"
#include "pism/util/Time.hh"
#include "NC3File.hh"

#include "pism/pism_config.hh"

#if (Pism_USE_PARALLEL_NETCDF4==1)
#include "NC4_Par.hh"
#endif

#if (Pism_USE_PNETCDF==1)
#include "PNCFile.hh"
#endif

#if (Pism_USE_PIO==1)
#include "ParallelIO.hh"
#endif

#include "CDI.hh"

#include "pism/util/error_handling.hh"
#include "pism/util/io/io_helpers.hh"

namespace pism {

struct File::Impl {
  MPI_Comm com;
  IO_Backend backend;
  io::NCFile::Ptr nc;
};

IO_Backend string_to_backend(const std::string &backend) {
  if (backend == "netcdf3") {
    return PISM_NETCDF3;
  }
  if (backend == "netcdf4_parallel") {
    return PISM_NETCDF4_PARALLEL;
  }
  if (backend == "pnetcdf") {
    return PISM_PNETCDF;
  }
  if (backend == "pio_pnetcdf") {
    return PISM_PIO_PNETCDF;
  }
  if (backend == "pio_netcdf") {
    return PISM_PIO_NETCDF;
  }
  if (backend == "pio_netcdf4c") {
    return PISM_PIO_NETCDF4C;
  }
  if (backend == "pio_netcdf4p") {
    return PISM_PIO_NETCDF4P;
  }
  if (backend == "cdi") {
    return PISM_CDI;
  }
  throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                "unknown or unsupported I/O backend: %s", backend.c_str());
}

// Chooses the best available I/O backend for reading from 'filename'.
static IO_Backend choose_backend(MPI_Comm com, const std::string &filename) {

  std::string format;
  {
    // This is the rank-0-only purely-serial mode of accessing NetCDF files, but it
    // supports all the kinds of NetCDF, so this is fine.
    io::NC3File file(com);

    file.open(filename, PISM_READONLY);
    format = file.get_format();
    file.close();
  }

  if (format == "netcdf4") {
#if (Pism_USE_PARALLEL_NETCDF4==1)
    return PISM_NETCDF4_PARALLEL;
#endif
  } else {
#if (Pism_USE_PNETCDF==1)
    return PISM_PNETCDF;
#endif
  }

  // this choice is appropriate for both NetCDF-3 and NetCDF-4
  return PISM_NETCDF3;
}

static io::NCFile::Ptr create_backend(MPI_Comm com, IO_Backend backend, int iosysid) {
  int size = 1;
  MPI_Comm_size(com, &size);

  if (backend == PISM_NETCDF3) {
    return io::NCFile::Ptr(new io::NC3File(com));
  }
#if (Pism_USE_PARALLEL_NETCDF4==1)
  if (backend == PISM_NETCDF4_PARALLEL) {
    return io::NCFile::Ptr(new io::NC4_Par(com));
  }
#endif
#if (Pism_USE_PNETCDF==1)
  if (backend == PISM_PNETCDF) {
    return io::NCFile::Ptr(new io::PNCFile(com));
  }
#endif
#if (Pism_USE_PIO==1)
  if (backend == PISM_PIO_PNETCDF or
      backend == PISM_PIO_NETCDF4P or
      backend == PISM_PIO_NETCDF4C or
      backend == PISM_PIO_NETCDF) {
    if (iosysid == -1) {
      throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                    "To use ParallelIO you have to pass iosysid to File");
    }
    return io::NCFile::Ptr(new io::ParallelIO(com, iosysid, backend));
  }
#else
  (void) iosysid;               // silence a compiler warning
#endif
//#if (Pism_USE_CDI==1)
  if (backend == PISM_CDI) {
    return io::NCFile::Ptr(new io::CDI(com));
  }
//#endif
  throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                "unknown or unsupported I/O backend: %d", backend);
}

File::File(MPI_Comm com, const std::string &filename, IO_Backend backend, IO_Mode mode,
           int iosysid, const std::map<std::string, int> &varsi, const std::vector<int>& gridIDs, int FileID)
  : m_impl(new Impl) {

  if (filename.empty()) {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                  "cannot open file: provided file name is empty");
  }

  if (backend == PISM_GUESS) {
    m_impl->backend = choose_backend(com, filename);
  } else {
    m_impl->backend = backend;
  }

  m_impl->com = com;
  m_impl->nc  = create_backend(m_impl->com, m_impl->backend, iosysid);
  this->set_gridIDs(gridIDs);
  this->open(filename, mode, varsi, FileID);
}

File::~File() {
  if (m_impl->nc and not filename().empty()) {
    try {
      // a file is still open, so we try to close it
      this->close();
    } catch (...) {
      // don't ever throw from here
      handle_fatal_errors(MPI_COMM_SELF);
    }
  }
  delete m_impl;
}

void File::set_gridIDs(const std::vector<int>& gridIDs) const {
 m_impl->nc->set_ncgridIDs(gridIDs); 
}

std::vector<int> File::get_gridIDs() const {
 return m_impl->nc->get_ncgridIDs();
}

int File::get_streamID() const {
  return m_impl->nc->get_ncstreamID();
}

int File::get_vlistID() const {
  return m_impl->nc->get_ncvlistID();
}

void File::set_dimatt() const {
  m_dimatt[std::string("x")] = std::string("not written");
  m_dimatt[std::string("y")] = std::string("not written");
  m_dimatt[std::string("zb")] = std::string("not written");
  m_dimatt[std::string("z")] = std::string("not written");
}

std::string File::get_dimatt_value(std::string &dim_name) const {
  return m_dimatt[dim_name];
}

void File::set_dimatt_value(std::string &dim_name, std::string &dimatt_value) const {
  m_dimatt[dim_name] = dimatt_value;
}

MPI_Comm File::com() const {
  return m_impl->com;
}

IO_Backend File::backend() const {
  return m_impl->backend;
}

void File::open(const std::string &filename, IO_Mode mode, const std::map<std::string, int> &varsi, int FileID) {
  try {

    // opening for reading
    if (mode == PISM_READONLY) {

      m_impl->nc->open(filename, mode, varsi);

    } else if (mode == PISM_READWRITE_CLOBBER or mode == PISM_READWRITE_MOVE) {

      if (mode == PISM_READWRITE_MOVE) {
        io::move_if_exists(m_impl->com, filename);
      } else {
        io::remove_if_exists(m_impl->com, filename);
      }

      m_impl->nc->create(filename, FileID);

      int old_fill;
      m_impl->nc->set_fill(PISM_NOFILL, old_fill);
    } else if (mode == PISM_READWRITE) {                      // mode == PISM_READWRITE

      m_impl->nc->open(filename, mode, varsi, FileID);

      int old_fill;
      m_impl->nc->set_fill(PISM_NOFILL, old_fill);
    } else {
      throw RuntimeError::formatted(PISM_ERROR_LOCATION, "invalid mode: %d", mode);
    }
  } catch (RuntimeError &e) {
    e.add_context("opening or creating \"" + filename + "\"");
    throw;
  }
}

void File::remove_attribute(const std::string &variable_name, const std::string &att_name) const {
  try {
    m_impl->nc->del_att(variable_name, att_name);
  } catch (RuntimeError &e) {
    e.add_context("deleting the attribute %s:%s", variable_name.c_str(), att_name.c_str());
    throw;
  }
}

void File::close() {
  try {
    m_impl->nc->close();
  } catch (RuntimeError &e) {
    e.add_context("closing \"" + filename() + "\"");
    throw;
  }
}

void File::sync() const {
  try {
    m_impl->nc->sync();
  } catch (RuntimeError &e) {
    e.add_context("synchronizing \"" + filename() + "\"");
    throw;
  }
}

void File::redef() const {
  try {
    m_impl->nc->redef();
  } catch (RuntimeError &e) {
    e.add_context("switching to define mode; file \"" + filename() + "\"");
    throw;
  }
}


void File::enddef() const {
  try {
    m_impl->nc->enddef();
  } catch (RuntimeError &e) {
    e.add_context("switching to data mode; file \"" + filename() + "\"");
    throw;
  }
}

std::string File::filename() const {
  return m_impl->nc->filename();
}


//! \brief Get the number of records. Uses the length of an unlimited dimension.
unsigned int File::nrecords() const {
  try {
    std::string dim;
    m_impl->nc->inq_unlimdim(dim);

    if (dim.empty()) {
      return 1;                 // one record
    } else {
      return this->dimension_length(dim);
    }
  } catch (RuntimeError &e) {
    e.add_context("getting the number of records in file \"" + filename() + "\"");
    throw;
  }
  return 0;                     // LCOV_EXCL_LINE
}

//! \brief Get the number of records of a certain variable. Uses the length of
//! an associated "time" dimension.
unsigned int File::nrecords(const std::string &name, const std::string &std_name,
                            units::System::Ptr unit_system) const {
  try {
    auto var = find_variable(name, std_name);

    if (not var.exists) {
      return 0;
    }

    auto dims = dimensions(var.name);

    for (auto d : dims) {
      AxisType dimtype = dimension_type(d, unit_system);

      if (dimtype == T_AXIS) {
        return this->dimension_length(d);
      }
    }

    return 1;                   // one record
  } catch (RuntimeError &e) {
    e.add_context("getting the number of records of variable '%s' ('%s') in '%s'",
                  name.c_str(), std_name.c_str(), filename().c_str());
    throw;
  }
  return 0;                     // LCOV_EXCL_LINE
}


//! \brief Find a variable using its standard name and/or short name.
/*!
 * Sets "result" to the short name found.
 */
VariableLookupData File::find_variable(const std::string &short_name, const std::string &std_name) const {
  VariableLookupData result;
  try {
    result.exists = false;

    if (not std_name.empty()) {

      int n_variables = nvariables();

      for (int j = 0; j < n_variables; ++j) {
        std::string
          name      = variable_name(j),
          attribute = read_text_attribute(name, "standard_name");

        if (attribute.empty()) {
          continue;
        }

        if (attribute == std_name) {
          if (not result.exists) {
            result.exists = true;
            result.found_using_standard_name = true;
            result.name = name;
          } else {
            throw RuntimeError::formatted(PISM_ERROR_LOCATION, "inconsistency in '%s': variables '%s' and '%s'\n"
                                          "have the same standard_name (%s)",
                                          filename().c_str(), result.name.c_str(),
                                          name.c_str(), attribute.c_str());
          }
        }

      } // end of the for loop
    } // end of if (not std_name.empty())

    if (not result.exists) {
      m_impl->nc->inq_varid(short_name, result.exists);
      if (result.exists) {
        result.name = short_name;
      } else {
        result.name = "";
      }

      result.found_using_standard_name = false;
    }

  } catch (RuntimeError &e) {
    e.add_context("searching for variable '%s' ('%s') in '%s'", short_name.c_str(), std_name.c_str(), filename().c_str());
    throw;
  }

  return result;
}

//! \brief Checks if a variable exists.
bool File::find_variable(const std::string &name) const {
  try {
    bool exists = false;
    m_impl->nc->inq_varid(name, exists);
    return exists;
  } catch (RuntimeError &e) {
    e.add_context("searching for variable '%s' in '%s'", name.c_str(), filename().c_str());
    throw;
  }
}

std::vector<std::string> File::dimensions(const std::string &variable_name) const {
  try {
    std::vector<std::string> result;
    m_impl->nc->inq_vardimid(variable_name, result);
    return result;
  } catch (RuntimeError &e) {
    e.add_context("getting dimensions of variable '%s' in '%s'", variable_name.c_str(),
                  filename().c_str());
    throw;
  }
}


//! \brief Checks if a dimension exists.
bool File::find_dimension(const std::string &name) const {
  try {
    bool exists = false;
    m_impl->nc->inq_dimid(name, exists);
    return exists;
  } catch (RuntimeError &e) {
    e.add_context("searching for dimension '%s' in '%s'", name.c_str(), filename().c_str());
    throw;
  }
}

//! \brief Get the length of a dimension.
/*!
 * Sets result to 0 if a dimension does not exist.
 */
unsigned int File::dimension_length(const std::string &name) const {
  try {
    if (find_dimension(name)) {
      unsigned int result = 0;
      m_impl->nc->inq_dimlen(name, result);
      return result;
    } else {
      return 0;
    }
  } catch (RuntimeError &e) {
    e.add_context("getting the length of dimension '%s' in '%s'", name.c_str(), filename().c_str());
    throw;
  }
}

//! \brief Get the "type" of a dimension.
/*!
 * The "type" is one of X_AXIS, Y_AXIS, Z_AXIS, T_AXIS.
 */
AxisType File::dimension_type(const std::string &name,
                              units::System::Ptr unit_system) const {
  try {
    if (not find_variable(name)) {
      throw RuntimeError(PISM_ERROR_LOCATION, "coordinate variable " + name + " is missing");
    }

    std::string
      axis          = read_text_attribute(name, "axis"),
      standard_name = read_text_attribute(name, "standard_name"),
      units         = read_text_attribute(name, "units");

    // check if it has units compatible with "seconds":

    units::Unit seconds(unit_system, "seconds");
    if (units::are_convertible(units::Unit(unit_system, units), seconds)) {
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
    if (axis == "T" or axis == "t") {
      return T_AXIS;
    } else if (axis == "X" or axis == "x") {
      return X_AXIS;
    } else if (axis == "Y" or axis == "y") {
      return Y_AXIS;
    } else if (axis == "Z" or axis == "z") {
      return Z_AXIS;
    }

    // check the variable name:
    if (name == "x" or name == "X" or
        name.find("x") == 0 or name.find("X") == 0) {
      return X_AXIS;
    } else if (name == "y" or name == "Y" or
               name.find("y") == 0 or name.find("Y") == 0) {
      return Y_AXIS;
    } else if (name == "z" or name == "Z" or
               name.find("z") == 0 or name.find("Z") == 0) {
      return Z_AXIS;
    } else if (name == "t" or name == "T" or name == "time" or
               name.find("t") == 0 or name.find("T") == 0) {
      return T_AXIS;
    }

    // we have no clue:
    return UNKNOWN_AXIS;
  } catch (RuntimeError &e) {
    e.add_context("getting the type of dimension '%s' in '%s'",
                  name.c_str(), filename().c_str());
    throw;
  }
  return UNKNOWN_AXIS;          // LCOV_EXCL_LINE
}

void File::define_dimension(const std::string &name, size_t length) const {
  try {
    m_impl->nc->def_dim(name, length);
  } catch (RuntimeError &e) {
    e.add_context("defining dimension '%s' in '%s'", name.c_str(), filename().c_str());
    throw;
  }
}

//! \brief Define a variable.
void File::define_variable(const std::string &name, IO_Type nctype, const std::vector<std::string> &dims) const {
  try {
    m_impl->nc->def_var(name, nctype, dims);

    // FIXME: I need to write and tune chunk_dimensions that would be called below before we use
    // this.
    //
    /*
    // if it's not a spatial variable, we're done
    if (dims.size() < 2) {
      return;
    }

    std::vector<size_t> dim_lengths;
    for (unsigned int k = 0; k < dims.size(); ++k) {
      dim_lengths.push_back(this->dimension_length(dims[k]));
    }

    std::vector<size_t> chunk_dims = chunk_dimensions(nctype, dim_lengths);

    m_impl->nc->def_var_chunking(name, chunk_dims);
    */

  } catch (RuntimeError &e) {
    e.add_context("defining variable '%s' in '%s'", name.c_str(), filename().c_str());
    throw;
  }
}

//! \brief Get dimension data (a coordinate variable).
std::vector<double>  File::read_dimension(const std::string &name) const {
  try {
    if (not find_variable(name)) {
      throw RuntimeError(PISM_ERROR_LOCATION, "coordinate variable not found");
    }

    unsigned int length = dimension_length(name);

    std::vector<double> result(length);

    read_variable(name, {0}, {length}, result.data());

    return result;
  } catch (RuntimeError &e) {
    e.add_context("reading dimension '%s' from '%s'", name.c_str(), filename().c_str());
    throw;
  }
}

//! \brief Append to the history global attribute.
/*!
 * Use write_attribute("PISM_GLOBAL", "history", ...) to overwrite "history".
 */
void File::append_history(const std::string &history) const {
  try {
    std::string old_history = read_text_attribute("PISM_GLOBAL", "history");
    redef();
    write_attribute("PISM_GLOBAL", "history", history + old_history);
  } catch (RuntimeError &e) {
    e.add_context("appending to the history attribute in \"" + filename() + "\"");
    throw;
  }
}

//! \brief Write a multiple-valued double attribute.
void File::write_attribute(const std::string &var_name, const std::string &att_name, IO_Type nctype,
                           const std::vector<double> &values) const {
  try {
    redef();
    m_impl->nc->put_att_double(var_name, att_name, nctype, values);
  } catch (RuntimeError &e) {
    e.add_context("writing double attribute '%s:%s' in '%s'",
                  var_name.c_str(), att_name.c_str(), filename().c_str());
    throw;
  }
}

//! \brief Write a text attribute.
void File::write_attribute(const std::string &var_name, const std::string &att_name,
                           const std::string &value) const {
  try {
    redef();
    // ensure that the string is null-terminated
    m_impl->nc->put_att_text(var_name, att_name, value + "\0");
  } catch (RuntimeError &e) {
    e.add_context("writing text attribute '%s:%s' in '%s'",
                  var_name.c_str(), att_name.c_str(), filename().c_str());
    throw;
  }
}

//! \brief Get a double attribute.
std::vector<double> File::read_double_attribute(const std::string &var_name, const std::string &att_name) const {
  try {
    auto att_type = attribute_type(var_name, att_name);

    // Give an understandable error message if a string attribute was found when
    // a number (or a list of numbers) was expected. (We've seen datasets with
    // "valid_min" stored as a string...)
    if (att_type == PISM_CHAR) {
      std::string tmp = read_text_attribute(var_name, att_name);

      throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                    "attribute %s is a string '%s'; expected a number or a list of numbers",
                                    att_name.c_str(), tmp.c_str());
    } else {
      // In this case att_type might be PISM_NAT (if an attribute does not
      // exist), but read_double_attribute can handle that.
      std::vector<double> result;
      m_impl->nc->get_att_double(var_name, att_name, result);
      return result;
    }
  } catch (RuntimeError &e) {
    e.add_context("reading double attribute '%s:%s' from '%s'",
                  var_name.c_str(), att_name.c_str(), filename().c_str());
    throw;
  }
}

//! \brief Get a text attribute.
std::string File::read_text_attribute(const std::string &var_name, const std::string &att_name) const {
  try {
    auto att_type = attribute_type(var_name, att_name);
    if (att_type != PISM_NAT and att_type != PISM_CHAR) {
      // attribute exists and is not a string
      throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                    "attribute %s is not a string", att_name.c_str());
    }

    std::string result;
    m_impl->nc->get_att_text(var_name, att_name, result);
    return result;
  } catch (RuntimeError &e) {
    e.add_context("reading text attribute '%s:%s' from %s", var_name.c_str(), att_name.c_str(), filename().c_str());
    throw;
  }
}

unsigned int File::nattributes(const std::string &var_name) const {
  try {
    int result = 0;
    m_impl->nc->inq_varnatts(var_name, result);
    return result;
  } catch (RuntimeError &e) {
    e.add_context("getting the number of attributes of variable '%s' in '%s'", var_name.c_str(), filename().c_str());
    throw;
  }
}


std::string File::attribute_name(const std::string &var_name, unsigned int n) const {
  try {
    std::string result;
    m_impl->nc->inq_attname(var_name, n, result);
    return result;
  } catch (RuntimeError &e) {
    e.add_context("getting the name of an attribute of variable '%s' in '%s'", var_name.c_str(), filename().c_str());
    throw;
  }
}


IO_Type File::attribute_type(const std::string &var_name, const std::string &att_name) const {
  try {
    IO_Type result;
    m_impl->nc->inq_atttype(var_name, att_name, result);
    return result;
  } catch (RuntimeError &e) {
    e.add_context("getting the type of an attribute of variable '%s' in '%s'", var_name.c_str(), filename().c_str());
    throw;
  }
}


void File::read_variable(const std::string &variable_name,
                           const std::vector<unsigned int> &start,
                           const std::vector<unsigned int> &count,
                          double *ip) const {
  try {
    m_impl->nc->get_vara_double(variable_name, start, count, ip);
  } catch (RuntimeError &e) {
    e.add_context("reading variable '%s' from '%s'", variable_name.c_str(), filename().c_str());
    throw;
  }
}


void File::write_variable(const std::string &variable_name,
                          const std::vector<unsigned int> &start,
                          const std::vector<unsigned int> &count,
                          const double *op) const {
  try {
    m_impl->nc->put_vara_double(variable_name, start, count, op);
  } catch (RuntimeError &e) {
    e.add_context("writing variable '%s' to '%s'", variable_name.c_str(), filename().c_str());
    throw;
  }
}


void File::write_distributed_array(const std::string &variable_name,
                                   const IceGrid &grid,
                                   unsigned int z_count,
                                   const double *input) const {
  try {
    unsigned int t_length = nrecords();
    assert(t_length > 0);

    m_impl->nc->write_darray(variable_name, grid, z_count, t_length - 1, input);
  } catch (RuntimeError &e) {
    e.add_context("writing distributed array '%s' to '%s'",
                  variable_name.c_str(), filename().c_str());
    throw;
  }
}


void File::read_variable_transposed(const std::string &variable_name,
                                    const std::vector<unsigned int> &start,
                                    const std::vector<unsigned int> &count,
                                    const std::vector<unsigned int> &imap, double *ip) const {
  try {
    m_impl->nc->get_varm_double(variable_name, start, count, imap, ip);
  } catch (RuntimeError &e) {
    e.add_context("reading variable '%s' from '%s'", variable_name.c_str(), filename().c_str());
    throw;
  }
}

unsigned int File::nvariables() const {
  int n_vars = 0;

  try {
    m_impl->nc->inq_nvars(n_vars);
  } catch (RuntimeError &e) {
    e.add_context("getting the number of variables in '%s'", filename().c_str());
    throw;
  }

  return n_vars;
}

std::string File::variable_name(unsigned int id) const {
  std::string result;
  try {
    m_impl->nc->inq_varname(id, result);
  } catch (RuntimeError &e) {
    e.add_context("getting the name of %d-th variable in '%s'", id, filename().c_str());
    throw;
  }

  return result;
}

void File::new_grid(int lengthx, int lengthy) const {
  try {
    m_impl->nc->create_grid(lengthx, lengthy);
  } catch (RuntimeError &e) {
    e.add_context("setting a new grid in '%s'", filename().c_str());
    throw;
  }
}

void File::new_timestep(int tsID) const {
  try {
    m_impl->nc->define_timestep(tsID);
  } catch (RuntimeError &e) {
    e.add_context("setting a new timestep in '%s'", filename().c_str());
    throw;
  }
}

void File::reference_date(double time) const {
  try {
    m_impl->nc->def_ref_date(time);
  } catch (RuntimeError &e) {
    e.add_context("setting reference date in '%s'", filename().c_str());
    throw;
  }
}

std::map<std::string, int> File::get_variables_map() const {
  return m_impl->nc->get_var_map();
}

void File::define_vlist() const {
    try {
    m_impl->nc->def_vlist();
  } catch (RuntimeError &e) {
    e.add_context("defining vlist in '%s'", filename().c_str());
    throw;
  }

}

void File::send_diagnostics(const std::set<std::string> &variables) const {
  try {
    m_impl->nc->set_diagvars(variables);
  } catch (RuntimeError &e) {
    e.add_context("setting diagvars in '%s'", filename().c_str());
    throw;
  }

}

void File::set_beforediag(bool value) const {
  m_impl->nc->set_bdiag(value);
}

} // end of namespace pism
