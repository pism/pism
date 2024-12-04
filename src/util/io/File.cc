// Copyright (C) 2012--2024 PISM Authors
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

#include <petscvec.h>
#include <set>

#include "pism/util/io/File.hh"
#include "pism/util/Grid.hh"
#include "pism/util/io/NC_Serial.hh"
#include "pism/util/io/NC4_Serial.hh"

#include "pism/pism_config.hh"

#if (Pism_USE_PARALLEL_NETCDF4==1)
#include "pism/util/io/NC4_Par.hh"
#endif

#if (Pism_USE_PNETCDF==1)
#include "pism/util/io/PNCFile.hh"
#endif

#include "pism/util/error_handling.hh"
#include "pism/util/io/io_helpers.hh"
#include "pism/util/io/IO_Flags.hh"
#include "pism/util/pism_utilities.hh"

namespace pism {

struct File::Impl {
  MPI_Comm com;
  std::shared_ptr<io::NCFile> nc;

  std::set<std::string> written_variables;
};

io::Backend string_to_backend(const std::string &backend) {
  std::map<std::string, io::Backend> backends =
    {
     {"netcdf3", io::PISM_NETCDF3},
     {"netcdf4_parallel", io::PISM_NETCDF4_PARALLEL},
     {"netcdf4_serial", io::PISM_NETCDF4_SERIAL},
     {"pnetcdf", io::PISM_PNETCDF},
  };

  if (backends.find(backend) != backends.end()) {
    return backends[backend];
  }

  throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                "unknown or unsupported I/O backend: %s",
                                backend.c_str());
}

static std::string backend_to_string(io::Backend backend) {
  std::map<io::Backend, std::string> backends =
    {
     {io::PISM_GUESS, "unknown"},
     {io::PISM_NETCDF3, "netcdf3"},
     {io::PISM_NETCDF4_PARALLEL, "netcdf4_parallel"},
     {io::PISM_NETCDF4_SERIAL, "netcdf4_serial"},
     {io::PISM_PNETCDF, "pnetcdf"}
  };

  return backends[backend];
}

// Chooses the best available I/O backend for reading from 'filename'.
static io::Backend choose_backend(MPI_Comm com, const std::string &filename) {

  std::string format;
  {
    // This is the rank-0-only purely-serial mode of accessing NetCDF files, but it
    // supports all the kinds of NetCDF, so this is fine.
    io::NC_Serial file(com);

    file.open(filename, io::PISM_READONLY);
    format = file.get_format();
    file.close();
  }

#if (Pism_USE_PARALLEL_NETCDF4==1)
  if (format == "netcdf4") {
    return io::PISM_NETCDF4_PARALLEL;
  }
#endif

#if (Pism_USE_PNETCDF==1)
  if (format != "netcdf4") {
    return io::PISM_PNETCDF;
  }
#endif

  // this choice is appropriate for both NetCDF-3 and NetCDF-4
  return io::PISM_NETCDF3;
}

static std::shared_ptr<io::NCFile> create_backend(MPI_Comm com, io::Backend backend) {
  int size = 1;
  MPI_Comm_size(com, &size);

  switch (backend) {

  case io::PISM_NETCDF3:
    return std::make_shared<io::NC_Serial>(com);

  case io::PISM_NETCDF4_SERIAL:
    return std::make_shared<io::NC4_Serial>(com);

  case io::PISM_NETCDF4_PARALLEL:
#if (Pism_USE_PARALLEL_NETCDF4 == 1)
    return std::make_shared<io::NC4_Par>(com);
#else
    break;
#endif

  case io::PISM_PNETCDF:
#if (Pism_USE_PNETCDF == 1)
    return std::make_shared<io::PNCFile>(com);
#else
    break;
#endif

  case io::PISM_GUESS:
    break;
  } // end of switch (backend)

  auto backend_name = backend_to_string(backend);

  throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                "unknown or unsupported I/O backend: %s",
                                backend_name.c_str());
}

File::File(MPI_Comm com, const std::string &filename, io::Backend backend, io::Mode mode)
  : m_impl(new Impl) {

  if (filename.empty()) {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                  "cannot open file: provided file name is empty");
  }

  if (backend == io::PISM_GUESS) {
    backend = choose_backend(com, filename);
  }

  m_impl->com = com;
  m_impl->nc  = create_backend(m_impl->com, backend);

  this->open(filename, mode);
}

File::~File() {
  if (m_impl->nc and not name().empty()) {
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

MPI_Comm File::com() const {
  return m_impl->com;
}

void File::set_compression_level(int level) const {
  m_impl->nc->set_compression_level(level);
}

void File::open(const std::string &filename, io::Mode mode) {
  try {

    // opening for reading
    if (mode == io::PISM_READONLY) {

      m_impl->nc->open(filename, mode);

    } else if (mode == io::PISM_READWRITE_CLOBBER or mode == io::PISM_READWRITE_MOVE) {

      if (mode == io::PISM_READWRITE_MOVE) {
        io::move_if_exists(m_impl->com, filename);
      } else {
        io::remove_if_exists(m_impl->com, filename);
      }

      m_impl->nc->create(filename);

      int old_fill;
      m_impl->nc->set_fill(io::PISM_NOFILL, old_fill);
    } else if (mode == io::PISM_READWRITE) {                      // mode == io::PISM_READWRITE

      m_impl->nc->open(filename, mode);

      int old_fill;
      m_impl->nc->set_fill(io::PISM_NOFILL, old_fill);
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
    e.add_context("closing \"" + name() + "\"");
    throw;
  }
}

void File::sync() const {
  try {
    m_impl->nc->sync();
  } catch (RuntimeError &e) {
    e.add_context("synchronizing \"" + name() + "\"");
    throw;
  }
}

void File::redef() const {
  try {
    m_impl->nc->redef();
  } catch (RuntimeError &e) {
    e.add_context("switching to define mode; file \"" + name() + "\"");
    throw;
  }
}


void File::enddef() const {
  try {
    m_impl->nc->enddef();
  } catch (RuntimeError &e) {
    e.add_context("switching to data mode; file \"" + name() + "\"");
    throw;
  }
}

std::string File::name() const {
  return m_impl->nc->filename();
}


//! \brief Get the number of records. Uses the length of an unlimited dimension.
unsigned int File::nrecords() const {
  try {
    std::string dim;
    m_impl->nc->inq_unlimdim(dim);

    if (dim.empty()) {
      return 1;                 // one record
    }

    return this->dimension_length(dim);
  } catch (RuntimeError &e) {
    e.add_context("getting the number of records in file \"" + name() + "\"");
    throw;
  }
  return 0;                     // LCOV_EXCL_LINE
}

//! \brief Get the number of records of a certain variable. Uses the length of
//! an associated "time" dimension.
unsigned int File::nrecords(const std::string &variable_name, const std::string &std_name,
                            units::System::Ptr unit_system) const {
  try {
    auto var = find_variable(variable_name, std_name);

    if (not var.exists) {
      return 0;
    }

    for (const auto &d : dimensions(var.name)) {
      if (dimension_type(d, unit_system) == T_AXIS) {
        return this->dimension_length(d);
      }
    }

    return 1;                   // one record
  } catch (RuntimeError &e) {
    e.add_context("getting the number of records of variable '%s' ('%s') in '%s'",
                  variable_name.c_str(), std_name.c_str(), name().c_str());
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
        std::string var_name      = variable_name(j);
        std::string attribute = read_text_attribute(var_name, "standard_name");

        if (attribute.empty()) {
          continue;
        }

        if (attribute == std_name) {
          if (not result.exists) {
            result.exists = true;
            result.name = var_name;
          } else {
            throw RuntimeError::formatted(PISM_ERROR_LOCATION, "inconsistency in '%s': variables '%s' and '%s'\n"
                                          "have the same standard_name (%s)",
                                          name().c_str(), result.name.c_str(),
                                          var_name.c_str(), attribute.c_str());
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
    }

  } catch (RuntimeError &e) {
    e.add_context("searching for variable '%s' ('%s') in '%s'", short_name.c_str(), std_name.c_str(), name().c_str());
    throw;
  }

  return result;
}

//! \brief Checks if a variable exists.
bool File::variable_exists(const std::string &variable_name) const {
  try {
    bool exists = false;
    m_impl->nc->inq_varid(variable_name, exists);
    return exists;
  } catch (RuntimeError &e) {
    e.add_context("searching for variable '%s' in '%s'", variable_name.c_str(),
                  name().c_str());
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
                  name().c_str());
    throw;
  }
}


//! \brief Checks if a dimension exists.
bool File::dimension_exists(const std::string &dimension_name) const {
  try {
    bool exists = false;
    m_impl->nc->inq_dimid(dimension_name, exists);
    return exists;
  } catch (RuntimeError &e) {
    e.add_context("searching for dimension '%s' in '%s'", dimension_name.c_str(),
                  name().c_str());
    throw;
  }
}

//! \brief Get the length of a dimension.
/*!
 * Sets result to 0 if a dimension does not exist.
 */
unsigned int File::dimension_length(const std::string &dimension_name) const {
  try {
    if (dimension_exists(dimension_name)) {
      unsigned int result = 0;
      m_impl->nc->inq_dimlen(dimension_name, result);
      return result;
    }

    return 0;
  } catch (RuntimeError &e) {
    e.add_context("getting the length of dimension '%s' in '%s'", dimension_name.c_str(),
                  name().c_str());
    throw;
  }
}

AxisType axis_type_from_string(const std::string &input) {
  if (input == "T" or input == "t") {
    return T_AXIS;
  }

  if (input == "X" or input == "x") {
    return X_AXIS;
  }

  if (input == "Y" or input == "y") {
    return Y_AXIS;
  }

  if (input == "Z" or input == "z") {
    return Z_AXIS;
  }

  return UNKNOWN_AXIS;
}

//! \brief Get the "type" of a dimension.
/*!
 * The "type" is one of X_AXIS, Y_AXIS, Z_AXIS, T_AXIS.
 */
AxisType File::dimension_type(const std::string &dimension_name,
                              units::System::Ptr unit_system) const {
  try {
    if (not variable_exists(dimension_name)) {
      throw RuntimeError(PISM_ERROR_LOCATION, "coordinate variable " + dimension_name + " is missing");
    }

    std::string
      axis          = read_text_attribute(dimension_name, "axis"),
      standard_name = read_text_attribute(dimension_name, "standard_name"),
      units         = read_text_attribute(dimension_name, "units");

    // check if it has units compatible with "seconds":

    units::Unit seconds(unit_system, "seconds");
    if (units::are_convertible(units::Unit(unit_system, units), seconds)) {
      return T_AXIS;
    }

    // check the standard_name attribute:
    if (standard_name == "time") {
      return T_AXIS;
    }

    if (standard_name == "projection_x_coordinate" or
        standard_name == "grid_longitude") {
      return X_AXIS;
    }

    if (standard_name == "projection_y_coordinate" or
        standard_name == "grid_latitude") {
      return Y_AXIS;
    }

    {
      AxisType tmp = axis_type_from_string(axis);
      if (tmp != UNKNOWN_AXIS) {
        return tmp;
      }
    }

    // check the variable name:
    if (member(dimension_name, {"x", "X", "rlon"}) or
        dimension_name.find('x') == 0 or dimension_name.find('X') == 0) {
      return X_AXIS;
    }

    if (member(dimension_name, {"y", "Y", "rlat"}) or
        dimension_name.find('y') == 0 or dimension_name.find('Y') == 0) {
      return Y_AXIS;
    }

    if (dimension_name == "z" or dimension_name == "Z" or
        dimension_name.find('z') == 0 or dimension_name.find('Z') == 0) {
      return Z_AXIS;
    }

    if (dimension_name == "t" or dimension_name == "T" or dimension_name == "time" or
        dimension_name.find('t') == 0 or dimension_name.find('T') == 0) {
      return T_AXIS;
    }

    // we have no clue:
    return UNKNOWN_AXIS;
  } catch (RuntimeError &e) {
    e.add_context("getting the type of dimension '%s' in '%s'",
                  dimension_name.c_str(), name().c_str());
    throw;
  }
  return UNKNOWN_AXIS;          // LCOV_EXCL_LINE
}

void File::define_dimension(const std::string &dimension_name, size_t length) const {
  try {
    m_impl->nc->def_dim(dimension_name, length);
  } catch (RuntimeError &e) {
    e.add_context("defining dimension '%s' in '%s'", dimension_name.c_str(),
                  name().c_str());
    throw;
  }
}

//! \brief Define a variable.
void File::define_variable(const std::string &variable_name, io::Type nctype,
                           const std::vector<std::string> &dims) const {
  try {
    m_impl->nc->def_var(variable_name, nctype, dims);

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
    e.add_context("defining variable '%s' in '%s'", variable_name.c_str(),
                  name().c_str());
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
    e.add_context("appending to the history attribute in \"" + name() + "\"");
    throw;
  }
}

//! \brief Write a multiple-valued double attribute.
void File::write_attribute(const std::string &var_name, const std::string &att_name, io::Type nctype,
                           const std::vector<double> &values) const {
  try {
    redef();
    m_impl->nc->put_att_double(var_name, att_name, nctype, values);
  } catch (RuntimeError &e) {
    e.add_context("writing double attribute '%s:%s' in '%s'",
                  var_name.c_str(), att_name.c_str(), name().c_str());
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
                  var_name.c_str(), att_name.c_str(), name().c_str());
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
    if (att_type == io::PISM_CHAR) {
      std::string tmp = read_text_attribute(var_name, att_name);

      throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                    "attribute %s is a string '%s'; expected a number or a list of numbers",
                                    att_name.c_str(), tmp.c_str());
    }

    // In this case att_type might be io::PISM_NAT (if an attribute does not
    // exist), but read_double_attribute can handle that.
    std::vector<double> result;
    m_impl->nc->get_att_double(var_name, att_name, result);
    return result;
  } catch (RuntimeError &e) {
    e.add_context("reading double attribute '%s:%s' from '%s'",
                  var_name.c_str(), att_name.c_str(), name().c_str());
    throw;
  }
}

//! \brief Get a text attribute.
std::string File::read_text_attribute(const std::string &var_name, const std::string &att_name) const {
  try {
    auto att_type = attribute_type(var_name, att_name);
    if (att_type != io::PISM_NAT and att_type != io::PISM_CHAR) {
      // attribute exists and is not a string
      throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                    "attribute %s is not a string", att_name.c_str());
    }

    std::string result;
    m_impl->nc->get_att_text(var_name, att_name, result);
    return result;
  } catch (RuntimeError &e) {
    e.add_context("reading text attribute '%s:%s' from %s", var_name.c_str(), att_name.c_str(), name().c_str());
    throw;
  }
}

unsigned int File::nattributes(const std::string &var_name) const {
  try {
    int result = 0;
    m_impl->nc->inq_varnatts(var_name, result);
    return result;
  } catch (RuntimeError &e) {
    e.add_context("getting the number of attributes of variable '%s' in '%s'", var_name.c_str(), name().c_str());
    throw;
  }
}


std::string File::attribute_name(const std::string &var_name, unsigned int n) const {
  try {
    std::string result;
    m_impl->nc->inq_attname(var_name, n, result);
    return result;
  } catch (RuntimeError &e) {
    e.add_context("getting the name of an attribute of variable '%s' in '%s'", var_name.c_str(), name().c_str());
    throw;
  }
}


io::Type File::attribute_type(const std::string &var_name, const std::string &att_name) const {
  try {
    io::Type result;
    m_impl->nc->inq_atttype(var_name, att_name, result);
    return result;
  } catch (RuntimeError &e) {
    e.add_context("getting the type of an attribute of variable '%s' in '%s'", var_name.c_str(), name().c_str());
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
    e.add_context("reading variable '%s' from '%s'", variable_name.c_str(), name().c_str());
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
    e.add_context("writing variable '%s' to '%s'", variable_name.c_str(), name().c_str());
    throw;
  }
}


void File::write_distributed_array(const std::string &variable_name,
                                   const Grid &grid,
                                   unsigned int z_count,
                                   bool time_dependent,
                                   const double *input) const {
  try {
    unsigned int t_length = nrecords();
    assert(t_length > 0);

    m_impl->nc->write_darray(variable_name, grid, z_count, time_dependent, t_length - 1, input);
  } catch (RuntimeError &e) {
    e.add_context("writing distributed array '%s' to '%s'",
                  variable_name.c_str(), name().c_str());
    throw;
  }
}

unsigned int File::nvariables() const {
  int n_vars = 0;

  try {
    m_impl->nc->inq_nvars(n_vars);
  } catch (RuntimeError &e) {
    e.add_context("getting the number of variables in '%s'", name().c_str());
    throw;
  }

  return n_vars;
}

std::string File::variable_name(unsigned int id) const {
  std::string result;
  try {
    m_impl->nc->inq_varname(id, result);
  } catch (RuntimeError &e) {
    e.add_context("getting the name of %d-th variable in '%s'", id, name().c_str());
    throw;
  }

  return result;
}

void File::set_variable_was_written(const std::string &name) const {
  m_impl->written_variables.insert(name);
}

bool File::get_variable_was_written(const std::string &name) const {
  return member(name, m_impl->written_variables);
}

} // end of namespace pism
