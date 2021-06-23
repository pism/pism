// Copyright (C) 2012, 2013, 2014, 2015, 2016, 2017, 2019, 2020, 2021 PISM Authors
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

#include "NC4File.hh"

// The following is a stupid kludge necessary to make NetCDF 4.x work in
// serial mode in an MPI program:
#ifndef MPI_INCLUDED
#define MPI_INCLUDED 1
#endif
#include <netcdf.h>

#include "pism_type_conversion.hh"
#include "pism/util/pism_utilities.hh"
#include "pism/util/error_handling.hh"

namespace pism {
namespace io {

//! \brief Prints an error message; for debugging.
static void check(const ErrorLocation &where, int return_code) {
  if (return_code != NC_NOERR) {
    throw RuntimeError(where, nc_strerror(return_code));
  }
}

NC4File::NC4File(MPI_Comm c, unsigned int compression_level)
  : NCFile(c), m_compression_level(compression_level) {
  // empty
}

NC4File::~NC4File() {
  // empty
}

// open/create/close

void NC4File::sync_impl() const {

  int stat = nc_sync(m_file_id); check(PISM_ERROR_LOCATION, stat);
}

void NC4File::close_impl() {
  int stat = nc_close(m_file_id); check(PISM_ERROR_LOCATION, stat);

  m_file_id = -1;
}

// redef/enddef
void NC4File::enddef_impl() const {
  int stat = nc_enddef(m_file_id); check(PISM_ERROR_LOCATION, stat);
}

void NC4File::redef_impl() const {
  int stat = nc_redef(m_file_id); check(PISM_ERROR_LOCATION, stat);
}

// dim
void NC4File::def_dim_impl(const std::string &name, size_t length, AxisType dim) const {
  (void) dim;

  int dimid = 0;

  int stat = nc_def_dim(m_file_id, name.c_str(), length, &dimid);
  check(PISM_ERROR_LOCATION, stat);
}

void NC4File::inq_dimid_impl(const std::string &dimension_name, bool &exists) const {
  int tmp = 0;

  int stat = nc_inq_dimid(m_file_id, dimension_name.c_str(), &tmp);

  if (stat == NC_NOERR) {
    exists = true;
  } else {
    exists = false;
  }
}

void NC4File::inq_dimlen_impl(const std::string &dimension_name, unsigned int &result) const {
  int dimid = -1;
  size_t len;

  int stat = nc_inq_dimid(m_file_id, dimension_name.c_str(), &dimid);
  check(PISM_ERROR_LOCATION, stat);

  stat = nc_inq_dimlen(m_file_id, dimid, &len); check(PISM_ERROR_LOCATION, stat);

  result = static_cast<unsigned int>(len);
}

void NC4File::inq_unlimdim_impl(std::string &result) const {
  int dimid = -1;
  std::vector<char> dimname(NC_MAX_NAME + 1, 0);

  int stat = nc_inq_unlimdim(m_file_id, &dimid); check(PISM_ERROR_LOCATION, stat);

  if (dimid == -1) {
    result.clear();
  } else {
    stat = nc_inq_dimname(m_file_id, dimid, dimname.data()); check(PISM_ERROR_LOCATION, stat);

    result = dimname.data();
  }
}

// var
void NC4File::def_var_impl(const std::string &name,
                           IO_Type nctype,
                           const std::vector<std::string> &dims) const {
  std::vector<int> dimids;
  int stat = 0, varid = -1;

  for (auto d : dims) {
    int dimid = -1;
    stat = nc_inq_dimid(m_file_id, d.c_str(), &dimid); check(PISM_ERROR_LOCATION, stat);
    dimids.push_back(dimid);
  }

  stat = nc_def_var(m_file_id, name.c_str(), pism_type_to_nc_type(nctype),
                    static_cast<int>(dims.size()), &dimids[0], &varid);
  check(PISM_ERROR_LOCATION, stat);

  // Compress 2D and 3D variables
  if (m_compression_level > 0 and dims.size() > 1) {
    stat = nc_inq_varid(m_file_id, name.c_str(), &varid); check(PISM_ERROR_LOCATION, stat);
    stat = nc_def_var_deflate(m_file_id, varid, 0, 1, m_compression_level);

    // The NetCDF version used by PISM may not support compression.
    if (stat == NC_EINVAL) {
      stat = NC_NOERR;
    }

    check(PISM_ERROR_LOCATION, stat);
  }
}

void NC4File::def_var_chunking_impl(const std::string &name,
                                   std::vector<size_t> &dimensions) const {
  int stat = 0, varid = 0;

  stat = nc_inq_varid(m_file_id, name.c_str(), &varid);
  check(PISM_ERROR_LOCATION, stat);

  stat = nc_def_var_chunking(m_file_id, varid, NC_CHUNKED, &dimensions[0]);
  check(PISM_ERROR_LOCATION, stat);
}

void NC4File::get_varm_double_impl(const std::string &variable_name,
                                  const std::vector<unsigned int> &start,
                                  const std::vector<unsigned int> &count,
                                  const std::vector<unsigned int> &imap, double *op) const {
  return this->get_put_var_double(variable_name,
                                  start, count, imap, op,
                                  true /*get*/,
                                  true /*transposed*/);
}

void NC4File::get_vara_double_impl(const std::string &variable_name,
                                  const std::vector<unsigned int> &start,
                                  const std::vector<unsigned int> &count,
                                  double *op) const {
  std::vector<unsigned int> dummy;
  return this->get_put_var_double(variable_name,
                                  start, count, dummy, op,
                                  true /*get*/,
                                  false /*not transposed*/);
}

void NC4File::put_vara_double_impl(const std::string &variable_name,
                                  const std::vector<unsigned int> &start,
                                  const std::vector<unsigned int> &count,
                                  const double *op) const {
  std::vector<unsigned int> dummy;
  return this->get_put_var_double(variable_name,
                                  start, count, dummy, const_cast<double*>(op),
                                  false /*put*/,
                                  false /*not transposed*/);
}


void NC4File::inq_nvars_impl(int &result) const {
  int stat = nc_inq_nvars(m_file_id, &result); check(PISM_ERROR_LOCATION, stat);
}

void NC4File::inq_vardimid_impl(const std::string &variable_name, std::vector<std::string> &result) const {
  int ndims = 0, varid = -1;
  std::vector<int> dimids;

  int stat = nc_inq_varid(m_file_id, variable_name.c_str(), &varid); check(PISM_ERROR_LOCATION, stat);

  stat = nc_inq_varndims(m_file_id, varid, &ndims); check(PISM_ERROR_LOCATION, stat);

  if (ndims == 0) {
    result.clear();
    return;
  }

  result.resize(ndims);
  dimids.resize(ndims);

  stat = nc_inq_vardimid(m_file_id, varid, &dimids[0]); check(PISM_ERROR_LOCATION, stat);

  for (int k = 0; k < ndims; ++k) {
    std::vector<char> name(NC_MAX_NAME + 1, 0);

    stat = nc_inq_dimname(m_file_id, dimids[k], name.data()); check(PISM_ERROR_LOCATION, stat);

    result[k] = name.data();
  }
}

int NC4File::get_varid(const std::string &variable_name) const {
  if (variable_name == "PISM_GLOBAL") {
    return NC_GLOBAL;
  }

  int result = 0;
  int stat = nc_inq_varid(m_file_id, variable_name.c_str(), &result);
  check(PISM_ERROR_LOCATION, stat);

  return result;
}

void NC4File::inq_varnatts_impl(const std::string &variable_name, int &result) const {
  int stat = nc_inq_varnatts(m_file_id, get_varid(variable_name), &result);
  check(PISM_ERROR_LOCATION, stat);
}

void NC4File::inq_varid_impl(const std::string &variable_name, bool &exists) const {
  int varid = -1;

  int stat = nc_inq_varid(m_file_id, variable_name.c_str(), &varid);

  exists = (stat == NC_NOERR);
}

void NC4File::inq_varname_impl(unsigned int j, std::string &result) const {
  std::vector<char> varname(NC_MAX_NAME + 1, 0);

  int stat = nc_inq_varname(m_file_id, j, varname.data());
  check(PISM_ERROR_LOCATION, stat);

  result = varname.data();
}

// att

void NC4File::get_att_double_impl(const std::string &variable_name,
                                  const std::string &att_name,
                                  std::vector<double> &result) const {
  int varid = get_varid(variable_name);

  // Read the attribute length:
  size_t attlen = 0;
  int stat = nc_inq_attlen(m_file_id, varid, att_name.c_str(), &attlen);

  int len = 0;
  if (stat == NC_NOERR) {
    len = static_cast<int>(attlen);
  } else if (stat == NC_ENOTATT) {
    len = 0;
  } else {
    check(PISM_ERROR_LOCATION, stat);
    len = 0;
  }

  if (len == 0) {
    result.clear();
    return;
  }

  result.resize(len);

  // Now read data and see if we succeeded:
  stat = nc_get_att_double(m_file_id, varid, att_name.c_str(), &result[0]);
  check(PISM_ERROR_LOCATION, stat);
}

// Get a text (character array) attribute on rank 0.
static void get_att_text(int ncid, int varid, const std::string &att_name,
                         std::string &result) {
  size_t attlen = 0;
  int stat = nc_inq_attlen(ncid, varid, att_name.c_str(), &attlen);
  if (stat != NC_NOERR) {
    result = "";
    return;
  }

  std::vector<char> buffer(attlen + 1, 0);
  stat = nc_get_att_text(ncid, varid, att_name.c_str(), buffer.data());

  result = (stat == NC_NOERR) ? buffer.data() : "";
}

// Get a string attribute on rank 0. In "string array" attributes array elements are
// concatenated using "," as the separator.
static void get_att_string(int ncid, int varid, const std::string &att_name,
                           std::string &result) {
  size_t attlen = 0;
  int stat = nc_inq_attlen(ncid, varid, att_name.c_str(), &attlen);
  if (stat != NC_NOERR) {
    result = "";
    return;
  }

  std::vector<char*> buffer(attlen + 1, 0);
  stat = nc_get_att_string(ncid, varid, att_name.c_str(), buffer.data());
  if (stat == NC_NOERR) {
    std::vector<std::string> strings(attlen);
    for (size_t k = 0; k < attlen; ++k) {
      strings[k] = buffer[k];
    }
    result = join(strings, ",");
  } else {
    result = "";
  }
  stat = nc_free_string(attlen, buffer.data());
  check(PISM_ERROR_LOCATION, stat);
}

void NC4File::get_att_text_impl(const std::string &variable_name,
                               const std::string &att_name, std::string &result) const {
  int varid = get_varid(variable_name);

  nc_type nctype;
  int stat = nc_inq_atttype(m_file_id, varid, att_name.c_str(), &nctype);

  if (stat == NC_NOERR) {
    if (nctype == NC_CHAR) {
      pism::io::get_att_text(m_file_id, varid, att_name, result);
    } else if (nctype == NC_STRING) {
      pism::io::get_att_string(m_file_id, varid, att_name, result);
    } else {
      result = "";
    }
  } else if (stat == NC_ENOTATT) {
    result = "";
  } else {
    check(PISM_ERROR_LOCATION, stat);
  }
}

void NC4File::del_att_impl(const std::string &variable_name, const std::string &att_name) const {
  int stat = nc_del_att(m_file_id, get_varid(variable_name), att_name.c_str());
  check(PISM_ERROR_LOCATION, stat);
}

void NC4File::put_att_double_impl(const std::string &variable_name,
                                  const std::string &att_name,
                                  IO_Type xtype,
                                  const std::vector<double> &data) const {
  int stat = nc_put_att_double(m_file_id, get_varid(variable_name), att_name.c_str(),
                               xtype, data.size(), &data[0]);
  check(PISM_ERROR_LOCATION, stat);
}

void NC4File::put_att_text_impl(const std::string &variable_name,
                                const std::string &att_name,
                                const std::string &value) const {
  int stat = nc_put_att_text(m_file_id, get_varid(variable_name), att_name.c_str(),
                             value.size(), value.c_str());
  check(PISM_ERROR_LOCATION, stat);
}

void NC4File::inq_attname_impl(const std::string &variable_name,
                               unsigned int n,
                               std::string &result) const {
  std::vector<char> name(NC_MAX_NAME + 1, 0);
  int stat = nc_inq_attname(m_file_id, get_varid(variable_name), n, name.data());
  check(PISM_ERROR_LOCATION, stat);

  result = name.data();
}

void NC4File::inq_atttype_impl(const std::string &variable_name,
                               const std::string &att_name,
                               IO_Type &result) const {
  nc_type tmp = NC_NAT;
  int stat = nc_inq_atttype(m_file_id, get_varid(variable_name), att_name.c_str(), &tmp);

  if (stat == NC_ENOTATT) {
    tmp = NC_NAT;
  } else {
    check(PISM_ERROR_LOCATION, stat);
  }

  result = nc_type_to_pism_type(tmp);
}

// misc

void NC4File::set_fill_impl(int fillmode, int &old_modep) const {
  int stat = nc_set_fill(m_file_id, fillmode, &old_modep); check(PISM_ERROR_LOCATION, stat);
}

void NC4File::set_access_mode(int, bool) const {
  // empty
}

void NC4File::get_put_var_double(const std::string &variable_name,
                                const std::vector<unsigned int> &start,
                                const std::vector<unsigned int> &count,
                                const std::vector<unsigned int> &imap_input,
                                double *op,
                                bool get,
                                bool transposed) const {
  int stat, varid, ndims = static_cast<int>(start.size());
  std::vector<unsigned int> imap = imap_input;

  if (not transposed) {
    imap.resize(ndims);
  }

  std::vector<size_t> nc_start(ndims), nc_count(ndims);
  std::vector<ptrdiff_t> nc_imap(ndims), nc_stride(ndims);

  stat = nc_inq_varid(m_file_id, variable_name.c_str(), &varid); check(PISM_ERROR_LOCATION, stat);

  for (int j = 0; j < ndims; ++j) {
    nc_start[j]  = start[j];
    nc_count[j]  = count[j];
    nc_imap[j]   = imap[j];
    nc_stride[j] = 1;
  }

  if (transposed) {

    set_access_mode(varid, transposed);

    if (get) {
      stat = nc_get_varm_double(m_file_id, varid,
                                &nc_start[0], &nc_count[0], &nc_stride[0], &nc_imap[0],
                                op); check(PISM_ERROR_LOCATION, stat);
    } else {
      stat = nc_put_varm_double(m_file_id, varid,
                                &nc_start[0], &nc_count[0], &nc_stride[0], &nc_imap[0],
                                op); check(PISM_ERROR_LOCATION, stat);
    }
  } else {

    set_access_mode(varid, transposed);

    if (get) {
      stat = nc_get_vara_double(m_file_id, varid,
                                &nc_start[0], &nc_count[0],
                                op); check(PISM_ERROR_LOCATION, stat);
    } else {
      stat = nc_put_vara_double(m_file_id, varid,
                                &nc_start[0], &nc_count[0],
                                op); check(PISM_ERROR_LOCATION, stat);
    }
  }
}

} // end of namespace io
} // end of namespace pism
