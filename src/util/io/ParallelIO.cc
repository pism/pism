/* Copyright (C) 2019 PISM Authors
 *
 * This file is part of PISM.
 *
 * PISM is free software; you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation; either version 3 of the License, or (at your option) any later
 * version.
 *
 * PISM is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License
 * along with PISM; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 */

#include <cassert>

// Why do I need this???
#define _NETCDF
#include <pio.h>

#include "ParallelIO.hh"

#include "pism/util/error_handling.hh"
#include "pism/util/io/pism_type_conversion.hh"


namespace pism {
namespace io {

static void check(const ErrorLocation &where, int return_code) {
  char message[PIO_MAX_NAME + 1];
  if (return_code != PIO_NOERR) {
    PIOc_strerror(return_code, message);
    throw RuntimeError(where, message);
  }
}

ParallelIO::ParallelIO(MPI_Comm com, int iosysid)
  : NCFile(com),
    m_iosysid(iosysid) {
  assert(iosysid != -1);
}

ParallelIO::~ParallelIO() {
  // empty
}

int ParallelIO::open_impl(const std::string &filename, IO_Mode mode) {
  int open_mode = mode == PISM_READONLY ? PIO_NOWRITE : PIO_WRITE;

  int stat = PIOc_open(m_iosysid, filename.c_str(), open_mode, &m_file_id);
  check(PISM_ERROR_LOCATION, stat);

  return stat;
}

int ParallelIO::create_impl(const std::string &filename) {

  // find the best available I/O type
  int iotype = 0;
  for (int t : {PIO_IOTYPE_PNETCDF, PIO_IOTYPE_NETCDF4P, PIO_IOTYPE_NETCDF4C, PIO_IOTYPE_NETCDF}) {
    if (PIOc_iotype_available(t)) {
      iotype = t;
      break;
    }
  }

  int mode = NC_CLOBBER;
  if (iotype == PIO_IOTYPE_PNETCDF) {
    mode |= NC_64BIT_DATA;
  }

  int stat = PIOc_createfile(m_iosysid, &m_file_id, &iotype, filename.c_str(),
                             mode);
  check(PISM_ERROR_LOCATION, stat);

  return stat;
}

int ParallelIO::sync_impl() const {

  int stat = PIOc_sync(m_file_id); check(PISM_ERROR_LOCATION, stat);

  return stat;
}

int ParallelIO::close_impl() {
  int stat = PIOc_closefile(m_file_id); check(PISM_ERROR_LOCATION, stat);
  return 0;
}

// redef/enddef
int ParallelIO::enddef_impl() const {
  int stat = PIOc_enddef(m_file_id); check(PISM_ERROR_LOCATION, stat);
  return 0;
}

int ParallelIO::redef_impl() const {
  int stat = PIOc_redef(m_file_id); check(PISM_ERROR_LOCATION, stat);
  return 0;
}

// dim
int ParallelIO::def_dim_impl(const std::string &name, size_t length) const {
  int dim_id = 0;
  int stat = PIOc_def_dim(m_file_id, name.c_str(), length, &dim_id); check(PISM_ERROR_LOCATION, stat);
  return 0;
}

int ParallelIO::inq_dimid_impl(const std::string &dimension_name, bool &exists) const {
  int tmp, stat;

  stat = PIOc_inq_dimid(m_file_id, dimension_name.c_str(), &tmp);

  if (stat == PIO_NOERR) {
    exists = true;
  } else {
    exists = false;
  }

  return 0;
}

int ParallelIO::inq_dimlen_impl(const std::string &dimension_name, unsigned int &result) const {
  int stat, dimid = -1;
  PIO_Offset len;

  stat = PIOc_inq_dimid(m_file_id, dimension_name.c_str(), &dimid); check(PISM_ERROR_LOCATION, stat);

  stat = PIOc_inq_dimlen(m_file_id, dimid, &len); check(PISM_ERROR_LOCATION, stat);

  result = static_cast<unsigned int>(len);

  return stat;
}

int ParallelIO::inq_unlimdim_impl(std::string &result) const {
  int stat, dimid;
  char dimname[PIO_MAX_NAME + 1];

  stat = PIOc_inq_unlimdim(m_file_id, &dimid); check(PISM_ERROR_LOCATION, stat);

  if (dimid == -1) {
    result.clear();
  } else {
    stat = PIOc_inq_dimname(m_file_id, dimid, dimname); check(PISM_ERROR_LOCATION, stat);

    result = dimname;
  }

  return stat;
}

// var
int ParallelIO::def_var_impl(const std::string &name, IO_Type nctype,
                             const std::vector<std::string> &dims) const {
  std::vector<int> dimids;
  int stat, varid;

  for (auto d : dims) {
    int dimid = -1;
    stat = PIOc_inq_dimid(m_file_id, d.c_str(), &dimid); check(PISM_ERROR_LOCATION, stat);
    dimids.push_back(dimid);
  }

  stat = PIOc_def_var(m_file_id, name.c_str(), pism_type_to_nc_type(nctype),
                      static_cast<int>(dims.size()), dimids.data(), &varid); check(PISM_ERROR_LOCATION, stat);

  return stat;
}

int ParallelIO::def_var_chunking_impl(const std::string &name,
                                      std::vector<size_t> &dimensions) const {
  (void) name;
  (void) dimensions;
  // FIXME
  return 0;
}

int ParallelIO::get_vara_double_impl(const std::string &variable_name,
                                     const std::vector<unsigned int> &start,
                                     const std::vector<unsigned int> &count,
                                     double *input) const {
  int stat, varid, ndims = static_cast<int>(start.size());

  std::vector<PIO_Offset>
    nc_start(ndims),
    nc_count(ndims);

  stat = PIOc_inq_varid(m_file_id, variable_name.c_str(), &varid); check(PISM_ERROR_LOCATION, stat);

  for (int j = 0; j < ndims; ++j) {
    nc_start[j] = start[j];
    nc_count[j] = count[j];
  }

  stat = PIOc_get_vara_double(m_file_id, varid, nc_start.data(), nc_count.data(), input);
  check(PISM_ERROR_LOCATION, stat);

  return stat;
}

int ParallelIO::put_vara_double_impl(const std::string &variable_name,
                                     const std::vector<unsigned int> &start,
                                     const std::vector<unsigned int> &count,
                                     const double *output) const {
  int stat, varid, ndims = static_cast<int>(start.size());

  std::vector<PIO_Offset>
    nc_start(ndims),
    nc_count(ndims);

  stat = PIOc_inq_varid(m_file_id, variable_name.c_str(), &varid); check(PISM_ERROR_LOCATION, stat);

  for (int j = 0; j < ndims; ++j) {
    nc_start[j] = start[j];
    nc_count[j] = count[j];
  }

  stat = PIOc_put_vara_double(m_file_id, varid, nc_start.data(), nc_count.data(), output);
  check(PISM_ERROR_LOCATION, stat);

  return stat;

}

int ParallelIO::get_varm_double_impl(const std::string &variable_name,
                                     const std::vector<unsigned int> &start,
                                     const std::vector<unsigned int> &count,
                                     const std::vector<unsigned int> &imap,
                                     double *input) const {
  (void) variable_name;
  (void) start;
  (void) count;
  (void) imap;
  (void) input;
  throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                "ParallelIO does not support transposed access");
}

int ParallelIO::inq_nvars_impl(int &result) const {
  int stat = PIOc_inq_nvars(m_file_id, &result); check(PISM_ERROR_LOCATION, stat);
  return stat;
}

int ParallelIO::inq_vardimid_impl(const std::string &variable_name, std::vector<std::string> &result) const {
  int stat, ndims, varid = -1;
  std::vector<int> dimids;

  stat = PIOc_inq_varid(m_file_id, variable_name.c_str(), &varid); check(PISM_ERROR_LOCATION, stat);

  stat = PIOc_inq_varndims(m_file_id, varid, &ndims); check(PISM_ERROR_LOCATION, stat);

  if (ndims == 0) {
    result.clear();
    return 0;
  }

  result.resize(ndims);
  dimids.resize(ndims);

  stat = PIOc_inq_vardimid(m_file_id, varid, dimids.data()); check(PISM_ERROR_LOCATION, stat);

  for (int k = 0; k < ndims; ++k) {
    char name[PIO_MAX_NAME];
    memset(name, 0, PIO_MAX_NAME);

    stat = PIOc_inq_dimname(m_file_id, dimids[k], name); check(PISM_ERROR_LOCATION, stat);

    result[k] = name;
  }

  return 0;
}

int ParallelIO::inq_varnatts_impl(const std::string &variable_name, int &result) const {
  int stat, varid = -1;

  if (variable_name == "PISM_GLOBAL") {
    varid = PIO_GLOBAL;
  } else {
    stat = PIOc_inq_varid(m_file_id, variable_name.c_str(), &varid); check(PISM_ERROR_LOCATION, stat);
  }

  stat = PIOc_inq_varnatts(m_file_id, varid, &result); check(PISM_ERROR_LOCATION, stat);

  return 0;
}

int ParallelIO::inq_varid_impl(const std::string &variable_name, bool &exists) const {
  int stat, flag = -1;

  stat = PIOc_inq_varid(m_file_id, variable_name.c_str(), &flag);

  if (stat == NC_NOERR) {
    flag = 1;
  } else {
    flag = 0;
  }

  exists = (flag == 1);

  return 0;
}

int ParallelIO::inq_varname_impl(unsigned int j, std::string &result) const {
  int stat;
  char varname[PIO_MAX_NAME];
  memset(varname, 0, PIO_MAX_NAME);

  stat = PIOc_inq_varname(m_file_id, j, varname); check(PISM_ERROR_LOCATION, stat);

  result = varname;

  return stat;
}

// att
int ParallelIO::get_att_double_impl(const std::string &variable_name,
                                    const std::string &att_name,
                                    std::vector<double> &result) const {
  int stat, len, varid = -1;
  PIO_Offset attlen;

  // Read the attribute length:
  if (variable_name == "PISM_GLOBAL") {
    varid = NC_GLOBAL;
  } else {
    stat = PIOc_inq_varid(m_file_id, variable_name.c_str(), &varid); check(PISM_ERROR_LOCATION, stat);
  }

  stat = PIOc_inq_attlen(m_file_id, varid, att_name.c_str(), &attlen);

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
    return 0;
  }

  result.resize(len);

  stat = PIOc_get_att_double(m_file_id, varid, att_name.c_str(), result.data());
  check(PISM_ERROR_LOCATION, stat);

  return 0;
}

int ParallelIO::get_att_text_impl(const std::string &variable_name, const std::string &att_name, std::string &result) const {
  int stat, len, varid = -1;

  // Read the attribute length:
  PIO_Offset attlen;

  if (variable_name == "PISM_GLOBAL") {
    varid = NC_GLOBAL;
  } else {
    stat = PIOc_inq_varid(m_file_id, variable_name.c_str(), &varid); check(PISM_ERROR_LOCATION, stat);
  }

  stat = PIOc_inq_attlen(m_file_id, varid, att_name.c_str(), &attlen);
  if (stat == NC_NOERR) {
    len = static_cast<int>(attlen);
  } else {
    len = 0;
  }

  // Allocate some memory or clear result.
  if (len == 0) {
    result.clear();
    return 0;
  }

  std::vector<char> str(len + 1, '\0');

  // Now read the string and see if we succeeded:
  stat = PIOc_get_att_text(m_file_id, varid, att_name.c_str(), str.data());
  check(PISM_ERROR_LOCATION, stat);

  result = str.data();

  return 0;

}

int ParallelIO::put_att_double_impl(const std::string &variable_name,
                                    const std::string &att_name,
                                    IO_Type xtype,
                                    const std::vector<double> &data) const {
  int stat = 0;

  int varid = -1;

  if (variable_name == "PISM_GLOBAL") {
    varid = NC_GLOBAL;
  } else {
    stat = PIOc_inq_varid(m_file_id, variable_name.c_str(), &varid);
    check(PISM_ERROR_LOCATION, stat);
  }

  stat = PIOc_put_att_double(m_file_id, varid, att_name.c_str(),
                             pism_type_to_nc_type(xtype), data.size(), data.data());
  check(PISM_ERROR_LOCATION, stat);

  return stat;
}

int ParallelIO::put_att_text_impl(const std::string &variable_name,
                                  const std::string &att_name,
                                  const std::string &value) const {
  int stat = 0, varid = -1;

  if (variable_name == "PISM_GLOBAL") {
    varid = NC_GLOBAL;
  } else {
    stat = PIOc_inq_varid(m_file_id, variable_name.c_str(), &varid);
    check(PISM_ERROR_LOCATION, stat);
  }

  stat = PIOc_put_att_text(m_file_id, varid, att_name.c_str(), value.size(), value.c_str());
  check(PISM_ERROR_LOCATION, stat);

  return stat;
}

int ParallelIO::inq_attname_impl(const std::string &variable_name,
                                 unsigned int n, std::string &result) const {
  int stat;
  char name[PIO_MAX_NAME];
  memset(name, 0, PIO_MAX_NAME);

  int varid = -1;

  if (variable_name == "PISM_GLOBAL") {
    varid = NC_GLOBAL;
  } else {
    stat = PIOc_inq_varid(m_file_id, variable_name.c_str(), &varid);
    check(PISM_ERROR_LOCATION, stat);
  }

  stat = PIOc_inq_attname(m_file_id, varid, n, name);
  check(PISM_ERROR_LOCATION, stat);

  result = name;

  return stat;
}

int ParallelIO::inq_atttype_impl(const std::string &variable_name,
                                 const std::string &att_name,
                                 IO_Type &result) const {
  int stat, varid = -1;
  nc_type tmp;

  if (variable_name == "PISM_GLOBAL") {
    varid = NC_GLOBAL;
  } else {
    stat = PIOc_inq_varid(m_file_id, variable_name.c_str(), &varid);
    check(PISM_ERROR_LOCATION, stat);
  }

  stat = PIOc_inq_atttype(m_file_id, varid, att_name.c_str(), &tmp);
  if (stat == PIO_ENOTATT) {
    tmp = NC_NAT;
  } else {
    check(PISM_ERROR_LOCATION, stat);
  }

  result = nc_type_to_pism_type(tmp);

  return 0;
}

// misc
int ParallelIO::set_fill_impl(int fillmode, int &old_modep) const {

  int stat = PIOc_set_fill(m_file_id, fillmode, &old_modep);
  check(PISM_ERROR_LOCATION, stat);

  return stat;

}

int ParallelIO::del_att_impl(const std::string &variable_name, const std::string &att_name) const {
  int stat, varid = -1;

  if (variable_name == "PISM_GLOBAL") {
    varid = NC_GLOBAL;
  } else {
    stat = PIOc_inq_varid(m_file_id, variable_name.c_str(), &varid);
    check(PISM_ERROR_LOCATION, stat);
  }

  stat = PIOc_del_att(m_file_id, varid, att_name.c_str());
  check(PISM_ERROR_LOCATION, stat);

  return 0;
}

} // end of namespace io
} // end of namespace pism
