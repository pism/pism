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
#include "pism/util/IceGrid.hh"
#include "File.hh"

namespace pism {
namespace io {

static void check(const ErrorLocation &where, int return_code) {
  char message[PIO_MAX_NAME + 1];
  if (return_code != PIO_NOERR) {
    PIOc_strerror(return_code, message);
    throw RuntimeError(where, message);
  }
}

IO_Backend ParallelIO::best_iotype(bool netcdf3) {
  if (PIOc_iotype_available(PIO_IOTYPE_PNETCDF) and netcdf3) {
    return PISM_PIO_PNETCDF;
  }
  if (PIOc_iotype_available(PIO_IOTYPE_NETCDF4P)) {
    return PISM_PIO_NETCDF4P;
  }
  if (PIOc_iotype_available(PIO_IOTYPE_NETCDF4C)) {
    return PISM_PIO_NETCDF4C;
  }
  // always available and supports all file formats
  return PISM_PIO_NETCDF;
}


ParallelIO::ParallelIO(MPI_Comm com, int iosysid, IO_Backend iotype)
  : NCFile(com),
    m_iosysid(iosysid) {
  assert(iosysid != -1);

  switch (iotype) {
  case PISM_PIO_PNETCDF:
    m_iotype = PIO_IOTYPE_PNETCDF;
    break;
  case PISM_PIO_NETCDF4P:
    m_iotype = PIO_IOTYPE_NETCDF4P;
    break;
  case PISM_PIO_NETCDF4C:
    m_iotype = PIO_IOTYPE_NETCDF4C;
    break;
  case PISM_PIO_NETCDF:
    m_iotype = PIO_IOTYPE_NETCDF;
    break;
  default:
    throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                  "invalid iotype in ParallelIO");
  }
}

ParallelIO::~ParallelIO() {
  // empty
}

void ParallelIO::open_impl(const std::string &filename, IO_Mode mode, const std::map<std::string, int> &varsi, int FileID) {
  int open_mode = mode == PISM_READONLY ? PIO_NOWRITE : PIO_WRITE;

  int stat = PIOc_open(m_iosysid, filename.c_str(), open_mode, &m_file_id);
  check(PISM_ERROR_LOCATION, stat);
}

void ParallelIO::create_impl(const std::string &filename, int FileID) {

  int mode = NC_CLOBBER;
  if (m_iotype == PIO_IOTYPE_PNETCDF) {
    mode |= NC_64BIT_DATA;
  }

  int stat = PIOc_createfile(m_iosysid, &m_file_id, &m_iotype, filename.c_str(),
                             mode);
  check(PISM_ERROR_LOCATION, stat);
}

void ParallelIO::sync_impl() const {

  int stat = PIOc_sync(m_file_id); check(PISM_ERROR_LOCATION, stat);
}

void ParallelIO::close_impl() {
  int stat = PIOc_closefile(m_file_id); check(PISM_ERROR_LOCATION, stat);
}

// redef/enddef
void ParallelIO::enddef_impl() const {
  int stat = PIOc_enddef(m_file_id); check(PISM_ERROR_LOCATION, stat);
}

void ParallelIO::redef_impl() const {
  int stat = PIOc_redef(m_file_id); check(PISM_ERROR_LOCATION, stat);
}

// dim
void ParallelIO::def_dim_impl(const std::string &name, size_t length) const {
  int dim_id = 0;
  int stat = PIOc_def_dim(m_file_id, name.c_str(), length, &dim_id);
  check(PISM_ERROR_LOCATION, stat);
}

void ParallelIO::inq_dimid_impl(const std::string &dimension_name, bool &exists) const {
  int tmp;

  int stat = PIOc_inq_dimid(m_file_id, dimension_name.c_str(), &tmp);

  if (stat == PIO_NOERR) {
    exists = true;
  } else {
    exists = false;
  }
}

void ParallelIO::inq_dimlen_impl(const std::string &dimension_name, unsigned int &result) const {
  int stat, dimid = -1;
  PIO_Offset len;

  stat = PIOc_inq_dimid(m_file_id, dimension_name.c_str(), &dimid); check(PISM_ERROR_LOCATION, stat);

  stat = PIOc_inq_dimlen(m_file_id, dimid, &len); check(PISM_ERROR_LOCATION, stat);

  result = static_cast<unsigned int>(len);
}

void ParallelIO::inq_unlimdim_impl(std::string &result) const {
  int stat = PIO_NOERR, dimid = -1;
  char dimname[PIO_MAX_NAME + 1];

  stat = PIOc_inq_unlimdim(m_file_id, &dimid); check(PISM_ERROR_LOCATION, stat);

  if (dimid == -1) {
    result.clear();
  } else {
    stat = PIOc_inq_dimname(m_file_id, dimid, dimname); check(PISM_ERROR_LOCATION, stat);

    result = dimname;
  }
}

// var
void ParallelIO::def_var_impl(const std::string &name, IO_Type nctype,
                             const std::vector<std::string> &dims) const {
  std::vector<int> dimids;
  int stat, varid;

  for (auto d : dims) {
    int dimid = -1;
    stat = PIOc_inq_dimid(m_file_id, d.c_str(), &dimid); check(PISM_ERROR_LOCATION, stat);
    dimids.push_back(dimid);
  }

  stat = PIOc_def_var(m_file_id, name.c_str(), pism_type_to_nc_type(nctype),
                      static_cast<int>(dims.size()), dimids.data(), &varid);
  check(PISM_ERROR_LOCATION, stat);
}

void ParallelIO::def_var_chunking_impl(const std::string &name,
                                      std::vector<size_t> &dimensions) const {
  (void) name;
  (void) dimensions;
  // FIXME
}

void ParallelIO::get_vara_double_impl(const std::string &variable_name,
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
}

void ParallelIO::put_vara_double_impl(const std::string &variable_name,
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
}

template<typename T>
std::vector<T> convert_data(const double *input, size_t length) {
  std::vector<T> buffer(length);
  for (size_t k = 0; k < length; ++k) {
    buffer[k] = static_cast<T>(input[k]);
  }
  return buffer;
}

void ParallelIO::write_darray_impl(const std::string &variable_name,
                                   const IceGrid &grid,
                                   unsigned int z_count,
                                   unsigned int record,
                                   const double *input) {

  int stat = 0, varid;

  stat = PIOc_inq_varid(m_file_id, variable_name.c_str(), &varid);
  check(PISM_ERROR_LOCATION, stat);

  int type = 0;
  stat = PIOc_inq_vartype(m_file_id, varid, &type);

  int decompid = grid.pio_io_decomposition(z_count, type);

  size_t length = grid.xm() * grid.ym() * z_count;

  switch (type) {
  case PIO_DOUBLE:
    // no conversion necessary
    stat = PIOc_put_vard_double(m_file_id, varid, decompid, (PIO_Offset)record, input);
    check(PISM_ERROR_LOCATION, stat);
    break;
  case PIO_FLOAT:
    {
      auto buffer = convert_data<float>(input, length);
      stat = PIOc_put_vard_float(m_file_id, varid, decompid, (PIO_Offset)record, buffer.data());
      check(PISM_ERROR_LOCATION, stat);
      break;
    }
  case PIO_INT:
    {
      auto buffer = convert_data<int>(input, length);
      stat = PIOc_put_vard_int(m_file_id, varid, decompid, (PIO_Offset)record, buffer.data());
      check(PISM_ERROR_LOCATION, stat);
      break;
    }
  default:
    throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                  "ParallelIO: type conversion is not supported");
  }
}


void ParallelIO::get_varm_double_impl(const std::string &variable_name,
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

void ParallelIO::inq_nvars_impl(int &result) const {
  int stat = PIOc_inq_nvars(m_file_id, &result); check(PISM_ERROR_LOCATION, stat);
}

void ParallelIO::inq_vardimid_impl(const std::string &variable_name, std::vector<std::string> &result) const {
  int stat, ndims, varid = -1;
  std::vector<int> dimids;

  stat = PIOc_inq_varid(m_file_id, variable_name.c_str(), &varid); check(PISM_ERROR_LOCATION, stat);

  stat = PIOc_inq_varndims(m_file_id, varid, &ndims); check(PISM_ERROR_LOCATION, stat);

  if (ndims == 0) {
    result.clear();
    return;
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
}

void ParallelIO::inq_varnatts_impl(const std::string &variable_name, int &result) const {
  int varid = get_varid(variable_name);

  int stat = PIO_NOERR;
  if (varid == PIO_GLOBAL) {
    stat = PIOc_inq_natts(m_file_id, &result);
  } else {
    stat = PIOc_inq_varnatts(m_file_id, varid, &result);
  }
  check(PISM_ERROR_LOCATION, stat);
}

void ParallelIO::inq_varid_impl(const std::string &variable_name, bool &exists) const {
  int stat, flag = -1;

  stat = PIOc_inq_varid(m_file_id, variable_name.c_str(), &flag);

  exists = (stat == PIO_NOERR);
}

void ParallelIO::inq_varname_impl(unsigned int j, std::string &result) const {
  int stat;
  char varname[PIO_MAX_NAME];
  memset(varname, 0, PIO_MAX_NAME);

  stat = PIOc_inq_varname(m_file_id, j, varname); check(PISM_ERROR_LOCATION, stat);

  result = varname;
}

// att
void ParallelIO::get_att_double_impl(const std::string &variable_name,
                                    const std::string &att_name,
                                    std::vector<double> &result) const {
  // Read the attribute length:
  int varid = get_varid(variable_name);

  PIO_Offset attlen = 0;
  int stat = PIOc_inq_attlen(m_file_id, varid, att_name.c_str(), &attlen);

  int len = 0;
  if (stat == NC_NOERR) {
    len = static_cast<int>(attlen);
  } else if (stat == NC_ENOTATT) {
    len = 0;
  } else {
    check(PISM_ERROR_LOCATION, stat);
  }

  if (len == 0) {
    result.clear();
    return;
  }

  result.resize(len);

  stat = PIOc_get_att_double(m_file_id, varid, att_name.c_str(), result.data());
  check(PISM_ERROR_LOCATION, stat);
}

void ParallelIO::get_att_text_impl(const std::string &variable_name,
                                   const std::string &att_name,
                                   std::string &result) const {
  PIO_Offset attlen;

  int varid = get_varid(variable_name);

  int stat = PIOc_inq_attlen(m_file_id, varid, att_name.c_str(), &attlen);

  int len = 0;
  if (stat == NC_NOERR) {
    len = static_cast<int>(attlen);
  } else {
    len = 0;
  }

  // Allocate some memory or clear result.
  if (len == 0) {
    result.clear();
    return;
  }

  // Now read the string and see if we succeeded:
  std::vector<char> str(len + 1, 0);
  stat = PIOc_get_att_text(m_file_id, varid, att_name.c_str(), str.data());
  check(PISM_ERROR_LOCATION, stat);

  result = str.data();
}

void ParallelIO::put_att_double_impl(const std::string &variable_name,
                                    const std::string &att_name,
                                    IO_Type xtype,
                                    const std::vector<double> &data) const {
  int stat = PIOc_put_att_double(m_file_id, get_varid(variable_name), att_name.c_str(),
                                 pism_type_to_nc_type(xtype), data.size(), data.data());
  check(PISM_ERROR_LOCATION, stat);
}

void ParallelIO::put_att_text_impl(const std::string &variable_name,
                                  const std::string &att_name,
                                  const std::string &value) const {
  int stat = PIOc_put_att_text(m_file_id, get_varid(variable_name), att_name.c_str(),
                               value.size(), value.c_str());
  check(PISM_ERROR_LOCATION, stat);
}

void ParallelIO::inq_attname_impl(const std::string &variable_name,
                                 unsigned int n, std::string &result) const {
  std::vector<char> name(PIO_MAX_NAME + 1, 0);

  int stat = PIOc_inq_attname(m_file_id, get_varid(variable_name), n, name.data());
  check(PISM_ERROR_LOCATION, stat);

  result = name.data();
}

void ParallelIO::inq_atttype_impl(const std::string &variable_name,
                                 const std::string &att_name,
                                 IO_Type &result) const {
  nc_type tmp;
  int stat = PIOc_inq_atttype(m_file_id, get_varid(variable_name), att_name.c_str(), &tmp);
  if (stat == PIO_ENOTATT) {
    tmp = NC_NAT;
  } else {
    check(PISM_ERROR_LOCATION, stat);
  }

  result = nc_type_to_pism_type(tmp);
}

// misc
void ParallelIO::set_fill_impl(int fillmode, int &old_modep) const {

  int stat = PIOc_set_fill(m_file_id, fillmode, &old_modep);
  check(PISM_ERROR_LOCATION, stat);
}

void ParallelIO::del_att_impl(const std::string &variable_name,
                              const std::string &att_name) const {
  int stat = PIOc_del_att(m_file_id, get_varid(variable_name), att_name.c_str());
  check(PISM_ERROR_LOCATION, stat);
}

int ParallelIO::get_varid(const std::string &variable_name) const {
  if (variable_name == "PISM_GLOBAL") {
    return NC_GLOBAL;
  } else {
    int varid = -2;
    int stat = PIOc_inq_varid(m_file_id, variable_name.c_str(), &varid);
    check(PISM_ERROR_LOCATION, stat);
    return varid;
  }
}

void ParallelIO::create_grid_impl(int lengthx, int lengthy) const {
  (void) lengthx;
  (void) lengthy;
}

void ParallelIO::define_timestep_impl(int tsID) const {
  (void) tsID;
}

void ParallelIO::write_timestep_impl() const {
}

void ParallelIO::def_ref_date_impl(double time) const {
  (void) time;
}

} // end of namespace io
} // end of namespace pism
