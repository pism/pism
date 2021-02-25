// Copyright (C) 2012, 2013, 2014, 2015, 2016, 2017, 2019 PISM Authors
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

#include <pnetcdf.h>
#include <sstream>
#include <string.h>

#include "PNCFile.hh"

#include "pism_type_conversion.hh" // has to go after pnetcdf.h
#include "pism_cdi_type_conversion.hh"

#include "pism/util/error_handling.hh"

namespace pism {
namespace io {

PNCFile::PNCFile(MPI_Comm c)
  : NCFile(c) {
  MPI_Info_create(&m_mpi_info);
}


PNCFile::~PNCFile() {
  MPI_Info_free(&m_mpi_info);
}

static void check(const ErrorLocation &where, int return_code) {
  if (return_code != NC_NOERR) {
    throw RuntimeError(where, ncmpi_strerror(return_code));
  }
}

void PNCFile::open_impl(const std::string &fname,
                        IO_Mode mode,
                        int FileID,
                        const std::map<std::string, int> &dimsa) {
  int stat;

  init_hints();

  int open_mode = mode == PISM_READONLY ? NC_NOWRITE : NC_WRITE;

  stat = ncmpi_open(m_com, fname.c_str(), open_mode, m_mpi_info, &m_file_id);
  check(PISM_ERROR_LOCATION, stat);
}


void PNCFile::create_impl(const std::string &fname, int FileID) {
  int stat;

  init_hints();

  stat = ncmpi_create(m_com, fname.c_str(), NC_CLOBBER|NC_64BIT_DATA,
                      m_mpi_info, &m_file_id); check(PISM_ERROR_LOCATION, stat);
}

void PNCFile::sync_impl() const {

  int stat = ncmpi_sync(m_file_id); check(PISM_ERROR_LOCATION, stat);
}


void PNCFile::close_impl() {
  int stat = ncmpi_close(m_file_id); check(PISM_ERROR_LOCATION, stat);

  m_file_id = -1;
}


void PNCFile::enddef_impl() const {

  int stat = ncmpi_enddef(m_file_id); check(PISM_ERROR_LOCATION, stat);
}


void PNCFile::redef_impl() const {

  int stat = ncmpi_redef(m_file_id); check(PISM_ERROR_LOCATION, stat);
}


void PNCFile::def_dim_impl(const std::string &name, size_t length, int dim) const {
  int dimid = 0, stat;

  stat = ncmpi_def_dim(m_file_id, name.c_str(), length, &dimid); check(PISM_ERROR_LOCATION, stat);
}


void PNCFile::inq_dimid_impl(const std::string &dimension_name, bool &exists) const {
  int tmp, stat;

  stat = ncmpi_inq_dimid(m_file_id, dimension_name.c_str(), &tmp);

  if (stat == NC_NOERR) {
    exists = true;
  } else {
    exists = false;
  }
}


void PNCFile::inq_dimlen_impl(const std::string &dimension_name, unsigned int &result) const {
  int stat, dimid = -1;
  MPI_Offset len;

  stat = ncmpi_inq_dimid(m_file_id, dimension_name.c_str(), &dimid); check(PISM_ERROR_LOCATION, stat);

  stat = ncmpi_inq_dimlen(m_file_id, dimid, &len); check(PISM_ERROR_LOCATION, stat);

  result = static_cast<unsigned int>(len);
}


void PNCFile::inq_unlimdim_impl(std::string &result) const {
  int stat = NC_NOERR, dimid = -1;
  char dimname[NC_MAX_NAME];

  stat = ncmpi_inq_unlimdim(m_file_id, &dimid); check(PISM_ERROR_LOCATION, stat);

  if (dimid == -1) {
    result.clear();
  } else {
    stat = ncmpi_inq_dimname(m_file_id, dimid, dimname); check(PISM_ERROR_LOCATION, stat);

    result = dimname;
  }
}

void PNCFile::def_var_impl(const std::string &name, IO_Type nctype, const std::vector<std::string> &dims) const {
  std::vector<int> dimids;
  int stat, varid;

  for (auto d : dims) {
    int dimid = -1;
    stat = ncmpi_inq_dimid(m_file_id, d.c_str(), &dimid); check(PISM_ERROR_LOCATION, stat);
    dimids.push_back(dimid);
  }

  stat = ncmpi_def_var(m_file_id, name.c_str(), pism_type_to_nc_type(nctype),
                       static_cast<int>(dims.size()), &dimids[0], &varid);
  check(PISM_ERROR_LOCATION, stat);
}


void PNCFile::get_vara_double_impl(const std::string &variable_name,
                                 const std::vector<unsigned int> &start,
                                 const std::vector<unsigned int> &count,
                                 double *ip) const {
  std::vector<unsigned int> dummy;
  return this->get_var_double(variable_name,
                              start, count, dummy, ip, false);
}


void PNCFile::put_vara_double_impl(const std::string &variable_name,
                                  const std::vector<unsigned int> &start,
                                  const std::vector<unsigned int> &count,
                                  const double *op) const {
  int stat, varid, ndims = static_cast<int>(start.size());

  std::vector<MPI_Offset> nc_start(ndims), nc_count(ndims), nc_stride(ndims);

  stat = ncmpi_inq_varid(m_file_id, variable_name.c_str(), &varid);
  check(PISM_ERROR_LOCATION, stat);

  for (int j = 0; j < ndims; ++j) {
    nc_start[j]  = start[j];
    nc_count[j]  = count[j];
    nc_stride[j] = 1;
  }

  stat = ncmpi_put_vara_double_all(m_file_id, varid,
                                   &nc_start[0], &nc_count[0],
                                   op); check(PISM_ERROR_LOCATION, stat);
}


void PNCFile::get_varm_double_impl(const std::string &variable_name,
                                 const std::vector<unsigned int> &start,
                                 const std::vector<unsigned int> &count,
                                 const std::vector<unsigned int> &imap,
                                 double *ip) const {
  return this->get_var_double(variable_name,
                              start, count, imap, ip, true);
}


void PNCFile::inq_nvars_impl(int &result) const {
  int stat;

  stat = ncmpi_inq_nvars(m_file_id, &result); check(PISM_ERROR_LOCATION, stat);
}


void PNCFile::inq_vardimid_impl(const std::string &variable_name, std::vector<std::string> &result) const {
  int stat, ndims, varid = -1;
  std::vector<int> dimids;

  stat = ncmpi_inq_varid(m_file_id, variable_name.c_str(), &varid); check(PISM_ERROR_LOCATION, stat);

  stat = ncmpi_inq_varndims(m_file_id, varid, &ndims); check(PISM_ERROR_LOCATION, stat);

  if (ndims == 0) {
    result.clear();
    return;
  }

  result.resize(ndims);
  dimids.resize(ndims);

  stat = ncmpi_inq_vardimid(m_file_id, varid, &dimids[0]); check(PISM_ERROR_LOCATION, stat);

  for (int k = 0; k < ndims; ++k) {
    char name[NC_MAX_NAME];
    memset(name, 0, NC_MAX_NAME);

    stat = ncmpi_inq_dimname(m_file_id, dimids[k], name); check(PISM_ERROR_LOCATION, stat);

    result[k] = name;
  }
}

int PNCFile::get_varid(const std::string &variable_name) const {
  if (variable_name == "PISM_GLOBAL") {
    return NC_GLOBAL;
  } else {
    int varid = -1;
    int stat = ncmpi_inq_varid(m_file_id, variable_name.c_str(), &varid);
    check(PISM_ERROR_LOCATION, stat);

    return varid;
  }
}

void PNCFile::inq_varnatts_impl(const std::string &variable_name, int &result) const {

  int stat = ncmpi_inq_varnatts(m_file_id, get_varid(variable_name), &result);
  check(PISM_ERROR_LOCATION, stat);
}


void PNCFile::inq_varid_impl(const std::string &variable_name, bool &exists) const {
  int stat, flag = -1;

  stat = ncmpi_inq_varid(m_file_id, variable_name.c_str(), &flag);

  exists = (stat == NC_NOERR);
}


void PNCFile::inq_varname_impl(unsigned int j, std::string &result) const {
  int stat;
  char varname[NC_MAX_NAME];
  memset(varname, 0, NC_MAX_NAME);

  stat = ncmpi_inq_varname(m_file_id, j, varname); check(PISM_ERROR_LOCATION, stat);

  result = varname;
}

void PNCFile::get_att_double_impl(const std::string &variable_name,
                                  const std::string &att_name,
                                  std::vector<double> &result) const {
  int len, varid = get_varid(variable_name);
  MPI_Offset attlen = 0;

  // Read the attribute length:
  int stat = ncmpi_inq_attlen(m_file_id, varid, att_name.c_str(), &attlen);

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

  // Now read data and broadcast stat to see if we succeeded:
  stat = ncmpi_get_att_double(m_file_id, varid, att_name.c_str(), &result[0]); check(PISM_ERROR_LOCATION, stat);

  // On error, print a message and stop.
  if (stat != NC_NOERR) {
    fprintf(stderr, "Error reading the %s attribute; (varid %d, NetCDF error %s)",
            att_name.c_str(), varid, ncmpi_strerror(stat));
  }
}


void PNCFile::get_att_text_impl(const std::string &variable_name, const std::string &att_name, std::string &result) const {
  int varid = get_varid(variable_name);

  // Read the attribute length:
  MPI_Offset attlen = 0;

  int stat = ncmpi_inq_attlen(m_file_id, varid, att_name.c_str(), &attlen);
  int len = (stat == NC_NOERR) ? static_cast<int>(attlen) : 0;

  // Allocate some memory or set result to NULL and return:
  if (len == 0) {
    result.clear();
    return;
  }

  std::vector<char>str(len + 1, 0);

  // Now read the string and see if we succeeded:
  stat = ncmpi_get_att_text(m_file_id, varid, att_name.c_str(), str.data());
  check(PISM_ERROR_LOCATION, stat);

  result = str.data();
}

void PNCFile::del_att_impl(const std::string &variable_name, const std::string &att_name) const {
  int stat = ncmpi_del_att(m_file_id, get_varid(variable_name), att_name.c_str());
  check(PISM_ERROR_LOCATION, stat);
}

void PNCFile::put_att_double_impl(const std::string &variable_name, const std::string &att_name, IO_Type nctype, const std::vector<double> &data) const {

  int stat = ncmpi_put_att_double(m_file_id, get_varid(variable_name), att_name.c_str(),
                                  pism_type_to_nc_type(nctype), data.size(), &data[0]);
  check(PISM_ERROR_LOCATION, stat);
}


void PNCFile::put_att_text_impl(const std::string &variable_name,
                                const std::string &att_name,
                                const std::string &value) const {
  int stat = ncmpi_put_att_text(m_file_id, get_varid(variable_name),
                                att_name.c_str(), value.size(), value.c_str());
  check(PISM_ERROR_LOCATION, stat);
}


void PNCFile::inq_attname_impl(const std::string &variable_name,
                               unsigned int n,
                               std::string &result) const {
  int stat;
  char name[NC_MAX_NAME];
  memset(name, 0, NC_MAX_NAME);

  int varid = get_varid(variable_name);

  stat = ncmpi_inq_attname(m_file_id, varid, n, name); check(PISM_ERROR_LOCATION, stat);

  result = name;
}


void PNCFile::inq_atttype_impl(const std::string &variable_name,
                               const std::string &att_name,
                               IO_Type &result) const {
  int varid = get_varid(variable_name);

  nc_type tmp = NC_NAT;
  int stat = ncmpi_inq_atttype(m_file_id, varid, att_name.c_str(), &tmp);
  if (stat == NC_ENOTATT) {
    tmp = NC_NAT;
  } else {
    check(PISM_ERROR_LOCATION, stat);
  }

  result = nc_type_to_pism_type(tmp);
}


void PNCFile::set_fill_impl(int fillmode, int &old_modep) const {
  int stat = ncmpi_set_fill(m_file_id, fillmode, &old_modep); check(PISM_ERROR_LOCATION, stat);
}


void PNCFile::get_var_double(const std::string &variable_name,
                            const std::vector<unsigned int> &start,
                            const std::vector<unsigned int> &count,
                            const std::vector<unsigned int> &imap_input, double *ip,
                            bool transposed) const {
  std::vector<unsigned int> imap = imap_input;
  int stat, varid, ndims = static_cast<int>(start.size());

  if (not transposed) {
    imap.resize(ndims);
  }

  std::vector<MPI_Offset> nc_start(ndims), nc_count(ndims),
    nc_imap(ndims), nc_stride(ndims);

  stat = ncmpi_inq_varid(m_file_id, variable_name.c_str(), &varid); check(PISM_ERROR_LOCATION, stat);

  for (int j = 0; j < ndims; ++j) {
    nc_start[j] = start[j];
    nc_count[j] = count[j];
    nc_imap[j]  = imap[j];
    nc_stride[j] = 1;
  }

  if (transposed) {
    stat = ncmpi_get_varm_double_all(m_file_id, varid,
                                     &nc_start[0], &nc_count[0], &nc_stride[0], &nc_imap[0],
                                     ip);
    check(PISM_ERROR_LOCATION, stat);
  } else {
    stat = ncmpi_get_vara_double_all(m_file_id, varid,
                                     &nc_start[0], &nc_count[0],
                                     ip);
    check(PISM_ERROR_LOCATION, stat);
  }
}

void PNCFile::init_hints() {

  for (auto hint : m_mpi_io_hints) {
    std::istringstream arg(hint);
    std::vector<std::string> words;
    std::string word;
    while (getline(arg, word, ':')) {
      words.push_back(word);
    }

    if (words.size() == 2) {
      // printf("Setting MPI I/O hint \"%s\" to \"%s\"...\n",
      //        words[0].c_str(), words[1].c_str());

      MPI_Info_set(m_mpi_info,
                   const_cast<char*>(words[0].c_str()),
                   const_cast<char*>(words[1].c_str()));
    } else {
      int rank = 0;
      MPI_Comm_rank(m_com, &rank);
      if (rank == 0) {
        printf("PISM WARNING: invalid MPI I/O hint: %s. Ignoring it...\n",
               hint.c_str());
      }
    }
  }
}

void PNCFile::create_grid_impl(int lengthx, int lengthy) const {
  (void) lengthx;
  (void) lengthy;
}

void PNCFile::define_timestep_impl(int tsID) const {
  (void) tsID;
}

void PNCFile::def_ref_date_impl(double time) const {
  (void) time;
}

} // end of namespace io
} // end of namespace pism
