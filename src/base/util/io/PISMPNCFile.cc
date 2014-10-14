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

#include <pnetcdf.h>
#include <sstream>
#include <string.h>

#include "PISMPNCFile.hh"

#include "error_handling.hh"

namespace pism {

#include "pism_type_conversion.hh" // has to go after pnetcdf.h

PNCFile::PNCFile(MPI_Comm c)
  : NCFile(c) {
  // MPI_Info_create(&mpi_info);
  mpi_info = MPI_INFO_NULL;
}


PNCFile::~PNCFile() {
  // MPI_Info_free(&mpi_info);
}

void PNCFile::check(int return_code) const {
  if (return_code != NC_NOERR) {
    throw RuntimeError(ncmpi_strerror(return_code));
  }
}

int PNCFile::integer_open_mode(IO_Mode input) const {
  if (input == PISM_READONLY) {
    return NC_NOWRITE;
  } else {
    return NC_WRITE;
  }
}

int PNCFile::open_impl(const std::string &fname, IO_Mode mode) {
  int stat;

  init_hints();

  m_filename = fname;

  int nc_mode = integer_open_mode(mode);
  stat = ncmpi_open(m_com, m_filename.c_str(), nc_mode, mpi_info, &m_file_id); check(stat);

  m_define_mode = false;

  return stat;
}


int PNCFile::create_impl(const std::string &fname) {
  int stat;

  init_hints();

  m_filename = fname;

  stat = ncmpi_create(m_com, m_filename.c_str(), NC_CLOBBER|NC_64BIT_OFFSET,
                      mpi_info, &m_file_id); check(stat);
  m_define_mode = true;

  return stat;
}


int PNCFile::close_impl() {
  int stat = ncmpi_close(m_file_id); check(stat);

  m_file_id = -1;

  m_filename.clear();

  return stat;
}


int PNCFile::enddef_impl() const {

  if (m_define_mode == false)
    return 0;

  int stat = ncmpi_enddef(m_file_id); check(stat);

  m_define_mode = false;

  return stat;
}


int PNCFile::redef_impl() const {

  if (m_define_mode == true)
    return 0;

  int stat = ncmpi_redef(m_file_id); check(stat);

  m_define_mode = true;

  return stat;
}


int PNCFile::def_dim_impl(const std::string &name, size_t length) const {
  int dimid = 0, stat;

  stat = ncmpi_def_dim(m_file_id, name.c_str(), length, &dimid); check(stat);

  return stat;
}


int PNCFile::inq_dimid_impl(const std::string &dimension_name, bool &exists) const {
  int tmp, stat;

  stat = ncmpi_inq_dimid(m_file_id, dimension_name.c_str(), &tmp);

  if (stat == NC_NOERR) {
    exists = true;
  } else {
    exists = false;
  }

  return 0;
}


int PNCFile::inq_dimlen_impl(const std::string &dimension_name, unsigned int &result) const {
  int stat, dimid = -1;
  MPI_Offset len;

  stat = ncmpi_inq_dimid(m_file_id, dimension_name.c_str(), &dimid);

  stat = ncmpi_inq_dimlen(m_file_id, dimid, &len); check(stat);

  result = static_cast<unsigned int>(len);

  return stat;
}


int PNCFile::inq_unlimdim_impl(std::string &result) const {
  int stat, dimid;
  char dimname[NC_MAX_NAME];

  stat = ncmpi_inq_unlimdim(m_file_id, &dimid); check(stat);

  if (dimid == -1) {
    result.clear();
  } else {
    stat = ncmpi_inq_dimname(m_file_id, dimid, dimname); check(stat);

    result = dimname;
  }

  return stat;
}

int PNCFile::inq_dimname_impl(int j, std::string &result) const {
  int stat;
  char dimname[NC_MAX_NAME];
  memset(dimname, 0, NC_MAX_NAME);

  stat = ncmpi_inq_dimname(m_file_id, j, dimname); check(stat);

  result = dimname;

  return stat;
}


int PNCFile::inq_ndims_impl(int &result) const {
  int stat;

  stat = ncmpi_inq_ndims(m_file_id, &result); check(stat);

  return stat;
}

int PNCFile::def_var_impl(const std::string &name, IO_Type nctype, const std::vector<std::string> &dims) const {
  std::vector<int> dimids;
  int stat, varid;

  std::vector<std::string>::const_iterator j;
  for (j = dims.begin(); j != dims.end(); ++j) {
    int dimid = -1;
    stat = ncmpi_inq_dimid(m_file_id, j->c_str(), &dimid); check(stat);
    dimids.push_back(dimid);
  }

  stat = ncmpi_def_var(m_file_id, name.c_str(), pism_type_to_nc_type(nctype),
                       static_cast<int>(dims.size()), &dimids[0], &varid); check(stat);

#if (PISM_DEBUG==1)
  if (stat != NC_NOERR) {
    fprintf(stderr, "def_var: filename = %s, var = %s, dims:", m_filename.c_str(),
            name.c_str());
    for (unsigned int k = 0; k < dims.size(); ++k) {
      fprintf(stderr, "%s(%d), ", dims[k].c_str(), dimids[k]);
    }
    fprintf(stderr, "\n");
  }
#endif

  return stat;
}


int PNCFile::get_vara_double_impl(const std::string &variable_name,
                                 const std::vector<unsigned int> &start,
                                 const std::vector<unsigned int> &count,
                                 double *ip) const {
  std::vector<unsigned int> dummy;
  return this->get_var_double(variable_name,
                              start, count, dummy, ip, false);
}


int PNCFile::put_vara_double_impl(const std::string &variable_name,
                                 const std::vector<unsigned int> &start,
                                 const std::vector<unsigned int> &count,
                                 const double *op) const {
  std::vector<unsigned int> dummy;
  return this->put_var_double(variable_name,
                              start, count, dummy, op, false);
}


int PNCFile::get_varm_double_impl(const std::string &variable_name,
                                 const std::vector<unsigned int> &start,
                                 const std::vector<unsigned int> &count,
                                 const std::vector<unsigned int> &imap,
                                 double *ip) const {
  return this->get_var_double(variable_name,
                              start, count, imap, ip, true);
}


int PNCFile::put_varm_double_impl(const std::string &variable_name,
                                 const std::vector<unsigned int> &start,
                                 const std::vector<unsigned int> &count,
                                 const std::vector<unsigned int> &imap,
                                 const double *op) const {
  return this->put_var_double(variable_name,
                              start, count, imap, op, true);
}


int PNCFile::inq_nvars_impl(int &result) const {
  int stat;

  stat = ncmpi_inq_nvars(m_file_id, &result); check(stat);

  return stat;
}


int PNCFile::inq_vardimid_impl(const std::string &variable_name, std::vector<std::string> &result) const {
  int stat, ndims, varid = -1;
  std::vector<int> dimids;

  stat = ncmpi_inq_varid(m_file_id, variable_name.c_str(), &varid); check(stat);

  stat = ncmpi_inq_varndims(m_file_id, varid, &ndims); check(stat);

  if (ndims == 0) {
    result.clear();
    return 0;
  }

  result.resize(ndims);
  dimids.resize(ndims);

  stat = ncmpi_inq_vardimid(m_file_id, varid, &dimids[0]); check(stat);

  for (int k = 0; k < ndims; ++k) {
    char name[NC_MAX_NAME];
    memset(name, 0, NC_MAX_NAME);

    stat = ncmpi_inq_dimname(m_file_id, dimids[k], name); check(stat);

    result[k] = name;
  }

  return 0;
}


int PNCFile::inq_varnatts_impl(const std::string &variable_name, int &result) const {
  int stat, varid = -1;

  if (variable_name == "PISM_GLOBAL") {
    varid = NC_GLOBAL;
  } else {
    stat = ncmpi_inq_varid(m_file_id, variable_name.c_str(), &varid); check(stat);
  }

  stat = ncmpi_inq_varnatts(m_file_id, varid, &result); check(stat);

  return 0;
}


int PNCFile::inq_varid_impl(const std::string &variable_name, bool &exists) const {
  int stat, flag = -1;

  stat = ncmpi_inq_varid(m_file_id, variable_name.c_str(), &flag);

  if (stat == NC_NOERR)
    flag = 1;
  else
    flag = 0;

  exists = (flag == 1);

  return 0;
}


int PNCFile::inq_varname_impl(unsigned int j, std::string &result) const {
  int stat;
  char varname[NC_MAX_NAME];
  memset(varname, 0, NC_MAX_NAME);

  stat = ncmpi_inq_varname(m_file_id, j, varname); check(stat);

  result = varname;

  return stat;
}

int PNCFile::inq_vartype_impl(const std::string &variable_name, IO_Type &result) const {
  int stat, varid;
  nc_type var_type;

  stat = ncmpi_inq_varid(m_file_id, variable_name.c_str(), &varid); check(stat);

  stat = ncmpi_inq_vartype(m_file_id, varid, &var_type); check(stat);

  result = nc_type_to_pism_type(var_type);

  return 0;
}


int PNCFile::get_att_double_impl(const std::string &variable_name, const std::string &att_name, std::vector<double> &result) const {
  int stat, len, varid = -1;
  MPI_Offset attlen;

  // Read the attribute length:
  if (variable_name == "PISM_GLOBAL") {
    varid = NC_GLOBAL;
  } else {
    stat = ncmpi_inq_varid(m_file_id, variable_name.c_str(), &varid); check(stat);
  }

  stat = ncmpi_inq_attlen(m_file_id, varid, att_name.c_str(), &attlen);

  if (stat == NC_NOERR)
    len = static_cast<int>(attlen);
  else if (stat == NC_ENOTATT)
    len = 0;
  else {
    check(stat);
    len = 0;
  }

  if (len == 0) {
    result.clear();
    return 0;
  }

  result.resize(len);

  // Now read data and broadcast stat to see if we succeeded:
  stat = ncmpi_get_att_double(m_file_id, varid, att_name.c_str(), &result[0]); check(stat);

  // On error, print a message and stop.
  if (stat != NC_NOERR) {
    fprintf(stderr, "Error reading the %s attribute; (varid %d, NetCDF error %s)",
            att_name.c_str(), varid, ncmpi_strerror(stat));
  }

  return 0;
}


int PNCFile::get_att_text_impl(const std::string &variable_name, const std::string &att_name, std::string &result) const {
  char *str = NULL;
  int stat, len, varid = -1;

  // Read the attribute length:
  MPI_Offset attlen;

  if (variable_name == "PISM_GLOBAL") {
    varid = NC_GLOBAL;
  } else {
    stat = ncmpi_inq_varid(m_file_id, variable_name.c_str(), &varid); check(stat);
  }

  stat = ncmpi_inq_attlen(m_file_id, varid, att_name.c_str(), &attlen);
  if (stat == NC_NOERR)
    len = static_cast<int>(attlen);
  else
    len = 0;

  // Allocate some memory or set result to NULL and return:
  if (len == 0) {
    result.clear();
    return 0;
  }

  str = new char[len + 1];
  memset(str, 0, len + 1);

  // Now read the string and see if we succeeded:
  stat = ncmpi_get_att_text(m_file_id, varid, att_name.c_str(), str);

  // On success, broadcast the string. On error, set str to "".
  if (stat != NC_NOERR) {
    fprintf(stderr, "Error reading the %s attribute; (variable %s, %s)",
            att_name.c_str(), variable_name.c_str(), ncmpi_strerror(stat));

    delete[] str;
    return stat;
  }

  result = str;

  delete[] str;
  return 0;
}


int PNCFile::put_att_double_impl(const std::string &variable_name, const std::string &att_name, IO_Type nctype, const std::vector<double> &data) const {
  int stat = 0;

  stat = redef(); check(stat);

  int varid = -1;

  if (variable_name == "PISM_GLOBAL") {
    varid = NC_GLOBAL;
  } else {
    stat = ncmpi_inq_varid(m_file_id, variable_name.c_str(), &varid); check(stat);
  }

  stat = ncmpi_put_att_double(m_file_id, varid, att_name.c_str(),
                              pism_type_to_nc_type(nctype), data.size(), &data[0]); check(stat);

  return stat;
}


int PNCFile::put_att_text_impl(const std::string &variable_name, const std::string &att_name, const std::string &value) const {
  int stat = 0, varid = -1;

  stat = redef(); check(stat);

  if (variable_name == "PISM_GLOBAL") {
    varid = NC_GLOBAL;
  } else {
    stat = ncmpi_inq_varid(m_file_id, variable_name.c_str(), &varid); check(stat);
  }

  stat = ncmpi_put_att_text(m_file_id, varid, att_name.c_str(), value.size(), value.c_str()); check(stat);

  return stat;
}


int PNCFile::inq_attname_impl(const std::string &variable_name, unsigned int n, std::string &result) const {
  int stat;
  char name[NC_MAX_NAME];
  memset(name, 0, NC_MAX_NAME);

  int varid = -1;

  if (variable_name == "PISM_GLOBAL") {
    varid = NC_GLOBAL;
  } else {
    stat = ncmpi_inq_varid(m_file_id, variable_name.c_str(), &varid); check(stat);
  }

  stat = ncmpi_inq_attname(m_file_id, varid, n, name); check(stat);

  result = name;

  return stat;
}


int PNCFile::inq_atttype_impl(const std::string &variable_name, const std::string &att_name, IO_Type &result) const {
  int stat, varid = -1;
  nc_type tmp;

  if (variable_name == "PISM_GLOBAL") {
    varid = NC_GLOBAL;
  } else {
    stat = ncmpi_inq_varid(m_file_id, variable_name.c_str(), &varid); check(stat);
  }

  stat = ncmpi_inq_atttype(m_file_id, varid, att_name.c_str(), &tmp);
  if (stat == NC_ENOTATT) {
    tmp = NC_NAT;
  } else {
    check(stat);
  }

  result = nc_type_to_pism_type(tmp);

  return 0;
}


int PNCFile::set_fill_impl(int fillmode, int &old_modep) const {
  int stat;

  stat = ncmpi_set_fill(m_file_id, fillmode, &old_modep); check(stat);

  return stat;
}


int PNCFile::get_var_double(const std::string &variable_name,
                            const std::vector<unsigned int> &start,
                            const std::vector<unsigned int> &count,
                            const std::vector<unsigned int> &imap_input, double *ip,
                            bool mapped) const {
  std::vector<unsigned int> imap = imap_input;
  int stat, varid, ndims = static_cast<int>(start.size());

#if (PISM_DEBUG==1)
  if (mapped) {
    if (start.size() != count.size() ||
        start.size() != imap.size()) {
      fprintf(stderr, "start, count and imap arrays have to have the same size\n");
      return NC_EINVAL;           // invalid argument error code
    }
  } else {
    if (start.size() != count.size()) {
      fprintf(stderr, "start and count arrays have to have the same size\n");
      return NC_EINVAL;           // invalid argument error code
    }
  }
#endif

  if (mapped == false)
    imap.resize(ndims);

  std::vector<MPI_Offset> nc_start(ndims), nc_count(ndims),
    nc_imap(ndims), nc_stride(ndims);

  stat = ncmpi_inq_varid(m_file_id, variable_name.c_str(), &varid); check(stat);

  for (int j = 0; j < ndims; ++j) {
    nc_start[j] = start[j];
    nc_count[j] = count[j];
    nc_imap[j]  = imap[j];
    nc_stride[j] = 1;
  }

  if (mapped) {
    stat = ncmpi_get_varm_double_all(m_file_id, varid,
                                     &nc_start[0], &nc_count[0], &nc_stride[0], &nc_imap[0],
                                     ip); check(stat);
  } else {
    stat = ncmpi_get_vara_double_all(m_file_id, varid,
                                     &nc_start[0], &nc_count[0],
                                     ip); check(stat);
  }

  return stat;
}


int PNCFile::put_var_double(const std::string &variable_name,
                            const std::vector<unsigned int> &start,
                            const std::vector<unsigned int> &count,
                            const std::vector<unsigned int> &imap_input, const double *op,
                            bool mapped) const {
  int stat, varid, ndims = static_cast<int>(start.size());
  std::vector<unsigned int> imap = imap_input;
#if (PISM_DEBUG==1)
  if (mapped) {
    if (start.size() != count.size() ||
        start.size() != imap.size()) {
      fprintf(stderr, "start, count and imap arrays have to have the same size\n");
      return NC_EINVAL;           // invalid argument error code
    }
  } else {
    if (start.size() != count.size()) {
      fprintf(stderr, "start and count arrays have to have the same size\n");
      return NC_EINVAL;           // invalid argument error code
    }
  }
#endif

  if (mapped == false)
    imap.resize(ndims);

  std::vector<MPI_Offset> nc_start(ndims), nc_count(ndims),
    nc_imap(ndims), nc_stride(ndims);

  stat = ncmpi_inq_varid(m_file_id, variable_name.c_str(), &varid); check(stat);

  for (int j = 0; j < ndims; ++j) {
    nc_start[j] = start[j];
    nc_count[j] = count[j];
    nc_imap[j]  = imap[j];
    nc_stride[j] = 1;
  }

  if (mapped) {
    stat = ncmpi_put_varm_double_all(m_file_id, varid,
                                     &nc_start[0], &nc_count[0], &nc_stride[0], &nc_imap[0],
                                     op); check(stat);
  } else {
    stat = ncmpi_put_vara_double_all(m_file_id, varid,
                                     &nc_start[0], &nc_count[0],
                                     op); check(stat);
  }

  return stat;
}


void PNCFile::init_hints() {

  std::vector<std::string>::iterator j = mpi_io_hints.begin();
  while (j != mpi_io_hints.end()) {
    std::istringstream arg(*j);
    std::vector<std::string> words;
    std::string word;
    while (getline(arg, word, ':'))
      words.push_back(word);

    if (words.size() == 2) {
      // printf("Setting MPI I/O hint \"%s\" to \"%s\"...\n",
      //        words[0].c_str(), words[1].c_str());

      MPI_Info_set(mpi_info,
                   const_cast<char*>(words[0].c_str()),
                   const_cast<char*>(words[1].c_str()));
    } else {
      int rank = 0;
      MPI_Comm_rank(m_com, &rank);
      if (rank == 0) {
        printf("PISM WARNING: invalid MPI I/O hint: %s. Ignoring it...\n",
               j->c_str());
      }
    }

    ++j;
  }

}

std::string PNCFile::get_format_impl() const {
  return "netcdf3";
}

} // end of namespace pism
