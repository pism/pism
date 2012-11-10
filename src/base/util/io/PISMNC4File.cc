// Copyright (C) 2012 PISM Authors
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

#include "PISMNC4File.hh"

extern "C" {
#include <netcdf_par.h>
}

#include <cstring>		// memset
#include <cstdio>		// stderr, fprintf

#include "pism_type_conversion.hh"

PISMNC4File::PISMNC4File(MPI_Comm c, int r)
  : PISMNCFile(c, r) {
  // empty
}

PISMNC4File::~PISMNC4File() {
  // empty
}

// open/create/close

int PISMNC4File::open(string fname, int mode) {
  MPI_Info info = MPI_INFO_NULL;
  int stat;

  filename = fname;

  stat = nc_open_par(filename.c_str(),
                     mode | NC_MPIIO,
                     com, info, &ncid); check(stat);

  define_mode = false;

  return stat;
}

int PISMNC4File::create(string fname) {
  MPI_Info info = MPI_INFO_NULL;
  int stat;

  filename = fname;

  stat = nc_create_par(filename.c_str(),
                       NC_NETCDF4 | NC_MPIIO,
                       com, info, &ncid); check(stat);
  define_mode = true;

  return stat;
}

int PISMNC4File::close() {
  int stat;

  stat = nc_close(ncid); check(stat);

  ncid = -1;

  filename.clear();

  return stat;
}

// redef/enddef
int PISMNC4File::enddef() const {

  if (define_mode == false)
    return 0;

  int stat = nc_enddef(ncid); check(stat);

  define_mode = false;

  return stat;
}

int PISMNC4File::redef() const {

  if (define_mode == true)
    return 0;

  int stat = nc_redef(ncid); check(stat);

  define_mode = true;

  return stat;
}

// dim
int PISMNC4File::def_dim(string name, size_t length) const {
  int dimid = 0, stat;

  stat = nc_def_dim(ncid, name.c_str(), length, &dimid); check(stat);

  return stat;
}

int PISMNC4File::inq_dimid(string dimension_name, bool &exists) const {
  int tmp, stat;

  stat = nc_inq_dimid(ncid, dimension_name.c_str(), &tmp);

  if (stat == NC_NOERR) {
    exists = true;
  } else {
    exists = false;
  }

  return 0;
}

int PISMNC4File::inq_dimlen(string dimension_name, unsigned int &result) const {
  int stat, dimid = -1;
  size_t len;

  stat = nc_inq_dimid(ncid, dimension_name.c_str(), &dimid);

  stat = nc_inq_dimlen(ncid, dimid, &len); check(stat);

  result = static_cast<unsigned int>(len);

  return stat;
}

int PISMNC4File::inq_unlimdim(string &result) const {
  int stat, dimid;
  char dimname[NC_MAX_NAME];

  stat = nc_inq_unlimdim(ncid, &dimid); check(stat);

  if (dimid == -1) {
    result.clear();
  } else {
    stat = nc_inq_dimname(ncid, dimid, dimname); check(stat);

    result = dimname;
  }

  return stat;
}

// var
int PISMNC4File::def_var(string name, PISM_IO_Type nctype, vector<string> dims) const {
  vector<int> dimids;
  int stat, varid;

  vector<string>::iterator j;
  for (j = dims.begin(); j != dims.end(); ++j) {
    int dimid = -1;
    stat = nc_inq_dimid(ncid, j->c_str(), &dimid); check(stat);
    dimids.push_back(dimid);
  }

  stat = nc_def_var(ncid, name.c_str(), pism_type_to_nc_type(nctype),
                    static_cast<int>(dims.size()), &dimids[0], &varid); check(stat);

#if (PISM_DEBUG==1)
  if (stat != NC_NOERR) {
    fprintf(stderr, "def_var: filename = %s, var = %s, dims:", filename.c_str(),
            name.c_str());
    for (unsigned int k = 0; k < dims.size(); ++k) {
      fprintf(stderr, "%s(%d), ", dims[k].c_str(), dimids[k]);
    }
    fprintf(stderr, "\n");
  }
#endif

  return stat;
}

int PISMNC4File::get_varm_double(string variable_name,
                                 vector<unsigned int> start,
                                 vector<unsigned int> count,
                                 vector<unsigned int> imap, double *op) const {
  return this->get_var_double(variable_name,
                              start, count, imap, op, true);
}

int PISMNC4File::get_vara_double(string variable_name,
                                 vector<unsigned int> start,
                                 vector<unsigned int> count,
                                 double *op) const {
  vector<unsigned int> dummy;
  return this->get_var_double(variable_name,
                              start, count, dummy, op, false);
}

int PISMNC4File::get_var_double(string variable_name,
                                vector<unsigned int> start,
                                vector<unsigned int> count,
                                vector<unsigned int> imap, double *ip,
                                bool mapped) const {
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

  vector<size_t> nc_start(ndims), nc_count(ndims);
  vector<ptrdiff_t> nc_imap(ndims), nc_stride(ndims);

  stat = nc_inq_varid(ncid, variable_name.c_str(), &varid); check(stat);

  for (int j = 0; j < ndims; ++j) {
    nc_start[j] = start[j];
    nc_count[j] = count[j];
    nc_imap[j]  = imap[j];
    nc_stride[j] = 1;
  }


  if (mapped) {
    // Use independent parallel access mode because it works. It would be
    // better to use collective mode, but I/O performance is ruined by
    // "mapping" anyway.

    stat = nc_var_par_access(ncid, varid, NC_INDEPENDENT); check(stat);

    stat = nc_get_varm_double(ncid, varid,
                              &nc_start[0], &nc_count[0], &nc_stride[0], &nc_imap[0],
                              ip); check(stat);
  } else {
    // Use collective parallel access mode because it is faster (and because it
    // works in this case).

    stat = nc_var_par_access(ncid, varid, NC_COLLECTIVE); check(stat);

    stat = nc_get_vara_double(ncid, varid,
                              &nc_start[0], &nc_count[0],
                              ip); check(stat);
  }

  return stat;
}

int PISMNC4File::put_varm_double(string variable_name,
                                 vector<unsigned int> start,
                                 vector<unsigned int> count,
                                 vector<unsigned int> imap, const double *op) const {
  return this->put_var_double(variable_name,
                              start, count, imap, op, true);
}

int PISMNC4File::put_vara_double(string variable_name,
                                 vector<unsigned int> start,
                                 vector<unsigned int> count,
                                 const double *op) const {
  vector<unsigned int> dummy;
  return this->put_var_double(variable_name,
                              start, count, dummy, op, false);
}

int PISMNC4File::put_var_double(string variable_name,
                                vector<unsigned int> start,
                                vector<unsigned int> count,
                                vector<unsigned int> imap, const double *op,
                                bool mapped) const {
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

  vector<size_t> nc_start(ndims), nc_count(ndims);
  vector<ptrdiff_t> nc_imap(ndims), nc_stride(ndims);

  stat = nc_inq_varid(ncid, variable_name.c_str(), &varid); check(stat);

  for (int j = 0; j < ndims; ++j) {
    nc_start[j] = start[j];
    nc_count[j] = count[j];
    nc_imap[j]  = imap[j];
    nc_stride[j] = 1;
  }

  if (mapped) {
    // Use independent parallel access mode because it works. It would be
    // better to use collective mode, but I/O performance is ruined by
    // "mapping" anyway.

    stat = nc_var_par_access(ncid, varid, NC_INDEPENDENT); check(stat);

    stat = nc_put_varm_double(ncid, varid,
                              &nc_start[0], &nc_count[0], &nc_stride[0], &nc_imap[0],
                              op); check(stat);
  } else {
    // Use collective parallel access mode because it is faster (and because it
    // works in this case).

    stat = nc_var_par_access(ncid, varid, NC_COLLECTIVE); check(stat);

    stat = nc_put_vara_double(ncid, varid,
                              &nc_start[0], &nc_count[0],
                              op); check(stat);
  }

  return stat;
}

int PISMNC4File::inq_nvars(int &result) const {
  int stat;

  stat = nc_inq_nvars(ncid, &result); check(stat);

  return stat;
}

int PISMNC4File::inq_vardimid(string variable_name, vector<string> &result) const {
  int stat, ndims, varid = -1;
  vector<int> dimids;

  stat = nc_inq_varid(ncid, variable_name.c_str(), &varid); check(stat);

  stat = nc_inq_varndims(ncid, varid, &ndims); check(stat);

  if (ndims == 0) {
    result.clear();
    return 0;
  }

  result.resize(ndims);
  dimids.resize(ndims);

  stat = nc_inq_vardimid(ncid, varid, &dimids[0]); check(stat);

  for (int k = 0; k < ndims; ++k) {
    char name[NC_MAX_NAME];
    memset(name, 0, NC_MAX_NAME);

    stat = nc_inq_dimname(ncid, dimids[k], name); check(stat);

    result[k] = name;
  }

  return 0;
}

int PISMNC4File::inq_varnatts(string variable_name, int &result) const {
  int stat, varid = -1;

  if (variable_name == "PISM_GLOBAL") {
    varid = NC_GLOBAL;
  } else {
    stat = nc_inq_varid(ncid, variable_name.c_str(), &varid); check(stat);
  }

  stat = nc_inq_varnatts(ncid, varid, &result); check(stat);

  return 0;
}

int PISMNC4File::inq_varid(string variable_name, bool &exists) const {
  int stat, flag = -1;

  stat = nc_inq_varid(ncid, variable_name.c_str(), &flag);

  if (stat == NC_NOERR)
    flag = 1;
  else
    flag = 0;

  exists = (flag == 1);

  return 0;
}

int PISMNC4File::inq_varname(unsigned int j, string &result) const {
  int stat;
  char varname[NC_MAX_NAME];
  memset(varname, 0, NC_MAX_NAME);

  stat = nc_inq_varname(ncid, j, varname); check(stat);

  result = varname;

  return stat;
}

// att

int PISMNC4File::get_att_double(string variable_name, string att_name, vector<double> &result) const {
  int stat, len, varid = -1;
  size_t attlen;

  // Read the attribute length:
  if (variable_name == "PISM_GLOBAL") {
    varid = NC_GLOBAL;
  } else {
    stat = nc_inq_varid(ncid, variable_name.c_str(), &varid); check(stat);
  }

  stat = nc_inq_attlen(ncid, varid, att_name.c_str(), &attlen);

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
  stat = nc_get_att_double(ncid, varid, att_name.c_str(), &result[0]); check(stat);

  // On error, print a message and stop.
  if (stat != NC_NOERR) {
    fprintf(stderr, "Error reading the %s attribute; (varid %d, NetCDF error %s)",
            att_name.c_str(), varid, nc_strerror(stat));
  }

  return 0;
}


int PISMNC4File::get_att_text(string variable_name, string att_name, string &result) const {
  char *str = NULL;
  int stat, len, varid = -1;

  // Read the attribute length:
  size_t attlen;

  if (variable_name == "PISM_GLOBAL") {
    varid = NC_GLOBAL;
  } else {
    stat = nc_inq_varid(ncid, variable_name.c_str(), &varid); check(stat);
  }

  stat = nc_inq_attlen(ncid, varid, att_name.c_str(), &attlen);
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
  stat = nc_get_att_text(ncid, varid, att_name.c_str(), str);

  // On success, broadcast the string. On error, set str to "".
  if (stat != NC_NOERR) {
    fprintf(stderr, "Error reading the %s attribute; (variable %s, %s)",
            att_name.c_str(), variable_name.c_str(), nc_strerror(stat));

    delete[] str;
    return stat;
  }

  result = str;

  delete[] str;
  return 0;
}

int PISMNC4File::put_att_double(string variable_name, string att_name, PISM_IO_Type xtype, vector<double> &data) const {
  int stat = 0;

  stat = redef(); check(stat);

  int varid = -1;

  if (variable_name == "PISM_GLOBAL") {
    varid = NC_GLOBAL;
  } else {
    stat = nc_inq_varid(ncid, variable_name.c_str(), &varid); check(stat);
  }

  stat = nc_put_att_double(ncid, varid, att_name.c_str(),
                           xtype, data.size(), &data[0]); check(stat);

  return stat;
}

int PISMNC4File::put_att_text(string variable_name, string att_name, string value) const {
  int stat = 0, varid = -1;

  stat = redef(); check(stat);

  if (variable_name == "PISM_GLOBAL") {
    varid = NC_GLOBAL;
  } else {
    stat = nc_inq_varid(ncid, variable_name.c_str(), &varid); check(stat);
  }

  stat = nc_put_att_text(ncid, varid, att_name.c_str(), value.size(), value.c_str()); check(stat);

  return stat;
}

int PISMNC4File::inq_attname(string variable_name, unsigned int n, string &result) const {
  int stat;
  char name[NC_MAX_NAME];
  memset(name, 0, NC_MAX_NAME);

  int varid = -1;

  if (variable_name == "PISM_GLOBAL") {
    varid = NC_GLOBAL;
  } else {
    stat = nc_inq_varid(ncid, variable_name.c_str(), &varid); check(stat);
  }

  stat = nc_inq_attname(ncid, varid, n, name); check(stat);

  result = name;

  return stat;
}

int PISMNC4File::inq_atttype(string variable_name, string att_name, PISM_IO_Type &result) const {
  int stat, varid = -1;
  nc_type tmp;

  if (variable_name == "PISM_GLOBAL") {
    varid = NC_GLOBAL;
  } else {
    stat = nc_inq_varid(ncid, variable_name.c_str(), &varid); check(stat);
  }

  stat = nc_inq_atttype(ncid, varid, att_name.c_str(), &tmp);
  if (stat == NC_ENOTATT) {
    tmp = NC_NAT;
  } else {
    check(stat);
  }

  result = nc_type_to_pism_type(tmp);

  return 0;
}

// misc

int PISMNC4File::set_fill(int fillmode, int &old_modep) const {
  int stat;

  stat = nc_set_fill(ncid, fillmode, &old_modep); check(stat);

  return stat;
}

