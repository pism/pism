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
#include <netcdf_par.h>

PISMNC4File::PISMNC4File(MPI_Comm c, int r)
  : PISMNCFile(c, r) {
  // empty
}

PISMNC4File::~PISMNC4File() {
  // empty
}

// open/create/close

int PISMNC4File::open(string fname, int mode) {

  filename = fname;

  int ierr = nc_open_par(filename.c_str(), mode, &ncid); check(ierr);

  define_mode = false;

  return ierr;
}

int PISMNC4File::create(string fname, int mode) {

  filename = fname;

  int ierr = nc_create_par(filename.c_str(), mode, &ncid); check(ierr);

  define_mode = true;

  return ierr;
}

int PISMNC4File::close() {
  int ierr;

  ierr = nc_close(ncid); check(ierr);

  ncid = -1;

  filename.clear();

  return ierr;
}

// redef/enddef
int PISMNC4File::enddef() const {

  if (define_mode == false)
    return 0;

  int ierr = nc_enddef(ncid); check(ierr);

  define_mode = false;

  return ierr;
}

int PISMNC4File::redef() const {

  if (define_mode == true)
    return 0;

  int ierr = nc_redef(ncid); check(ierr);

  define_mode = true;

  return ierr;
}

// dim
int PISMNC4File::def_dim(string name, size_t length) const {
  int dimid = 0, ierr;

  ierr = nc_def_dim(ncid, name.c_str(), length, &dimid); check(ierr);

  return ierr;
}

int PISMNC4File::inq_dimid(string dimension_name, bool &exists) const {
  int tmp, ierr;

  ierr = nc_inq_dimid(ncid, dimension_name.c_str(), &tmp);

  if (ierr == NC_NOERR) {
    exists = true;
  } else {
    exists = false;
  }

  return 0;
}

int PISMNC4File::inq_dimlen(string dimension_name, unsigned int &result) const {
  int ierr;
  size_t len;

  ierr = nc_inq_dimlen(ncid, dimension_name.c_str(), &len); check(ierr);

  result = static_cast<unsigned int>(len);

  return ierr;
}

int PISMNC4File::inq_unlimdim(string &result) const {
  int ierr, dimid;
  char dimname[NC_MAX_NAME];

  ierr = nc_inq_unlimdim(ncid, &dimid); check(ierr);

  if (dimid == -1) {
    result.clear();
  } else {
    ierr = nc_inq_dimname(ncid, dimid, dimname); check(ierr);

    result = dimname;
  }

  return ierr;
}

// var
int PISMNC4File::def_var(string name, nc_type nctype, vector<string> dims) const {
  vector<int> dimids;
  int ierr, varid;

  vector<string>::iterator j;
  for (j = dims.begin(); j != dims.end(); ++j) {
    int dimid;
    ierr = nc_inq_dimid(ncid, j->c_str(), &dimid); check(ierr);
    dimids.push_back(dimid);
  }

  ierr = nc_def_var(ncid, name.c_str(), nctype,
                    static_cast<int>(dims.size()), &dimids[0], &varid); check(ierr);

  return ierr;
}

int PISMNC4File::get_varm_double(string variable_name,
                                 vector<unsigned int> start,
                                 vector<unsigned int> count,
                                 vector<unsigned int> imap, double *ip) const {
  int ierr, varid, ndims = static_cast<int>(start.size());

  vector<size_t> nc_start(ndims), nc_count(ndims);
  vector<ptrdiff_t> nc_imap(ndims), nc_stride(ndims);

  ierr = nc_inq_varid(ncid, variable_name.c_str(), &varid); check(ierr);

  for (int j = 0; j < ndims; ++j) {
    nc_start[j] = start[j];
    nc_count[j] = count[j];
    nc_imap[j]  = imap[j];
    nc_stride[j] = 1;
  }

  ierr = nc_get_varm_double(ncid, varid,
                            &nc_start[0], &nc_count[0], &nc_stride[0], &nc_imap[0],
                            ip); check(ierr);

  return ierr;
}

int PISMNC4File::put_varm_double(string variable_name,
                                 vector<unsigned int> start,
                                 vector<unsigned int> count,
                                 vector<unsigned int> imap, const double *op) const {
  int ierr, varid, ndims = static_cast<int>(start.size());

  vector<size_t> nc_start(ndims), nc_count(ndims);
  vector<ptrdiff_t> nc_imap(ndims), nc_stride(ndims);

  ierr = nc_inq_varid(ncid, variable_name.c_str(), &varid); check(ierr);

  for (int j = 0; j < ndims; ++j) {
    nc_start[j] = start[j];
    nc_count[j] = count[j];
    nc_imap[j]  = imap[j];
    nc_stride[j] = 1;
  }

  ierr = nc_put_varm_double(ncid, varid,
                            &nc_start[0], &nc_count[0], &nc_stride[0], &nc_imap[0],
                            op); check(ierr);

  return ierr;
}

int PISMNC4File::inq_nvars(int &result) const {
  int ierr;

  ierr = nc_inq_nvars(ncid, &result); check(ierr);

  return ierr;
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

    if (rank == 0) {
      stat = nc_inq_dimname(ncid, dimids[k], name); check(stat);
    }

    result[k] = name;
  }

  return 0;
}

int PISMNC4File::inq_varnatts(string variable_name, int &result) const {
  int stat, varid = -1;

  if (variable_name == "NC_GLOBAL") {
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
  if (variable_name == "NC_GLOBAL") {
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

  if (variable_name == "NC_GLOBAL") {
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

int PISMNC4File::put_att_double(string variable_name, string att_name, nc_type xtype, vector<double> &data) const {
  int stat = 0;

  stat = redef(); check(stat);

  int varid = -1;

  if (variable_name == "NC_GLOBAL") {
    varid = NC_GLOBAL;
  } else {
    stat = nc_inq_varid(ncid, variable_name.c_str(), &varid); check(stat);
  }

  stat = nc_put_att_double(ncid, varid, att_name.c_str(),
                           nctype, data.size(), &data[0]); check(stat);

  return stat;
}

int PISMNC4File::put_att_text(string variable_name, string att_name, string value) const {
  int stat = 0, varid = -1;

  stat = redef(); check(stat);

  if (variable_name == "NC_GLOBAL") {
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

  if (variable_name == "NC_GLOBAL") {
    varid = NC_GLOBAL;
  } else {
    stat = nc_inq_varid(ncid, variable_name.c_str(), &varid); check(stat);
  }

  stat = nc_inq_attname(ncid, varid, n, name); check(stat);

  result = name;

  return stat;
}

int PISMNC4File::inq_atttype(string variable_name, string att_name, nc_type &result) const {
  int stat, tmp, varid = -1;

  if (variable_name == "NC_GLOBAL") {
    varid = NC_GLOBAL;
  } else {
    stat = nc_inq_varid(ncid, variable_name.c_str(), &varid); check(stat);
  }

  stat = nc_inq_atttype(ncid, varid, att_name.c_str(), &result); check(stat);

  return 0;
}

// misc

int PISMNC4File::set_fill(int fillmode, int &old_modep) const {
  int stat;

  stat = nc_set_fill(ncid, fillmode, &old_modep); check(stat);

  return stat;
}

