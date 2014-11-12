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

#include "PISMNC4_Quilt.hh"
#include <assert.h>
#include "pism_const.hh"

// The following is a stupid kludge necessary to make NetCDF 4.x work in
// serial mode in an MPI program:
#ifndef MPI_INCLUDED
#define MPI_INCLUDED 1
#endif
#include <netcdf.h>

namespace pism {

static std::string patch_filename(std::string input, int mpi_rank) {
  char tmp[TEMPORARY_STRING_LENGTH];

  snprintf(tmp, TEMPORARY_STRING_LENGTH, "%04d", mpi_rank);

  std::string::size_type n = input.find("RANK");
  if (n != std::string::npos) {
    snprintf(tmp, TEMPORARY_STRING_LENGTH, "%04d", mpi_rank);
    input.replace(n, 4, tmp);
  } else {
    snprintf(tmp, TEMPORARY_STRING_LENGTH, "-rank%04d", mpi_rank);
    input = pism_filename_add_suffix(input, tmp, "");
  }

  return input;
}

int NC4_Quilt::integer_open_mode(IO_Mode input) const {
  if (input == PISM_READONLY) {
    return NC_NOWRITE;
  } else {
    return NC_WRITE;
  }
}

int NC4_Quilt::open_impl(const std::string &fname, IO_Mode mode) {
  int stat, rank = 0;
  MPI_Comm_rank(m_com, &rank);

  m_filename = patch_filename(fname, rank);

  int nc_mode = integer_open_mode(mode);
  stat = nc_open(m_filename.c_str(), nc_mode, &m_file_id); check(stat);

  m_define_mode = false;

  return global_stat(stat);
}


int NC4_Quilt::create_impl(const std::string &fname) {
  int stat = 0, rank = 0;

  MPI_Comm_rank(m_com, &rank);
  m_filename = patch_filename(fname, rank);

  stat = nc_create(m_filename.c_str(), NC_NETCDF4, &m_file_id); check(stat);

  m_define_mode = true;

  return global_stat(stat);
}

int NC4_Quilt::close_impl() {
  int stat;

  stat = nc_close(m_file_id); check(stat);

  m_file_id = -1;

  m_filename.clear();

  return global_stat(stat);
}

int NC4_Quilt::def_dim_impl(const std::string &name, size_t length) const {
  int stat;
  size_t length_local = 0;

  if (name == "x") {
    assert(m_xm > 0);
    length_local = m_xm;
  } else if (name == "y") {
    assert(m_ym > 0);
    length_local = m_ym;
  }

  if (length_local > 0) {
    stat = this->def_dim_impl(name + suffix, length_local); check(stat);
  }

  stat = NC4File::def_dim_impl(name, length); check(stat);

  return global_stat(stat);
}

int NC4_Quilt::def_var_impl(const std::string &name, IO_Type nctype,
                       const std::vector<std::string> &dims_input) const {
  std::vector<std::string> dims = dims_input;
  int stat = 0, rank = 0;
  MPI_Comm_rank(m_com, &rank);

  if (name == "x" || name == "y") {
    std::vector<std::string> dims_local;
    dims_local.push_back(name + suffix);
    stat = this->def_var_impl(name + suffix, nctype, dims_local); check(stat);

    std::vector<double> buffer(1, name == "x" ? m_xs : m_ys);
    stat = this->put_att_double_impl(name + suffix, "patch_offset", PISM_INT,
                                     buffer); check(stat);

    int size;
    MPI_Comm_size(m_com, &size);
    buffer[0] = rank;
    stat = this->put_att_double_impl(name + suffix, "mpi_rank", PISM_INT, buffer); check(stat);
    buffer[0] = size;
    stat = this->put_att_double_impl(name + suffix, "mpi_size", PISM_INT, buffer); check(stat);
  }

  // Replace "x|y" with "x|y" + suffix (for 2D and 3D variables).
  for (unsigned int j = 0; dims.size() > 1 && j < dims.size(); ++j) {
    if (dims[j] == "x" || dims[j] == "y") {
      dims[j] = dims[j] + suffix;
    }
  }

  stat = NC4File::def_var_impl(name, nctype, dims); check(stat);

  return global_stat(stat);
}

int NC4_Quilt::put_att_double_impl(const std::string &name, const std::string &att_name,
                                      IO_Type xtype, const std::vector<double> &data) const {
  int stat;

  if (name == "x" || name == "y") {
    stat = this->put_att_double_impl(name + suffix, att_name, xtype, data); check(stat);
  }

  stat = NC4File::put_att_double_impl(name, att_name, xtype, data); check(stat);

  return global_stat(stat);
}


int NC4_Quilt::put_att_text_impl(const std::string &name, const std::string &att_name, const std::string &value) const {
  int stat;

  if (name == "x" || name == "y") {
    stat = this->put_att_text_impl(name + suffix, att_name, value); check(stat);
  }

  stat = NC4File::put_att_text_impl(name, att_name, value); check(stat);

  return global_stat(stat);
}


void NC4_Quilt::correct_start_and_count(const std::string &name,
                                            std::vector<unsigned int> &start,
                                            std::vector<unsigned int> &count) const {
  std::vector<std::string> dim_names;

  this->inq_vardimid(name, dim_names);

  for (unsigned int j = 0; j < dim_names.size(); ++j) {
    if (dim_names[j] == "x" + suffix) {
      assert(m_xm > 0);
      start[j] -= m_xs;
      count[j] = PetscMin(count[j], m_xm);
    }

    if (dim_names[j] == "y" + suffix) {
      assert(m_ym > 0);
      start[j] -= m_ys;
      count[j] = PetscMin(count[j], m_ym);
    }
  }

}

int NC4_Quilt::get_put_var_double(const std::string &name,
                                  const std::vector<unsigned int> &start_input,
                                  const std::vector<unsigned int> &count_input,
                                  const std::vector<unsigned int> &imap_input, double *data,
                                  bool get,
                                  bool mapped) const {
  int stat;
  std::vector<unsigned int> start = start_input;
  std::vector<unsigned int> count = count_input;
  std::vector<unsigned int> imap = imap_input;
  if (get) {
    if (!(name == "x" || name == "y")) {
      correct_start_and_count(name, start, count);
    }

    stat = NC4File::get_put_var_double(name, start, count, imap, data, get, mapped); check(stat);
  } else {
    if (name == "x" || name == "y") {
      stat = NC4File::get_put_var_double(name, start, count, imap, data, get, mapped); check(stat);

      double *tmp;
      if (name == "x") {
        assert(m_xm > 0);
        start[0] = 0;
        count[0] = m_xm;
        tmp = &data[m_xs];
      } else {
        assert(m_ym > 0);
        start[0] = 0;
        count[0] = m_ym;
        tmp = &data[m_ys];
      }

      stat = NC4File::get_put_var_double(name + suffix, start, count, imap,
                                         tmp, get, mapped); check(stat);
      return global_stat(stat);
    }

    correct_start_and_count(name, start, count);

    stat = NC4File::get_put_var_double(name, start, count, imap,
                                       data, get, mapped); check(stat);
  }

  return global_stat(stat);
}

int NC4_Quilt::global_stat(int stat) const {
  int tmp;

  MPI_Allreduce(&stat, &tmp, 1, MPI_INT, MPI_SUM, m_com);

  return tmp != 0;
}

int NC4_Quilt::move_if_exists_impl(const std::string &file, int /*rank_to_use*/) {
  int stat = 0, rank = 0;
  MPI_Comm_rank(m_com, &rank);

  stat = NCFile::move_if_exists_impl(patch_filename(file, rank), rank);

  return global_stat(stat);
}

} // end of namespace pism
