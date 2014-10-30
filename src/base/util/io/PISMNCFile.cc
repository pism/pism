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

#include "PISMNCFile.hh"

#include <cstdio>               // fprintf, stderr, rename, remove
#include "pism_const.hh"

// The following is a stupid kludge necessary to make NetCDF 4.x work in
// serial mode in an MPI program:
#ifndef MPI_INCLUDED
#define MPI_INCLUDED 1
#endif
#include <netcdf.h>

namespace pism {

NCFile::NCFile(MPI_Comm c)
  : m_com(c) {
  m_file_id = -1;
  m_define_mode = false;
  m_xs = m_xm = m_ys = m_ym = -1;
}

NCFile::~NCFile() {
  // empty
}

std::string NCFile::get_filename() const {
  return m_filename;
}

std::string NCFile::get_format() const {
  return this->get_format_impl();
}

int NCFile::put_att_double_impl(const std::string &variable_name, const std::string &att_name, IO_Type nctype, double value) const {
  std::vector<double> tmp(1);
  tmp[0] = value;
  return put_att_double(variable_name, att_name, nctype, tmp);
}

//! \brief Prints an error message; for debugging.
void NCFile::check(int return_code) const {
  if (return_code != NC_NOERR) {
    fprintf(stderr, "NC_ERR: %s\n", nc_strerror(return_code));
  }
}

void NCFile::set_local_extent_impl(unsigned int xs, unsigned int xm,
                                   unsigned int ys, unsigned int ym) const {
  m_xs = xs;
  m_xm = xm;
  m_ys = ys;
  m_ym = ym;
}

//! \brief Moves the file aside (file.nc -> file.nc~).
/*!
 * Note: only processor 0 does the renaming.
 */
int NCFile::move_if_exists_impl(const std::string &file_to_move, int rank_to_use) {
  int stat, rank = 0;
  MPI_Comm_rank(m_com, &rank);

  if (rank == rank_to_use) {
    bool exists = false;

    // Check if the file exists:
    if (FILE *f = fopen(file_to_move.c_str(), "r")) {
      fclose(f);
      exists = true;
    } else {
      exists = false;
    }

    if (exists) {
      std::string tmp = file_to_move + "~";

      stat = rename(file_to_move.c_str(), tmp.c_str());
      if (stat != 0) {
        printf("PISM ERROR: can't move '%s' to '%s'.\n", file_to_move.c_str(), tmp.c_str());
        return stat;
      }

      if (getVerbosityLevel() >= 2) {
        printf("PISM WARNING: output file '%s' already exists. Moving it to '%s'.\n",
               file_to_move.c_str(), tmp.c_str());
      }

    }

  }

  return 0;
}

//! \brief Check if a file is present are remove it.
/*!
 * Note: only processor 0 does the job.
 */
int NCFile::remove_if_exists_impl(const std::string &file_to_remove, int rank_to_use) {
  int stat, rank = 0;
  MPI_Comm_rank(m_com, &rank);

  if (rank == rank_to_use) {
    bool exists = false;

    // Check if the file exists:
    if (FILE *f = fopen(file_to_remove.c_str(), "r")) {
      fclose(f);
      exists = true;
    } else {
      exists = false;
    }

    if (exists) {
      stat = remove(file_to_remove.c_str());
      if (stat != 0) {
        printf("PISM ERROR: can't remove '%s'.\n", file_to_remove.c_str());
        return stat;
      }

      if (getVerbosityLevel() >= 2) {
        printf("PISM WARNING: output file '%s' already exists. Deleting it...\n",
               file_to_remove.c_str());
      }

    }

  }

  return 0;
}

int NCFile::open(const std::string &filename, IO_Mode mode) {
  return this->open_impl(filename, mode);
}

int NCFile::create(const std::string &filename) {
  return this->create_impl(filename);
}

int NCFile::close() {
  return this->close_impl();
}

int NCFile::enddef() const {
  return this->enddef_impl();
}

int NCFile::redef() const {
  return this->redef_impl();
}

int NCFile::def_dim(const std::string &name, size_t length) const {
  return this->def_dim_impl(name,length);
}

int NCFile::inq_dimid(const std::string &dimension_name, bool &exists) const {
  return this->inq_dimid_impl(dimension_name,exists);
}

int NCFile::inq_dimlen(const std::string &dimension_name, unsigned int &result) const {
  return this->inq_dimlen_impl(dimension_name,result);
}

int NCFile::inq_unlimdim(std::string &result) const {
  return this->inq_unlimdim_impl(result);
}

int NCFile::inq_dimname(int j, std::string &result) const {
  return this->inq_dimname_impl(j,result);
}

int NCFile::inq_ndims(int &result) const {
  return this->inq_ndims_impl(result);
}

int NCFile::def_var(const std::string &name, IO_Type nctype,
                    const std::vector<std::string> &dims) const {
  return this->def_var_impl(name, nctype, dims);
}

int NCFile::get_vara_double(const std::string &variable_name,
                            const std::vector<unsigned int> &start,
                            const std::vector<unsigned int> &count,
                            double *ip) const {
  return this->get_vara_double_impl(variable_name, start, count, ip);
}

int NCFile::put_vara_double(const std::string &variable_name,
                            const std::vector<unsigned int> &start,
                            const std::vector<unsigned int> &count,
                            const double *op) const {
  return this->put_vara_double_impl(variable_name, start, count, op);
}

int NCFile::get_varm_double(const std::string &variable_name,
                            const std::vector<unsigned int> &start,
                            const std::vector<unsigned int> &count,
                            const std::vector<unsigned int> &imap,
                            double *ip) const {
  return this->get_varm_double_impl(variable_name, start, count, imap, ip);
}

int NCFile::put_varm_double(const std::string &variable_name,
                            const std::vector<unsigned int> &start,
                            const std::vector<unsigned int> &count,
                            const std::vector<unsigned int> &imap,
                            const double *op) const {
  return this->put_varm_double_impl(variable_name, start, count, imap, op);
}

int NCFile::inq_nvars(int &result) const {
  return this->inq_nvars_impl(result);
}

int NCFile::inq_vardimid(const std::string &variable_name, std::vector<std::string> &result) const {
  return this->inq_vardimid_impl(variable_name, result);
}

int NCFile::inq_varnatts(const std::string &variable_name, int &result) const {
  return this->inq_varnatts_impl(variable_name, result);
}

int NCFile::inq_varid(const std::string &variable_name, bool &result) const {
  return this->inq_varid_impl(variable_name, result);
}

int NCFile::inq_varname(unsigned int j, std::string &result) const {
  return this->inq_varname_impl(j, result);
}

int NCFile::inq_vartype(const std::string &variable_name, IO_Type &result) const {
  return this->inq_vartype_impl(variable_name, result);
}

int NCFile::get_att_double(const std::string &variable_name, const std::string &att_name, std::vector<double> &result) const {
  return this->get_att_double_impl(variable_name, att_name, result);
}

int NCFile::get_att_text(const std::string &variable_name, const std::string &att_name, std::string &result) const {
  return this->get_att_text_impl(variable_name, att_name, result);
}

int NCFile::put_att_double(const std::string &variable_name, const std::string &att_name, IO_Type xtype, const std::vector<double> &data) const {
  return this->put_att_double_impl(variable_name, att_name, xtype, data);
}

int NCFile::put_att_double(const std::string &variable_name, const std::string &att_name, IO_Type xtype, double value) const {
  return this->put_att_double_impl(variable_name, att_name, xtype, value);
}

int NCFile::put_att_text(const std::string &variable_name, const std::string &att_name, const std::string &value) const {
  return this->put_att_text_impl(variable_name, att_name, value);
}

int NCFile::inq_attname(const std::string &variable_name, unsigned int n, std::string &result) const {
  return this->inq_attname_impl(variable_name, n, result);
}

int NCFile::inq_atttype(const std::string &variable_name, const std::string &att_name, IO_Type &result) const {
  return this->inq_atttype_impl(variable_name, att_name, result);
}

int NCFile::set_fill(int fillmode, int &old_modep) const {
  return this->set_fill_impl(fillmode, old_modep);
}

void NCFile::set_local_extent(unsigned int xs, unsigned int xm,
                              unsigned int ys, unsigned int ym) const {
  return this->set_local_extent_impl(xs, xm, ys, ym);
}

int NCFile::move_if_exists(const std::string &filename, int rank_to_use) {
  return this->move_if_exists_impl(filename, rank_to_use);
}

int NCFile::remove_if_exists(const std::string &filename, int rank_to_use) {
  return this->remove_if_exists_impl(filename, rank_to_use);
}

} // end of namespace pism
