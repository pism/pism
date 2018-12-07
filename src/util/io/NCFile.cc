// Copyright (C) 2012, 2013, 2014, 2015, 2016, 2017 PISM Authors
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

#include "NCFile.hh"

#include <cstdio>               // fprintf, stderr, rename, remove
#include "pism/util/pism_utilities.hh"
#include "pism/util/error_handling.hh"

// The following is a stupid kludge necessary to make NetCDF 4.x work in
// serial mode in an MPI program:
#ifndef MPI_INCLUDED
#define MPI_INCLUDED 1
#endif
#include <netcdf.h>

namespace pism {
namespace io {

NCFile::NCFile(MPI_Comm c)
  : m_com(c) {
  m_file_id = -1;
  m_define_mode = false;
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
  put_att_double(variable_name, att_name, nctype, tmp);
  return 0;
}

//! \brief Prints an error message; for debugging.
void NCFile::check(const ErrorLocation &where, int return_code) const {
  if (return_code != NC_NOERR) {
    throw RuntimeError(where, nc_strerror(return_code));
  }
}

//! \brief Moves the file aside (file.nc -> file.nc~).
/*!
 * Note: only processor 0 does the renaming.
 */
int NCFile::move_if_exists_impl(const std::string &file_to_move, int rank_to_use) {
  int stat = 0, rank = 0;
  MPI_Comm_rank(m_com, &rank);
  std::string backup_filename = file_to_move + "~";

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
      stat = rename(file_to_move.c_str(), backup_filename.c_str());
      if (stat != 0) {
        fprintf(stderr, "PISM ERROR: can't move '%s' to '%s'.\n", file_to_move.c_str(), backup_filename.c_str());
      }

    }

  } // end of "if (rank == rank_to_use)"

  int global_stat = 0;
  MPI_Allreduce(&stat, &global_stat, 1, MPI_INT, MPI_SUM, m_com);

  return global_stat;
}

//! \brief Check if a file is present are remove it.
/*!
 * Note: only processor 0 does the job.
 */
int NCFile::remove_if_exists_impl(const std::string &file_to_remove, int rank_to_use) {
  int stat = 0, rank = 0;
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
        fprintf(stderr, "PISM ERROR: can't remove '%s'.\n", file_to_remove.c_str());
      }
    }
  } // end of "if (rank == rank_to_use)"

  int global_stat = 0;
  MPI_Allreduce(&stat, &global_stat, 1, MPI_INT, MPI_SUM, m_com);

  return global_stat;
}

int NCFile::def_var_chunking_impl(const std::string &name,
                                  std::vector<size_t> &dimensions) const {
  (void) name;
  (void) dimensions;
  // the default implementation does nothing
  return 0;
}


void NCFile::open(const std::string &filename, IO_Mode mode) {
  int stat = this->open_impl(filename, mode); check(PISM_ERROR_LOCATION, stat);
  m_filename = filename;
  m_define_mode = false;
}

void NCFile::create(const std::string &filename) {
  int stat = this->create_impl(filename); check(PISM_ERROR_LOCATION, stat);
  m_filename = filename;
  m_define_mode = true;
}

void NCFile::close() {
  int stat = this->close_impl(); check(PISM_ERROR_LOCATION, stat);
  m_filename.clear();
}

void NCFile::enddef() const {
  if (m_define_mode) {
    int stat = this->enddef_impl(); check(PISM_ERROR_LOCATION, stat);
    m_define_mode = false;
  }
}

void NCFile::redef() const {
  if (not m_define_mode) {
    int stat = this->redef_impl(); check(PISM_ERROR_LOCATION, stat);
    m_define_mode = true;
  }
}

void NCFile::def_dim(const std::string &name, size_t length) const {
  int stat = this->def_dim_impl(name,length); check(PISM_ERROR_LOCATION, stat);
}

void NCFile::inq_dimid(const std::string &dimension_name, bool &exists) const {
  int stat = this->inq_dimid_impl(dimension_name,exists); check(PISM_ERROR_LOCATION, stat);
}

void NCFile::inq_dimlen(const std::string &dimension_name, unsigned int &result) const {
  int stat = this->inq_dimlen_impl(dimension_name,result); check(PISM_ERROR_LOCATION, stat);
}

void NCFile::inq_unlimdim(std::string &result) const {
  int stat = this->inq_unlimdim_impl(result); check(PISM_ERROR_LOCATION, stat);
}

void NCFile::inq_dimname(int j, std::string &result) const {
  int stat = this->inq_dimname_impl(j,result); check(PISM_ERROR_LOCATION, stat);
}

void NCFile::inq_ndims(int &result) const {
  int stat = this->inq_ndims_impl(result); check(PISM_ERROR_LOCATION, stat);
}

void NCFile::def_var(const std::string &name, IO_Type nctype,
                    const std::vector<std::string> &dims) const {
  int stat = this->def_var_impl(name, nctype, dims); check(PISM_ERROR_LOCATION, stat);
}

void NCFile::def_var_chunking(const std::string &name,
                              std::vector<size_t> &dimensions) const {
  int stat = this->def_var_chunking_impl(name, dimensions); check(PISM_ERROR_LOCATION, stat);
}


void NCFile::get_vara_double(const std::string &variable_name,
                            const std::vector<unsigned int> &start,
                            const std::vector<unsigned int> &count,
                            double *ip) const {
  int stat = this->get_vara_double_impl(variable_name, start, count, ip); check(PISM_ERROR_LOCATION, stat);
}

void NCFile::put_vara_double(const std::string &variable_name,
                            const std::vector<unsigned int> &start,
                            const std::vector<unsigned int> &count,
                            const double *op) const {
  int stat = this->put_vara_double_impl(variable_name, start, count, op); check(PISM_ERROR_LOCATION, stat);
}

void NCFile::get_varm_double(const std::string &variable_name,
                            const std::vector<unsigned int> &start,
                            const std::vector<unsigned int> &count,
                            const std::vector<unsigned int> &imap,
                            double *ip) const {
  int stat = this->get_varm_double_impl(variable_name, start, count, imap, ip); check(PISM_ERROR_LOCATION, stat);
}

void NCFile::put_varm_double(const std::string &variable_name,
                            const std::vector<unsigned int> &start,
                            const std::vector<unsigned int> &count,
                            const std::vector<unsigned int> &imap,
                            const double *op) const {
  int stat = this->put_varm_double_impl(variable_name, start, count, imap, op); check(PISM_ERROR_LOCATION, stat);
}

void NCFile::inq_nvars(int &result) const {
  int stat = this->inq_nvars_impl(result); check(PISM_ERROR_LOCATION, stat);
}

void NCFile::inq_vardimid(const std::string &variable_name, std::vector<std::string> &result) const {
  int stat = this->inq_vardimid_impl(variable_name, result); check(PISM_ERROR_LOCATION, stat);
}

void NCFile::inq_varnatts(const std::string &variable_name, int &result) const {
  int stat = this->inq_varnatts_impl(variable_name, result); check(PISM_ERROR_LOCATION, stat);
}

void NCFile::inq_varid(const std::string &variable_name, bool &result) const {
  int stat = this->inq_varid_impl(variable_name, result); check(PISM_ERROR_LOCATION, stat);
}

void NCFile::inq_varname(unsigned int j, std::string &result) const {
  int stat = this->inq_varname_impl(j, result); check(PISM_ERROR_LOCATION, stat);
}

void NCFile::inq_vartype(const std::string &variable_name, IO_Type &result) const {
  int stat = this->inq_vartype_impl(variable_name, result); check(PISM_ERROR_LOCATION, stat);
}

void NCFile::get_att_double(const std::string &variable_name, const std::string &att_name, std::vector<double> &result) const {
  int stat = this->get_att_double_impl(variable_name, att_name, result); check(PISM_ERROR_LOCATION, stat);
}

void NCFile::get_att_text(const std::string &variable_name, const std::string &att_name, std::string &result) const {
  int stat = this->get_att_text_impl(variable_name, att_name, result); check(PISM_ERROR_LOCATION, stat);
}

void NCFile::put_att_double(const std::string &variable_name, const std::string &att_name, IO_Type xtype, const std::vector<double> &data) const {
  int stat = this->put_att_double_impl(variable_name, att_name, xtype, data); check(PISM_ERROR_LOCATION, stat);
}

void NCFile::put_att_double(const std::string &variable_name, const std::string &att_name, IO_Type xtype, double value) const {
  int stat = this->put_att_double_impl(variable_name, att_name, xtype, value); check(PISM_ERROR_LOCATION, stat);
}

void NCFile::put_att_text(const std::string &variable_name, const std::string &att_name, const std::string &value) const {
  int stat = this->put_att_text_impl(variable_name, att_name, value); check(PISM_ERROR_LOCATION, stat);
}

void NCFile::inq_attname(const std::string &variable_name, unsigned int n, std::string &result) const {
  int stat = this->inq_attname_impl(variable_name, n, result); check(PISM_ERROR_LOCATION, stat);
}

void NCFile::inq_atttype(const std::string &variable_name, const std::string &att_name, IO_Type &result) const {
  int stat = this->inq_atttype_impl(variable_name, att_name, result); check(PISM_ERROR_LOCATION, stat);
}

void NCFile::set_fill(int fillmode, int &old_modep) const {
  int stat = this->set_fill_impl(fillmode, old_modep); check(PISM_ERROR_LOCATION, stat);
}

void NCFile::move_if_exists(const std::string &filename, int rank_to_use) {
  int stat = this->move_if_exists_impl(filename, rank_to_use); check(PISM_ERROR_LOCATION, stat);
}

void NCFile::remove_if_exists(const std::string &filename, int rank_to_use) {
  int stat = this->remove_if_exists_impl(filename, rank_to_use); check(PISM_ERROR_LOCATION, stat);
}

} // end of namespace io
} // end of namespace pism
