// Copyright (C) 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019, 2020, 2023 PISM Authors
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
#include "pism/util/Grid.hh"

// The following is a stupid kludge necessary to make NetCDF 4.x work in
// serial mode in an MPI program:
#ifndef MPI_INCLUDED
#define MPI_INCLUDED 1
#endif
#include <netcdf.h>

namespace pism {
namespace io {

NCFile::NCFile(MPI_Comm c)
  : m_com(c), m_file_id(-1), m_define_mode(false) {
}

std::string NCFile::filename() const {
  return m_filename;
}

void NCFile::set_compression_level(int level) const {
  set_compression_level_impl(level);
}

void NCFile::set_compression_level_impl(int level) const {
  (void) level;
  // the default implementation does nothing
}

void NCFile::def_var_chunking_impl(const std::string &name,
                                   std::vector<size_t> &dimensions) const {
  (void) name;
  (void) dimensions;
  // the default implementation does nothing
}


void NCFile::open(const std::string &filename, io::Mode mode) {
  this->open_impl(filename, mode);
  m_filename = filename;
  m_define_mode = false;
}

void NCFile::create(const std::string &filename) {
  this->create_impl(filename);
  m_filename = filename;
  m_define_mode = true;
}

void NCFile::sync() const {
  enddef();
  this->sync_impl();
}

void NCFile::close() {
  this->close_impl();
  m_filename.clear();
  m_file_id = -1;
}

void NCFile::enddef() const {
  if (m_define_mode) {
    this->enddef_impl();
    m_define_mode = false;
  }
}

void NCFile::redef() const {
  if (not m_define_mode) {
    this->redef_impl();
    m_define_mode = true;
  }
}

void NCFile::def_dim(const std::string &name, size_t length) const {
  redef();
  this->def_dim_impl(name, length);
}

void NCFile::inq_dimid(const std::string &dimension_name, bool &exists) const {
  this->inq_dimid_impl(dimension_name,exists);
}

void NCFile::inq_dimlen(const std::string &dimension_name, unsigned int &result) const {
  this->inq_dimlen_impl(dimension_name,result);
}

void NCFile::inq_unlimdim(std::string &result) const {
  this->inq_unlimdim_impl(result);
}

void NCFile::def_var(const std::string &name, io::Type nctype,
                    const std::vector<std::string> &dims) const {
  redef();
  this->def_var_impl(name, nctype, dims);
}

void NCFile::def_var_chunking(const std::string &name,
                              std::vector<size_t> &dimensions) const {
  this->def_var_chunking_impl(name, dimensions);
}


void NCFile::get_vara_double(const std::string &variable_name,
                            const std::vector<unsigned int> &start,
                            const std::vector<unsigned int> &count,
                            double *ip) const {
#if (Pism_DEBUG==1)
  if (start.size() != count.size()) {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                  "start and count arrays have to have the same size");
  }
#endif

  enddef();
  this->get_vara_double_impl(variable_name, start, count, ip);
}

void NCFile::put_vara_double(const std::string &variable_name,
                            const std::vector<unsigned int> &start,
                            const std::vector<unsigned int> &count,
                            const double *op) const {
#if (Pism_DEBUG==1)
  if (start.size() != count.size()) {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                  "start and count arrays have to have the same size");
  }
#endif

  enddef();
  this->put_vara_double_impl(variable_name, start, count, op);
}


void NCFile::write_darray(const std::string &variable_name,
                          const Grid &grid,
                          unsigned int z_count,
                          bool time_dependent,
                          unsigned int record,
                          const double *input) {
  enddef();
  this->write_darray_impl(variable_name, grid, z_count, time_dependent, record, input);
}

/*!
 * The default implementation computes start and count and calls put_vara_double()
 */
void NCFile::write_darray_impl(const std::string &variable_name,
                               const Grid &grid,
                               unsigned int z_count,
                               bool time_dependent,
                               unsigned int record,
                               const double *input) {
  std::vector<std::string> dims;
  this->inq_vardimid(variable_name, dims);

  std::vector<unsigned int> start, count;

  // time
  if (time_dependent) {
    start.push_back(record);
    count.push_back(1);
  }

  // y
  start.push_back(grid.ys());
  count.push_back(grid.ym());

  // x
  start.push_back(grid.xs());
  count.push_back(grid.xm());

  // z (these are not used when writing 2D fields)
  start.push_back(0);
  count.push_back(z_count);

  this->put_vara_double(variable_name, start, count, input);
}


void NCFile::get_varm_double(const std::string &variable_name,
                            const std::vector<unsigned int> &start,
                            const std::vector<unsigned int> &count,
                            const std::vector<unsigned int> &imap,
                            double *ip) const {

#if (Pism_DEBUG==1)
  if (start.size() != count.size() or
      start.size() != imap.size()) {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                  "start, count and imap arrays have to have the same size");
  }
#endif

  enddef();
  this->get_varm_double_impl(variable_name, start, count, imap, ip);
}

void NCFile::inq_nvars(int &result) const {
  this->inq_nvars_impl(result);
}

void NCFile::inq_vardimid(const std::string &variable_name, std::vector<std::string> &result) const {
  this->inq_vardimid_impl(variable_name, result);
}

void NCFile::inq_varnatts(const std::string &variable_name, int &result) const {
  this->inq_varnatts_impl(variable_name, result);
}

void NCFile::inq_varid(const std::string &variable_name, bool &result) const {
  this->inq_varid_impl(variable_name, result);
}

void NCFile::inq_varname(unsigned int j, std::string &result) const {
  this->inq_varname_impl(j, result);
}

void NCFile::get_att_double(const std::string &variable_name,
                            const std::string &att_name,
                            std::vector<double> &result) const {
  this->get_att_double_impl(variable_name, att_name, result);
}

void NCFile::get_att_text(const std::string &variable_name,
                          const std::string &att_name,
                          std::string &result) const {
  this->get_att_text_impl(variable_name, att_name, result);
}

void NCFile::put_att_double(const std::string &variable_name,
                            const std::string &att_name,
                            io::Type xtype,
                            const std::vector<double> &data) const {
  this->put_att_double_impl(variable_name, att_name, xtype, data);
}

void NCFile::put_att_text(const std::string &variable_name,
                          const std::string &att_name,
                          const std::string &value) const {
  this->put_att_text_impl(variable_name, att_name, value);
}

void NCFile::inq_attname(const std::string &variable_name,
                         unsigned int n,
                         std::string &result) const {
  this->inq_attname_impl(variable_name, n, result);
}

void NCFile::inq_atttype(const std::string &variable_name,
                         const std::string &att_name,
                         io::Type &result) const {
  this->inq_atttype_impl(variable_name, att_name, result);
}

void NCFile::set_fill(int fillmode, int &old_modep) const {
  redef();
  this->set_fill_impl(fillmode, old_modep);
}

void NCFile::del_att(const std::string &variable_name, const std::string &att_name) const {
  this->del_att_impl(variable_name, att_name);
}

} // end of namespace io
} // end of namespace pism
