// Copyright (C) 2012, 2013, 2014, 2015, 2016, 2017, 2019, 2020 PISM Authors
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

#include "NC_Serial.hh"

// The following is a stupid kludge necessary to make NetCDF 4.x work in
// serial mode in an MPI program:
#ifndef MPI_INCLUDED
#define MPI_INCLUDED 1
#endif
#include <netcdf.h>

#include <cstdio>               // stderr, fprintf

#include "pism/util/pism_utilities.hh" // join
#include "pism/util/error_handling.hh"

#include "pism_type_conversion.hh" // This has to be included *after* netcdf.h.
#include "pism_cdi_type_conversion.hh"

namespace pism {
namespace io {

//! \brief Prints an error message; for debugging.
static void check(const ErrorLocation &where, int return_code) {
  if (return_code != NC_NOERR) {
    throw RuntimeError(where, nc_strerror(return_code));
  }
}

//! call MPI_Abort() if a NetCDF call failed
static void check_and_abort(MPI_Comm com, const ErrorLocation &where, int return_code) {
  if (return_code != NC_NOERR) {
    fprintf(stderr, "%s:%d: %s\n", where.filename, where.line_number, nc_strerror(return_code));
    MPI_Abort(com, -1);
  }
}

NC_Serial::NC_Serial(MPI_Comm c)
  : NCFile(c), m_rank(0) {
  MPI_Comm_rank(m_com, &m_rank);
}

NC_Serial::~NC_Serial() {
  if (m_file_id >= 0) {
    if (m_rank == 0) {
      nc_close(m_file_id);
      fprintf(stderr, "NC_Serial::~NC_Serial: NetCDF file %s is still open\n",
              m_filename.c_str());
    }
    m_file_id = -1;
  }
}

void NC_Serial::set_compression_level_impl(int level) const {
  (void) level;
  // NetCDF-3 does not support compression.
}

// open/create/close
void NC_Serial::open_impl(const std::string &fname, IO_Mode mode) {
  int stat = NC_NOERR;

  int open_mode = mode == PISM_READONLY ? NC_NOWRITE : NC_WRITE;

  if (m_rank == 0) {
    stat = nc_open(fname.c_str(), open_mode, &m_file_id);
  }

  MPI_Barrier(m_com);
  MPI_Bcast(&m_file_id, 1, MPI_INT, 0, m_com);
  MPI_Bcast(&stat, 1, MPI_INT, 0, m_com);

  check(PISM_ERROR_LOCATION, stat);
}

//! \brief Create a NetCDF file.
void NC_Serial::create_impl(const std::string &fname) {
  int stat = NC_NOERR;

  if (m_rank == 0) {
    stat = nc_create(fname.c_str(), NC_CLOBBER | NC_64BIT_OFFSET, &m_file_id);
  }

  MPI_Barrier(m_com);
  MPI_Bcast(&m_file_id, 1, MPI_INT, 0, m_com);
  MPI_Bcast(&stat, 1, MPI_INT, 0, m_com);

  check(PISM_ERROR_LOCATION, stat);
}

//! \brief Close a NetCDF file.
void NC_Serial::close_impl() {
  int stat = NC_NOERR;

  if (m_rank == 0) {
    stat = nc_close(m_file_id);
  }

  m_file_id = -1;

  MPI_Barrier(m_com);
  MPI_Bcast(&stat, 1, MPI_INT, 0, m_com);

  check(PISM_ERROR_LOCATION, stat);
}


void NC_Serial::sync_impl() const {
  int stat = NC_NOERR;

  if (m_rank == 0) {
    stat = nc_sync(m_file_id);
  }

  MPI_Barrier(m_com);
  MPI_Bcast(&stat, 1, MPI_INT, 0, m_com);

  check(PISM_ERROR_LOCATION, stat);
}


//! \brief Exit define mode.
void NC_Serial::enddef_impl() const {
  int stat = NC_NOERR;

  int header_size = 200 * 1024;

  if (m_rank == 0) {
    stat = nc__enddef(m_file_id, header_size, 4, 0, 4);
  }

  MPI_Barrier(m_com);
  MPI_Bcast(&stat, 1, MPI_INT, 0, m_com);

  check(PISM_ERROR_LOCATION, stat);
}

//! \brief Enter define mode.
void NC_Serial::redef_impl() const {
  int stat = NC_NOERR;

  if (m_rank == 0) {
    stat = nc_redef(m_file_id);
  }

  MPI_Barrier(m_com);
  MPI_Bcast(&stat, 1, MPI_INT, 0, m_com);

  check(PISM_ERROR_LOCATION, stat);
}


//! \brief Define a dimension.
void NC_Serial::def_dim_impl(const std::string &name, size_t length) const {
  int stat = NC_NOERR;

  if (m_rank == 0) {
    int dimid;
    stat = nc_def_dim(m_file_id, name.c_str(), length, &dimid);
  }

  MPI_Barrier(m_com);
  MPI_Bcast(&stat, 1, MPI_INT, 0, m_com);

  check(PISM_ERROR_LOCATION, stat);
}

void NC_Serial::inq_dimid_impl(const std::string &dimension_name, bool &exists) const {
  int stat, flag = -1;

  if (m_rank == 0) {
    stat = nc_inq_dimid(m_file_id, dimension_name.c_str(), &flag);

    flag = (stat == NC_NOERR) ? 1 : 0;
  }
  MPI_Barrier(m_com);
  MPI_Bcast(&flag, 1, MPI_INT, 0, m_com);

  exists = (flag == 1);
}


//! \brief Get a dimension length.
void NC_Serial::inq_dimlen_impl(const std::string &dimension_name, unsigned int &result) const {
  int stat = NC_NOERR;

  if (m_rank == 0) {
    int dimid;
    size_t length;

    stat = nc_inq_dimid(m_file_id, dimension_name.c_str(), &dimid);

    if (stat == NC_NOERR) {
      stat = nc_inq_dimlen(m_file_id, dimid, &length);
      result = static_cast<unsigned int>(length);
    }
  }

  MPI_Barrier(m_com);
  MPI_Bcast(&result, 1, MPI_UNSIGNED, 0, m_com);
  MPI_Bcast(&stat,   1, MPI_INT,      0, m_com);

  check(PISM_ERROR_LOCATION, stat);
}

//! \brief Get an unlimited dimension.
void NC_Serial::inq_unlimdim_impl(std::string &result) const {
  int stat = NC_NOERR;
  std::vector<char> dimname(NC_MAX_NAME + 1, 0);

  if (m_rank == 0) {
    int dimid;
    stat = nc_inq_unlimdim(m_file_id, &dimid);

    // nc_inq_unlimdim() sets dimid to -1 if there is no unlimited dimension
    if (stat == NC_NOERR and dimid != -1) {
      stat = nc_inq_dimname(m_file_id, dimid, dimname.data());
    } else {
      // leave dimname empty
      stat = NC_NOERR;
    }
  }

  MPI_Barrier(m_com);

  MPI_Bcast(&stat,   1, MPI_INT, 0, m_com);
  MPI_Bcast(dimname.data(), NC_MAX_NAME, MPI_CHAR, 0, m_com);

  check(PISM_ERROR_LOCATION, stat);

  result = dimname.data();
}

//! \brief Define a variable.
void NC_Serial::def_var_impl(const std::string &name,
                           IO_Type nctype,
                           const std::vector<std::string> &dims) const {
  int stat = NC_NOERR;

  if (m_rank == 0) {
    std::vector<int> dimids;
    int varid;

    for (auto d : dims) {
      int dimid;
      stat = nc_inq_dimid(m_file_id, d.c_str(), &dimid);
      if (stat != NC_NOERR) {
        break;
      }
      dimids.push_back(dimid);
    }

    if (stat == NC_NOERR) {
      stat = nc_def_var(m_file_id, name.c_str(), pism_type_to_nc_type(nctype),
                        static_cast<int>(dims.size()), &dimids[0], &varid);
    }
  }

  MPI_Barrier(m_com);
  MPI_Bcast(&stat, 1, MPI_INT, 0, m_com);

  check(PISM_ERROR_LOCATION, stat);
}

void NC_Serial::get_varm_double_impl(const std::string &variable_name,
                                 const std::vector<unsigned int> &start,
                                 const std::vector<unsigned int> &count,
                                 const std::vector<unsigned int> &imap, double *op) const {
  return this->get_var_double(variable_name,
                              start, count, imap, op, true);
}

void NC_Serial::get_vara_double_impl(const std::string &variable_name,
                                 const std::vector<unsigned int> &start,
                                 const std::vector<unsigned int> &count,
                                 double *op) const {
  std::vector<unsigned int> dummy;
  return this->get_var_double(variable_name,
                              start, count, dummy, op, false);
}

//! \brief Get variable data.
void NC_Serial::get_var_double(const std::string &variable_name,
                            const std::vector<unsigned int> &start_input,
                            const std::vector<unsigned int> &count_input,
                            const std::vector<unsigned int> &imap_input, double *ip,
                            bool transposed) const {
  std::vector<unsigned int> start = start_input;
  std::vector<unsigned int> count = count_input;
  std::vector<unsigned int> imap = imap_input;
  const int start_tag = 1,
    count_tag = 2,
    data_tag =  3,
    imap_tag =  4,
    chunk_size_tag = 5;
  int stat = NC_NOERR, com_size, ndims = static_cast<int>(start.size());
  std::vector<double> processor_0_buffer;
  MPI_Status mpi_stat;
  unsigned int local_chunk_size = 1,
    processor_0_chunk_size = 0;

  if (not transposed) {
    imap.resize(ndims);
  }

  // get the size of the communicator
  MPI_Comm_size(m_com, &com_size);

  // compute the size of a local chunk
  for (int k = 0; k < ndims; ++k) {
    local_chunk_size *= count[k];
  }

  // compute the maximum and send it to processor 0; this is the size of the
  // buffer processor 0 will need
  MPI_Reduce(&local_chunk_size, &processor_0_chunk_size, 1, MPI_UNSIGNED, MPI_MAX, 0, m_com);

  // now we need to send start, count and imap data to processor 0 and receive data
  if (m_rank == 0) {
    // Note: this could be optimized: if processor_0_chunk_size <=
    // max(local_chunk_size) we can avoid allocating this buffer. The inner for
    // loop will have to be re-ordered, though.
    processor_0_buffer.resize(processor_0_chunk_size);

    // MPI calls below require C datatypes (so that we don't have to worry
    // about sizes of size_t and ptrdiff_t), so we make local copies of start,
    // count, and imap to use in the nc_get_varm_double() call.
    std::vector<size_t> nc_start(ndims), nc_count(ndims);
    std::vector<ptrdiff_t> nc_imap(ndims), nc_stride(ndims);
    int varid;

    stat = nc_inq_varid(m_file_id, variable_name.c_str(), &varid);
    check_and_abort(m_com, PISM_ERROR_LOCATION, stat);

    for (int r = 0; r < com_size; ++r) {

      if (r != 0) {
        // Note: start, count, imap, and local_chunk_size on processor zero are
        // used *before* they get overwritten by these calls
        MPI_Recv(&start[0],         ndims, MPI_UNSIGNED, r, start_tag,      m_com, &mpi_stat);
        MPI_Recv(&count[0],         ndims, MPI_UNSIGNED, r, count_tag,      m_com, &mpi_stat);
        MPI_Recv(&imap[0],          ndims, MPI_UNSIGNED, r, imap_tag,       m_com, &mpi_stat);
        MPI_Recv(&local_chunk_size, 1,     MPI_UNSIGNED, r, chunk_size_tag, m_com, &mpi_stat);
      }

      // This for loop uses start, count and imap passed in as arguments when r
      // == 0. For r > 0 they are overwritten by MPI_Recv calls above.
      for (int k = 0; k < ndims; ++k) {
        nc_start[k]  = start[k];
        nc_count[k]  = count[k];
        nc_imap[k]   = imap[k];
        nc_stride[k] = 1;       // fill with ones; this way it works even with
                                // NetCDF versions with a bug affecting the
                                // stride == NULL case.
      }

      if (transposed) {
        stat = nc_get_varm_double(m_file_id, varid, &nc_start[0], &nc_count[0], &nc_stride[0], &nc_imap[0],
                                  &processor_0_buffer[0]);
      } else {
        stat = nc_get_vara_double(m_file_id, varid, &nc_start[0], &nc_count[0],
                                  &processor_0_buffer[0]);
      }
      check_and_abort(m_com, PISM_ERROR_LOCATION, stat);

      if (r != 0) {
        MPI_Send(&processor_0_buffer[0], local_chunk_size, MPI_DOUBLE, r, data_tag, m_com);
      } else {
        for (unsigned int k = 0; k < local_chunk_size; ++k) {
          ip[k] = processor_0_buffer[k];
        }
      }
    } // end of the for loop

  } else {
    MPI_Send(&start[0],          ndims, MPI_UNSIGNED, 0, start_tag,      m_com);
    MPI_Send(&count[0],          ndims, MPI_UNSIGNED, 0, count_tag,      m_com);
    MPI_Send(&imap[0],           ndims, MPI_UNSIGNED, 0, imap_tag,       m_com);
    MPI_Send(&local_chunk_size,  1,     MPI_UNSIGNED, 0, chunk_size_tag, m_com);

    MPI_Recv(ip, local_chunk_size, MPI_DOUBLE, 0, data_tag, m_com, &mpi_stat);
  }
}

void NC_Serial::put_vara_double_impl(const std::string &variable_name,
                                 const std::vector<unsigned int> &start_input,
                                 const std::vector<unsigned int> &count_input,
                                 const double *op) const {
  // make copies of start and count so that we can use them in MPI_Recv() calls below
  std::vector<unsigned int> start = start_input;
  std::vector<unsigned int> count = count_input;
  const int start_tag = 1,
    count_tag = 2,
    data_tag =  3,
    chunk_size_tag = 4;
  int stat = NC_NOERR, com_size = 0, ndims = static_cast<int>(start.size());
  std::vector<double> processor_0_buffer;
  MPI_Status mpi_stat;
  unsigned int local_chunk_size = 1,
    processor_0_chunk_size = 0;

  // get the size of the communicator
  MPI_Comm_size(m_com, &com_size);

  // compute the size of a local chunk
  for (int k = 0; k < ndims; ++k) {
    local_chunk_size *= count[k];
  }

  // compute the maximum and send it to processor 0; this is the size of the
  // buffer processor 0 will need
  MPI_Reduce(&local_chunk_size, &processor_0_chunk_size, 1, MPI_UNSIGNED, MPI_MAX, 0, m_com);

  // now we need to send start and count data to processor 0 and receive data
  if (m_rank == 0) {
    processor_0_buffer.resize(processor_0_chunk_size);

    // MPI calls below require C datatypes (so that we don't have to worry about sizes of
    // size_t and ptrdiff_t), so we make local copies of start and count to use in the
    // nc_get_vara_double() call.
    std::vector<size_t> nc_start(ndims), nc_count(ndims);
    std::vector<ptrdiff_t> nc_stride(ndims);
    int varid;

    stat = nc_inq_varid(m_file_id, variable_name.c_str(), &varid);
    check_and_abort(m_com, PISM_ERROR_LOCATION, stat);

    for (int r = 0; r < com_size; ++r) {

      if (r != 0) {
        // Note: start, count, and local_chunk_size on processor zero are used *before*
        // they get overwritten by these calls
        MPI_Recv(&start[0],         ndims, MPI_UNSIGNED, r, start_tag,      m_com, &mpi_stat);
        MPI_Recv(&count[0],         ndims, MPI_UNSIGNED, r, count_tag,      m_com, &mpi_stat);
        MPI_Recv(&local_chunk_size, 1,     MPI_UNSIGNED, r, chunk_size_tag, m_com, &mpi_stat);

        MPI_Recv(&processor_0_buffer[0], local_chunk_size, MPI_DOUBLE, r, data_tag, m_com, &mpi_stat);
      } else {
        for (unsigned int k = 0; k < local_chunk_size; ++k) {
          processor_0_buffer[k] = op[k];
        }
      }

      // This for loop uses start and count passed in as arguments when r == 0. For r > 0
      // they are overwritten by MPI_Recv calls above.
      for (int k = 0; k < ndims; ++k) {
        nc_start[k]  = start[k];
        nc_count[k]  = count[k];
        nc_stride[k] = 1;       // fill with ones; this way it works even with
                                // NetCDF versions with a bug affecting the
                                // stride == NULL case.
      }

      stat = nc_put_vara_double(m_file_id, varid, &nc_start[0], &nc_count[0],
                                &processor_0_buffer[0]);
      check_and_abort(m_com, PISM_ERROR_LOCATION, stat);
    } // end of the for loop
  } else {
    MPI_Send(&start[0],          ndims, MPI_UNSIGNED, 0, start_tag,      m_com);
    MPI_Send(&count[0],          ndims, MPI_UNSIGNED, 0, count_tag,      m_com);
    MPI_Send(&local_chunk_size,  1,     MPI_UNSIGNED, 0, chunk_size_tag, m_com);

    MPI_Send(const_cast<double*>(op), local_chunk_size, MPI_DOUBLE, 0, data_tag, m_com);
  }
}

//! \brief Get the number of variables.
void NC_Serial::inq_nvars_impl(int &result) const {
  int stat = NC_NOERR;

  if (m_rank == 0) {
    stat = nc_inq_nvars(m_file_id, &result);
  }
  MPI_Barrier(m_com);

  MPI_Bcast(&stat,   1, MPI_INT, 0, m_com);
  check(PISM_ERROR_LOCATION, stat);

  MPI_Bcast(&result, 1, MPI_INT, 0, m_com);
}

//! \brief Get dimensions a variable depends on.
void NC_Serial::inq_vardimid_impl(const std::string &variable_name,
                                std::vector<std::string> &result) const {
  int stat, ndims, varid = -1;
  std::vector<int> dimids;

  if (m_rank == 0) {
    stat = nc_inq_varid(m_file_id, variable_name.c_str(), &varid);

    if (stat == NC_NOERR) {
      stat = nc_inq_varndims(m_file_id, varid, &ndims);
    }
  }

  MPI_Bcast(&stat,   1, MPI_INT, 0, m_com);
  check(PISM_ERROR_LOCATION, stat);

  MPI_Bcast(&ndims, 1, MPI_INT, 0, m_com);

  if (ndims == 0) {
    result.clear();
    return;
  }

  result.resize(ndims);
  dimids.resize(ndims);

  if (m_rank == 0) {
    stat = nc_inq_vardimid(m_file_id, varid, &dimids[0]);
  }

  MPI_Bcast(&stat,   1, MPI_INT, 0, m_com);
  check(PISM_ERROR_LOCATION, stat);

  MPI_Barrier(m_com);

  for (int k = 0; k < ndims; ++k) {
    std::vector<char> name(NC_MAX_NAME + 1, 0);

    if (m_rank == 0) {
      stat = nc_inq_dimname(m_file_id, dimids[k], name.data());
    }

    MPI_Bcast(&stat,   1, MPI_INT, 0, m_com);
    check(PISM_ERROR_LOCATION, stat);

    MPI_Barrier(m_com);
    MPI_Bcast(name.data(), name.size(), MPI_CHAR, 0, m_com);

    result[k] = name.data();
  }
}

//! \brief Get the number of attributes of a variable.
/*!
 * Use "PISM_GLOBAL" as the "variable_name" to get the number of global attributes.
 */
void NC_Serial::inq_varnatts_impl(const std::string &variable_name, int &result) const {
  int stat = NC_NOERR;

  if (m_rank == 0) {
    int varid = get_varid(variable_name);

    if (varid >= NC_GLOBAL) {
      stat = nc_inq_varnatts(m_file_id, varid, &result);
    } else {
      stat = varid;             // LCOV_EXCL_LINE
    }
  }
  MPI_Barrier(m_com);

  MPI_Bcast(&stat, 1, MPI_INT, 0, m_com);
  check(PISM_ERROR_LOCATION, stat);

  MPI_Bcast(&result, 1, MPI_INT, 0, m_com);
}

//! \brief Finds a variable and sets the "exists" flag.
void NC_Serial::inq_varid_impl(const std::string &variable_name, bool &exists) const {
  int stat, flag = -1;

  if (m_rank == 0) {
    stat = nc_inq_varid(m_file_id, variable_name.c_str(), &flag);

    flag = (stat == NC_NOERR) ? 1 : 0;
  }
  MPI_Barrier(m_com);
  MPI_Bcast(&flag, 1, MPI_INT, 0, m_com);

  exists = (flag == 1);
}

void NC_Serial::inq_varname_impl(unsigned int j, std::string &result) const {
  int stat = NC_NOERR;
  std::vector<char> varname(NC_MAX_NAME + 1, 0);

  if (m_rank == 0) {
    stat = nc_inq_varname(m_file_id, j, varname.data());
  }

  MPI_Barrier(m_com);

  MPI_Bcast(&stat,   1, MPI_INT, 0, m_com);
  MPI_Bcast(varname.data(), NC_MAX_NAME, MPI_CHAR, 0, m_com);

  check(PISM_ERROR_LOCATION, stat);

  result = varname.data();
}

//! \brief Gets a double attribute.
/*!
 * Use "PISM_GLOBAL" as the "variable_name" to get the number of global attributes.
 */
void NC_Serial::get_att_double_impl(const std::string &variable_name,
                                  const std::string &att_name,
                                  std::vector<double> &result) const {
  int stat = NC_NOERR, len = 0;

  int varid = get_varid(variable_name);

  // Read and broadcast the attribute length:
  if (m_rank == 0) {
    size_t attlen = 0;

    if (varid >= NC_GLOBAL) {
      stat = nc_inq_attlen(m_file_id, varid, att_name.c_str(), &attlen);
    } else {
      stat = varid;             // LCOV_EXCL_LINE
    }

    if (stat == NC_NOERR) {
      len = static_cast<int>(attlen);
    } else {
      len = 0;
    }
  }
  MPI_Bcast(&len, 1, MPI_INT, 0, m_com);

  if (len == 0) {
    result.clear();
    return;
  }

  result.resize(len);

  // Now read data and broadcast stat to see if we succeeded:
  if (m_rank == 0) {
    stat = nc_get_att_double(m_file_id, varid, att_name.c_str(), &result[0]);
  }
  MPI_Bcast(&stat, 1, MPI_INT, 0, m_com);

  check(PISM_ERROR_LOCATION, stat);

  // Broadcast data
  MPI_Bcast(&result[0], len, MPI_DOUBLE, 0, m_com);
}

// Get a text (character array) attribute on rank 0.
static int get_att_text(int ncid, int varid, const std::string &att_name,
                        std::string &result) {
  int stat = NC_NOERR;

  size_t attlen = 0;
  stat = nc_inq_attlen(ncid, varid, att_name.c_str(), &attlen);
  if (stat != NC_NOERR) {
    result = "";
    return 0;
  }

  std::vector<char> buffer(attlen + 1, 0);
  stat = nc_get_att_text(ncid, varid, att_name.c_str(), &buffer[0]);
  if (stat == NC_NOERR) {
    result = &buffer[0];
  } else {
    result = "";
  }

  return 0;
}

// Get a string attribute on rank 0. In "string array" attributes array elements are concatenated
// using "," as the separator.
static int get_att_string(int ncid, int varid, const std::string &att_name,
                          std::string &result) {
  int stat = NC_NOERR;

  size_t attlen = 0;
  stat = nc_inq_attlen(ncid, varid, att_name.c_str(), &attlen);
  if (stat != NC_NOERR) {
    result = "";
    return 0;
  }

  std::vector<char*> buffer(attlen + 1, 0);
  stat = nc_get_att_string(ncid, varid, att_name.c_str(), &buffer[0]);
  if (stat == NC_NOERR) {
    std::vector<std::string> strings(attlen);
    for (size_t k = 0; k < attlen; ++k) {
      strings[k] = buffer[k];
    }
    result = join(strings, ",");
  } else {
    result = "";
  }
  stat = nc_free_string(attlen, &buffer[0]);

  return stat;
}


//! \brief Gets a text attribute.
/*!
 * Use "PISM_GLOBAL" as the "variable_name" to get the number of global attributes.
 */
void NC_Serial::get_att_text_impl(const std::string &variable_name,
                                const std::string &att_name, std::string &result) const {
  int stat = NC_NOERR;

  // Read and broadcast the attribute length:
  if (m_rank == 0) {

    int varid = get_varid(variable_name);

    if (varid >= NC_GLOBAL) {
      nc_type nctype = NC_NAT;
      stat = nc_inq_atttype(m_file_id, varid, att_name.c_str(), &nctype);

      if (stat == NC_NOERR) {
        switch (nctype) {
        case NC_CHAR:
          stat = pism::io::get_att_text(m_file_id, varid, att_name, result);
          break;
        case NC_STRING:
          stat = pism::io::get_att_string(m_file_id, varid, att_name, result);
          break;
        default:
          result = "";
          stat = NC_NOERR;
        }
      } else if (stat == NC_ENOTATT) {
        result = "";
        stat = NC_NOERR;
      }
    } else {
      stat = varid;             // LCOV_EXCL_LINE
    }
  }
  MPI_Bcast(&stat, 1, MPI_INT, 0, m_com);
  check(PISM_ERROR_LOCATION, stat);

  int len = result.size();
  MPI_Bcast(&len, 1, MPI_INT, 0, m_com);

  result.resize(len);
  MPI_Bcast(&result[0], len, MPI_CHAR, 0, m_com);
}


//! \brief Writes a double attribute.
/*!
 * Use "PISM_GLOBAL" as the "variable_name" to get the number of global attributes.
 */
void NC_Serial::put_att_double_impl(const std::string &variable_name, const std::string &att_name,
                               IO_Type nctype, const std::vector<double> &data) const {
  int stat = NC_NOERR;

  if (m_rank == 0) {
    int varid = get_varid(variable_name);

    if (varid >= NC_GLOBAL) {
      stat = nc_put_att_double(m_file_id, varid, att_name.c_str(),
                               pism_type_to_nc_type(nctype), data.size(), &data[0]);
    } else {
      stat = varid;             // LCOV_EXCL_LINE
    }
  }

  MPI_Barrier(m_com);
  MPI_Bcast(&stat, 1, MPI_INT, 0, m_com);

  check(PISM_ERROR_LOCATION, stat);
}



//! \brief Writes a text attribute.
/*!
 * Use "PISM_GLOBAL" as the "variable_name" to get the number of global attributes.
 */
void NC_Serial::put_att_text_impl(const std::string &variable_name, const std::string &att_name,
                               const std::string &value) const {
  int stat = NC_NOERR;

  if (m_rank == 0) {
    int varid = get_varid(variable_name);

    if (varid >= NC_GLOBAL) {
      stat = nc_put_att_text(m_file_id, varid, att_name.c_str(), value.size(), value.c_str());
    } else {
      stat = varid;             // LCOV_EXCL_LINE
    }
  }

  MPI_Barrier(m_com);
  MPI_Bcast(&stat, 1, MPI_INT, 0, m_com);

  check(PISM_ERROR_LOCATION, stat);
}

//! \brief Gets the name of a numbered attribute.
/*!
 * Use "PISM_GLOBAL" as the "variable_name" to get the number of global attributes.
 */
void NC_Serial::inq_attname_impl(const std::string &variable_name, unsigned int n, std::string &result) const {
  int stat = NC_NOERR;
  std::vector<char> name(NC_MAX_NAME + 1, 0);

  if (m_rank == 0) {
    int varid = get_varid(variable_name);

    if (varid >= NC_GLOBAL) {
      stat = nc_inq_attname(m_file_id, varid, n, name.data()); check(PISM_ERROR_LOCATION, stat);
    } else {
      stat = varid;             // LCOV_EXCL_LINE
    }
  }
  MPI_Barrier(m_com);
  MPI_Bcast(name.data(), NC_MAX_NAME, MPI_CHAR, 0, m_com);
  MPI_Bcast(&stat, 1, MPI_INT, 0, m_com);

  check(PISM_ERROR_LOCATION, stat);

  result = name.data();
}

//! \brief Gets the type of an attribute.
/*!
 * Use "PISM_GLOBAL" as the "variable_name" to get the number of global attributes.
 */
void NC_Serial::inq_atttype_impl(const std::string &variable_name, const std::string &att_name, IO_Type &result) const {
  int stat, tmp;

  if (m_rank == 0) {
    int varid = get_varid(variable_name);

    if (varid >= NC_GLOBAL) {
      // In NetCDF 3.6.x nc_type is an enum; in 4.x it is 'typedef int'.
      nc_type nctype = NC_NAT;
      stat = nc_inq_atttype(m_file_id, varid, att_name.c_str(), &nctype);
      if (stat == NC_ENOTATT) {
        tmp = NC_NAT;
        stat = NC_NOERR;
      } else {
        tmp = static_cast<int>(nctype);
      }
    } else {
      stat = varid;             // LCOV_EXCL_LINE
    }
  }
  MPI_Barrier(m_com);
  MPI_Bcast(&tmp, 1, MPI_INT, 0, m_com);

  MPI_Bcast(&stat, 1, MPI_INT, 0, m_com);
  check(PISM_ERROR_LOCATION, stat);

  result = nc_type_to_pism_type(tmp);
}


//! \brief Sets the fill mode.
void NC_Serial::set_fill_impl(int fillmode, int &old_modep) const {
  int stat = NC_NOERR;

  if (m_rank == 0) {
    stat = nc_set_fill(m_file_id, fillmode, &old_modep);
  }

  MPI_Barrier(m_com);
  MPI_Bcast(&old_modep, 1, MPI_INT, 0, m_com);
  MPI_Bcast(&stat, 1, MPI_INT, 0, m_com);

  check(PISM_ERROR_LOCATION, stat);
}

std::string NC_Serial::get_format() const {
  int format;

  if (m_rank == 0) {
    int stat = nc_inq_format(m_file_id, &format); check(PISM_ERROR_LOCATION, stat);
  }
  MPI_Barrier(m_com);
  MPI_Bcast(&format, 1, MPI_INT, 0, m_com);

  switch(format) {
  case NC_FORMAT_CLASSIC:
  case NC_FORMAT_64BIT_OFFSET:
    return "netcdf3";
  case NC_FORMAT_64BIT_DATA:
    return "cdf5";
  case NC_FORMAT_NETCDF4:
  case NC_FORMAT_NETCDF4_CLASSIC:
  default:
    return "netcdf4";
  }
}

void NC_Serial::del_att_impl(const std::string &variable_name, const std::string &att_name) const {
  int stat = NC_NOERR;

  if (m_rank == 0) {
    int varid = get_varid(variable_name);

    if (varid >= NC_GLOBAL) {
      stat = nc_del_att(m_file_id, varid, att_name.c_str());
    }
  }

  MPI_Barrier(m_com);
  MPI_Bcast(&stat, 1, MPI_INT, 0, m_com);

  check(PISM_ERROR_LOCATION, stat);
}

/*!
 * return the varid corresponding to a variable.
 *
 * If the value returned is NC_GLOBAL or greater, it is a varid, otherwise it is an error
 * code.
 */
int NC_Serial::get_varid(const std::string &variable_name) const {
  if (variable_name == "PISM_GLOBAL") {
    return NC_GLOBAL;
  }

  if (m_rank == 0) {
    int varid = -2;
    int stat = nc_inq_varid(m_file_id, variable_name.c_str(), &varid);

    if (stat == NC_NOERR) {
      return varid;
    } else {
      return stat;
    }
  } else {
    return -2;                  // this value will not be used
  }
}

void NC3File::create_grid_impl(int lengthx, int lengthy) const {
  (void) lengthx;
  (void) lengthy;
}

void NC3File::define_timestep_impl(int tsID) const {
  (void) tsID;
}

void NC3File::write_timestep_impl() const {
}

} // end of namespace io
} // end of namespace pism
