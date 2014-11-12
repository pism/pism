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

#include "PISMNC3File.hh"

// The following is a stupid kludge necessary to make NetCDF 4.x work in
// serial mode in an MPI program:
#ifndef MPI_INCLUDED
#define MPI_INCLUDED 1
#endif
#include <netcdf.h>
#include <cstring>              // memset
#include <cstdio>               // stderr, fprintf

namespace pism {

#include "pism_type_conversion.hh" // This has to be included *after* netcdf.h.

NC3File::NC3File(MPI_Comm c)
  : NCFile(c) {
  MPI_Comm_rank(c, &m_rank);
}

NC3File::~NC3File() {
  if (m_file_id >= 0) {
    if (m_rank == 0) {
      nc_close(m_file_id);
      fprintf(stderr, "NC3File::~NC3File: NetCDF file %s is still open\n",
              m_filename.c_str());
    }
    m_file_id = -1;
  }
}


int NC3File::integer_open_mode(IO_Mode input) const {
  if (input == PISM_READONLY) {
    return NC_NOWRITE;
  } else {
    return NC_WRITE;
  }
}

// open/create/close
int NC3File::open_impl(const std::string &fname, IO_Mode mode) {
  int stat;

  m_filename = fname;

  if (m_rank == 0) {
    int nc_mode = integer_open_mode(mode);
    stat = nc_open(m_filename.c_str(), nc_mode, &m_file_id);
  }

  MPI_Barrier(m_com);
  MPI_Bcast(&m_file_id, 1, MPI_INT, 0, m_com);
  MPI_Bcast(&stat, 1, MPI_INT, 0, m_com);

  m_define_mode = false;

  return stat;
}

//! \brief Create a NetCDF file.
int NC3File::create_impl(const std::string &fname) {
  int stat;

  m_filename = fname;

  if (m_rank == 0) {
    stat = nc_create(m_filename.c_str(), NC_CLOBBER|NC_64BIT_OFFSET, &m_file_id);
  }

  MPI_Barrier(m_com);
  MPI_Bcast(&m_file_id, 1, MPI_INT, 0, m_com);
  MPI_Bcast(&stat, 1, MPI_INT, 0, m_com);

  m_define_mode = true;

  return stat;
}

//! \brief Close a NetCDF file.
int NC3File::close_impl() {
  int stat;

  if (m_rank == 0) {
    stat = nc_close(m_file_id);
    m_file_id = -1;
  }

  MPI_Barrier(m_com);
  MPI_Bcast(&m_file_id, 1, MPI_INT, 0, m_com);
  MPI_Bcast(&stat, 1, MPI_INT, 0, m_com);

  m_filename.clear();

  return stat;
}


//! \brief Exit define mode.
int NC3File::enddef_impl() const {
  int stat;

  if (m_define_mode == false) {
    return 0;
  }

  if (m_rank == 0) {
    //! 50000 (below) means that we allocate ~50Kb for metadata in NetCDF files
    //! created by PISM.
    stat = nc__enddef(m_file_id, 50000, 4, 0, 4); check(stat);
  }

  MPI_Barrier(m_com);
  MPI_Bcast(&stat, 1, MPI_INT, 0, m_com);

  m_define_mode = false;

  return stat;
}

//! \brief Enter define mode.
int NC3File::redef_impl() const {
  int stat;

  if (m_define_mode == true) {
    return 0;
  }

  if (m_rank == 0) {
    stat = nc_redef(m_file_id);
  }

  MPI_Barrier(m_com);
  MPI_Bcast(&stat, 1, MPI_INT, 0, m_com);

  m_define_mode = true;

  return stat;
}


//! \brief Define a dimension.
int NC3File::def_dim_impl(const std::string &name, size_t length) const {
  int stat;

  if (m_rank == 0) {
    int dimid;
    stat = nc_def_dim(m_file_id, name.c_str(), length, &dimid); check(stat);
  }

  MPI_Barrier(m_com);
  MPI_Bcast(&stat, 1, MPI_INT, 0, m_com);

  return stat;
}

int NC3File::inq_dimid_impl(const std::string &dimension_name, bool &exists) const {
  int stat, flag = -1;

  if (m_rank == 0) {
    stat = nc_inq_dimid(m_file_id, dimension_name.c_str(), &flag);

    if (stat == NC_NOERR) {
      flag = 1;
    } else {
      flag = 0;
    }

  }
  MPI_Barrier(m_com);
  MPI_Bcast(&flag, 1, MPI_INT, 0, m_com);

  exists = (flag == 1);

  return 0;
}


//! \brief Get a dimension length.
int NC3File::inq_dimlen_impl(const std::string &dimension_name, unsigned int &result) const {
  int stat;

  if (m_rank == 0) {
    int dimid;
    size_t length;

    stat = nc_inq_dimid(m_file_id, dimension_name.c_str(), &dimid); check(stat);
    if (stat == NC_NOERR) {
      stat = nc_inq_dimlen(m_file_id, dimid, &length); check(stat);
      result = static_cast<unsigned int>(length);
    }
  }

  MPI_Barrier(m_com);
  MPI_Bcast(&result, 1, MPI_UNSIGNED, 0, m_com);
  MPI_Bcast(&stat,   1, MPI_INT,      0, m_com);

  return stat;
}

//! \brief Get an unlimited dimension.
int NC3File::inq_unlimdim_impl(std::string &result) const {
  int stat;
  char dimname[NC_MAX_NAME];
  memset(dimname, 0, NC_MAX_NAME);

  if (m_rank == 0) {
    int dimid;
    stat = nc_inq_unlimdim(m_file_id, &dimid); check(stat);

    if (dimid != -1) {
      stat = nc_inq_dimname(m_file_id, dimid, dimname); check(stat);
    }
  }

  MPI_Barrier(m_com);

  MPI_Bcast(&stat,   1, MPI_INT, 0, m_com);
  MPI_Bcast(dimname, NC_MAX_NAME, MPI_CHAR, 0, m_com);

  result = dimname;

  return stat;
}

int NC3File::inq_dimname_impl(int j, std::string &result) const {
  int stat;
  char dimname[NC_MAX_NAME];
  memset(dimname, 0, NC_MAX_NAME);

  if (m_rank == 0) {
    stat = nc_inq_dimname(m_file_id, j, dimname); check(stat);
  }

  MPI_Barrier(m_com);

  MPI_Bcast(&stat,   1, MPI_INT, 0, m_com);
  MPI_Bcast(dimname, NC_MAX_NAME, MPI_CHAR, 0, m_com);

  result = dimname;

  return stat;
}


int NC3File::inq_ndims_impl(int &result) const {
  int stat;

  if (m_rank == 0) {
    stat = nc_inq_ndims(m_file_id, &result); check(stat);
  }

  MPI_Barrier(m_com);
  MPI_Bcast(&result, 1, MPI_INT, 0, m_com);
  MPI_Bcast(&stat,   1, MPI_INT,      0, m_com);

  return stat;
}


//! \brief Define a variable.
int NC3File::def_var_impl(const std::string &name, IO_Type nctype, const std::vector<std::string> &dims) const {
  int stat;

  if (m_rank == 0) {
    std::vector<int> dimids;
    int varid;

    std::vector<std::string>::const_iterator j;
    for (j = dims.begin(); j != dims.end(); ++j) {
      int dimid;
      stat = nc_inq_dimid(m_file_id, j->c_str(), &dimid); check(stat);
      dimids.push_back(dimid);
    }

    stat = nc_def_var(m_file_id, name.c_str(), pism_type_to_nc_type(nctype),
                      static_cast<int>(dims.size()), &dimids[0], &varid); check(stat);
  }

  MPI_Barrier(m_com);
  MPI_Bcast(&stat,   1, MPI_INT, 0, m_com);

  return stat;
}

int NC3File::get_varm_double_impl(const std::string &variable_name,
                                 const std::vector<unsigned int> &start,
                                 const std::vector<unsigned int> &count,
                                 const std::vector<unsigned int> &imap, double *op) const {
  return this->get_var_double(variable_name,
                              start, count, imap, op, true);
}

int NC3File::get_vara_double_impl(const std::string &variable_name,
                                 const std::vector<unsigned int> &start,
                                 const std::vector<unsigned int> &count,
                                 double *op) const {
  std::vector<unsigned int> dummy;
  return this->get_var_double(variable_name,
                              start, count, dummy, op, false);
}

//! \brief Get variable data.
int NC3File::get_var_double(const std::string &variable_name,
                            const std::vector<unsigned int> &start_input,
                            const std::vector<unsigned int> &count_input,
                            const std::vector<unsigned int> &imap_input, double *ip,
                            bool mapped) const {
  std::vector<unsigned int> start = start_input;
  std::vector<unsigned int> count = count_input;
  std::vector<unsigned int> imap = imap_input;
  const int start_tag = 1,
    count_tag = 2,
    data_tag =  3,
    imap_tag =  4,
    chunk_size_tag = 5;
  int stat = 0, com_size, ndims = static_cast<int>(start.size());
  double *processor_0_buffer = NULL;
  MPI_Status mpi_stat;
  unsigned int local_chunk_size = 1,
    processor_0_chunk_size = 0;

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

  if (mapped == false) {
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
    processor_0_buffer = new double[processor_0_chunk_size];

    // MPI calls below require C datatypes (so that we don't have to worry
    // about sizes of size_t and ptrdiff_t), so we make local copies of start,
    // count, and imap to use in the nc_get_varm_double() call.
    std::vector<size_t> nc_start(ndims), nc_count(ndims);
    std::vector<ptrdiff_t> nc_imap(ndims), nc_stride(ndims);
    int varid;

    stat = nc_inq_varid(m_file_id, variable_name.c_str(), &varid); check(stat);

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

      if (mapped) {
        stat = nc_get_varm_double(m_file_id, varid, &nc_start[0], &nc_count[0], &nc_stride[0], &nc_imap[0],
                                  processor_0_buffer); check(stat);
      } else {
        stat = nc_get_vara_double(m_file_id, varid, &nc_start[0], &nc_count[0],
                                  processor_0_buffer); check(stat);
      }

      if (r != 0) {
        MPI_Send(processor_0_buffer, local_chunk_size, MPI_DOUBLE, r, data_tag, m_com);
      } else {
        for (unsigned int k = 0; k < local_chunk_size; ++k) {
          ip[k] = processor_0_buffer[k];
        }
      }

    } // end of the for loop

    delete[] processor_0_buffer;
  } else {
    MPI_Send(&start[0],          ndims, MPI_UNSIGNED, 0, start_tag,      m_com);
    MPI_Send(&count[0],          ndims, MPI_UNSIGNED, 0, count_tag,      m_com);
    MPI_Send(&imap[0],           ndims, MPI_UNSIGNED, 0, imap_tag,       m_com);
    MPI_Send(&local_chunk_size,  1,     MPI_UNSIGNED, 0, chunk_size_tag, m_com);

    MPI_Recv(ip, local_chunk_size, MPI_DOUBLE, 0, data_tag, m_com, &mpi_stat);
  }

  return stat;
}

int NC3File::put_varm_double_impl(const std::string &variable_name,
                                 const std::vector<unsigned int> &start,
                                 const std::vector<unsigned int> &count,
                                 const std::vector<unsigned int> &imap, const double *op) const {
  return this->put_var_double(variable_name,
                              start, count, imap, op, true);
}

int NC3File::put_vara_double_impl(const std::string &variable_name,
                                 const std::vector<unsigned int> &start,
                                 const std::vector<unsigned int> &count,
                                 const double *op) const {
  std::vector<unsigned int> dummy;
  return this->put_var_double(variable_name,
                              start, count, dummy, op, false);
}


//! \brief Put variable data (mapped).
int NC3File::put_var_double(const std::string &variable_name,
                            const std::vector<unsigned int> &start_input,
                            const std::vector<unsigned int> &count_input,
                            const std::vector<unsigned int> &imap_input, const double *op,
                            bool mapped) const {
  std::vector<unsigned int> start = start_input;
  std::vector<unsigned int> count = count_input;
  std::vector<unsigned int> imap = imap_input;
  const int start_tag = 1,
    count_tag = 2,
    data_tag =  3,
    imap_tag =  4,
    chunk_size_tag = 5;
  int stat = 0, com_size = 0, ndims = static_cast<int>(start.size());
  double *processor_0_buffer = NULL;
  MPI_Status mpi_stat;
  unsigned int local_chunk_size = 1,
    processor_0_chunk_size = 0;

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

  if (mapped == false) {
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
    processor_0_buffer = new double[processor_0_chunk_size];

    // MPI calls below require C datatypes (so that we don't have to worry
    // about sizes of size_t and ptrdiff_t), so we make local copies of start,
    // count, and imap to use in the nc_get_varm_double() call.
    std::vector<size_t> nc_start(ndims), nc_count(ndims);
    std::vector<ptrdiff_t> nc_imap(ndims), nc_stride(ndims);
    int varid;

    stat = nc_inq_varid(m_file_id, variable_name.c_str(), &varid); check(stat);

    for (int r = 0; r < com_size; ++r) {

      if (r != 0) {
        // Note: start, count, imap, and local_chunk_size on processor zero are
        // used *before* they get overwritten by these calls
        MPI_Recv(&start[0],         ndims, MPI_UNSIGNED, r, start_tag,      m_com, &mpi_stat);
        MPI_Recv(&count[0],         ndims, MPI_UNSIGNED, r, count_tag,      m_com, &mpi_stat);
        MPI_Recv(&imap[0],          ndims, MPI_UNSIGNED, r, imap_tag,       m_com, &mpi_stat);
        MPI_Recv(&local_chunk_size, 1,     MPI_UNSIGNED, r, chunk_size_tag, m_com, &mpi_stat);

        MPI_Recv(processor_0_buffer, local_chunk_size, MPI_DOUBLE, r, data_tag, m_com, &mpi_stat);
      } else {
        for (unsigned int k = 0; k < local_chunk_size; ++k) {
          processor_0_buffer[k] = op[k];
        }
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

      if (mapped) {
        stat = nc_put_varm_double(m_file_id, varid, &nc_start[0], &nc_count[0], &nc_stride[0], &nc_imap[0],
                                  processor_0_buffer); check(stat);
      } else {
        stat = nc_put_vara_double(m_file_id, varid, &nc_start[0], &nc_count[0],
                                  processor_0_buffer); check(stat);
      }

      if (stat != NC_NOERR) {
        fprintf(stderr, "NetCDF call nc_put_var?_double failed with return code %d, '%s'\n",
                stat, nc_strerror(stat));
        fprintf(stderr, "while writing '%s' to '%s'\n",
                variable_name.c_str(), m_filename.c_str());

        for (int k = 0; k < ndims; ++k) {
          fprintf(stderr, "start[%d] = %d\n", k, start[k]);
        }

        for (int k = 0; k < ndims; ++k) {
          fprintf(stderr, "count[%d] = %d\n", k, count[k]);
        }

        for (int k = 0; k < ndims; ++k) {
          fprintf(stderr, "imap[%d] = %d\n", k, imap[k]);
        }
      }

    } // end of the for loop

    delete[] processor_0_buffer;
  } else {
    MPI_Send(&start[0],          ndims, MPI_UNSIGNED, 0, start_tag,      m_com);
    MPI_Send(&count[0],          ndims, MPI_UNSIGNED, 0, count_tag,      m_com);
    MPI_Send(&imap[0],           ndims, MPI_UNSIGNED, 0, imap_tag,       m_com);
    MPI_Send(&local_chunk_size,  1,     MPI_UNSIGNED, 0, chunk_size_tag, m_com);

    MPI_Send(const_cast<double*>(op), local_chunk_size, MPI_DOUBLE, 0, data_tag, m_com);
  }

  return stat;
}

//! \brief Get the number of variables.
int NC3File::inq_nvars_impl(int &result) const {
  int stat;

  if (m_rank == 0) {
    stat = nc_inq_nvars(m_file_id, &result); check(stat);
  }
  MPI_Barrier(m_com);
  MPI_Bcast(&result, 1, MPI_INT, 0, m_com);

  return 0;
}

//! \brief Get dimensions a variable depends on.
int NC3File::inq_vardimid_impl(const std::string &variable_name, std::vector<std::string> &result) const {
  int stat, ndims, varid = -1;
  std::vector<int> dimids;

  if (m_rank == 0) {
    stat = nc_inq_varid(m_file_id, variable_name.c_str(), &varid); check(stat);

    stat = nc_inq_varndims(m_file_id, varid, &ndims); check(stat);
  }
  MPI_Bcast(&ndims, 1, MPI_INT, 0, m_com);

  if (ndims == 0) {
    result.clear();
    return 0;
  }

  result.resize(ndims);
  dimids.resize(ndims);

  if (m_rank == 0) {
    stat = nc_inq_vardimid(m_file_id, varid, &dimids[0]); check(stat);
  }

  MPI_Barrier(m_com);

  for (int k = 0; k < ndims; ++k) {
    char name[NC_MAX_NAME];
    memset(name, 0, NC_MAX_NAME);

    if (m_rank == 0) {
      stat = nc_inq_dimname(m_file_id, dimids[k], name); check(stat);
    }

    MPI_Barrier(m_com);
    MPI_Bcast(name, NC_MAX_NAME, MPI_CHAR, 0, m_com);

    result[k] = name;
  }

  return 0;
}

//! \brief Get the number of attributes of a variable.
/*!
 * Use "PISM_GLOBAL" as the "variable_name" to get the number of global attributes.
 */
int NC3File::inq_varnatts_impl(const std::string &variable_name, int &result) const {
  int stat;

  if (m_rank == 0) {
    int varid = -1;

    if (variable_name == "PISM_GLOBAL") {
      varid = NC_GLOBAL;
    } else {
      stat = nc_inq_varid(m_file_id, variable_name.c_str(), &varid); check(stat);
    }

    stat = nc_inq_varnatts(m_file_id, varid, &result); check(stat);
  }
  MPI_Barrier(m_com);
  MPI_Bcast(&result, 1, MPI_INT, 0, m_com);

  return 0;
}

//! \brief Finds a variable and sets the "exists" flag.
int NC3File::inq_varid_impl(const std::string &variable_name, bool &exists) const {
  int stat, flag = -1;

  if (m_rank == 0) {
    stat = nc_inq_varid(m_file_id, variable_name.c_str(), &flag);

    if (stat == NC_NOERR) {
      flag = 1;
    } else {
      flag = 0;
    }

  }
  MPI_Barrier(m_com);
  MPI_Bcast(&flag, 1, MPI_INT, 0, m_com);

  exists = (flag == 1);

  return 0;
}

int NC3File::inq_varname_impl(unsigned int j, std::string &result) const {
  int stat;
  char varname[NC_MAX_NAME];
  memset(varname, 0, NC_MAX_NAME);

  if (m_rank == 0) {
    stat = nc_inq_varname(m_file_id, j, varname); check(stat);
  }

  MPI_Barrier(m_com);

  MPI_Bcast(&stat,   1, MPI_INT, 0, m_com);
  MPI_Bcast(varname, NC_MAX_NAME, MPI_CHAR, 0, m_com);

  result = varname;

  return stat;
}

int NC3File::inq_vartype_impl(const std::string &variable_name, IO_Type &result) const {
  int stat, tmp;

  if (m_rank == 0) {
    nc_type var_type;
    stat = nc_inq_varid(m_file_id, variable_name.c_str(), &tmp); check(stat);
    stat = nc_inq_vartype(m_file_id, tmp, &var_type); check(stat);

    tmp = var_type;
  }

  MPI_Barrier(m_com);

  MPI_Bcast(&stat,   1, MPI_INT, 0, m_com);
  MPI_Bcast(&tmp,    1, MPI_INT, 0, m_com);

  result = nc_type_to_pism_type(tmp);

  return stat;
}



//! \brief Gets a double attribute.
/*!
 * Use "PISM_GLOBAL" as the "variable_name" to get the number of global attributes.
 */
int NC3File::get_att_double_impl(const std::string &variable_name, const std::string &att_name, std::vector<double> &result) const {
  int stat, len, varid = -1;

  // Read and broadcast the attribute length:
  if (m_rank == 0) {
    size_t attlen;

    if (variable_name == "PISM_GLOBAL") {
      varid = NC_GLOBAL;
    } else {
      stat = nc_inq_varid(m_file_id, variable_name.c_str(), &varid); check(stat);
    }

    stat = nc_inq_attlen(m_file_id, varid, att_name.c_str(), &attlen);

    if (stat == NC_NOERR) {
      len = static_cast<int>(attlen);
    } else if (stat == NC_ENOTATT) {
      len = 0;
    } else {
      check(stat);
      len = 0;
    }
  }
  MPI_Bcast(&len, 1, MPI_INT, 0, m_com);

  if (len == 0) {
    result.clear();
    return 0;
  }

  result.resize(len);

  // Now read data and broadcast stat to see if we succeeded:
  if (m_rank == 0) {
    stat = nc_get_att_double(m_file_id, varid, att_name.c_str(), &result[0]); check(stat);
  }
  MPI_Bcast(&stat, 1, MPI_INT, 0, m_com);

  // On success, broadcast the data. On error, stop.
  if (stat == NC_NOERR) {
    MPI_Bcast(&result[0], len, MPI_DOUBLE, 0, m_com);
  } else {
    fprintf(stderr, "Error reading the %s attribute; (varid %d, NetCDF error %s)",
            att_name.c_str(), varid, nc_strerror(stat));
  }

  return 0;
}

//! \brief Gets a text attribute.
/*!
 * Use "PISM_GLOBAL" as the "variable_name" to get the number of global attributes.
 */
int NC3File::get_att_text_impl(const std::string &variable_name, const std::string &att_name, std::string &result) const {
  char *str = NULL;
  int stat, len, varid = -1;

  // Read and broadcast the attribute length:
  if (m_rank == 0) {
    size_t attlen;

    if (variable_name == "PISM_GLOBAL") {
      varid = NC_GLOBAL;
    } else {
      stat = nc_inq_varid(m_file_id, variable_name.c_str(), &varid); check(stat);
    }

    stat = nc_inq_attlen(m_file_id, varid, att_name.c_str(), &attlen);
    if (stat == NC_NOERR) {
      len = static_cast<int>(attlen);
    } else {
      len = 0;
    }
  }
  MPI_Bcast(&len, 1, MPI_INT, 0, m_com);

  // Allocate some memory or clear result and return:
  if (len == 0) {
    result.clear();
    return 0;
  }
  str = new char[len + 1];
  // Zealously clear the string, so that we don't risk moving unitialized bytes
  // over MPI (because Valgrind can't tell the difference between these
  // harmless bytes and potential memory errors)
  memset(str, 0, len + 1);

  // Now read the string and broadcast stat to see if we succeeded:
  if (m_rank == 0) {
    stat = nc_get_att_text(m_file_id, varid, att_name.c_str(), str);
  }
  MPI_Bcast(&stat, 1, MPI_INT, 0, m_com);

  // On success, broadcast the string. On error, set str to "".
  if (stat == NC_NOERR) {
    MPI_Bcast(str, len + 1, MPI_CHAR, 0, m_com);
  } else {
    strcpy(str, "");
  }

  result = str;

  delete[] str;
  return 0;
}

//! \brief Writes a double attribute.
/*!
 * Use "PISM_GLOBAL" as the "variable_name" to get the number of global attributes.
 */
int NC3File::put_att_double_impl(const std::string &variable_name, const std::string &att_name,
                               IO_Type nctype, const std::vector<double> &data) const {

  int stat = 0;

  stat = redef_impl(); check(stat);

  if (m_rank == 0) {
    int varid = -1;

    if (variable_name == "PISM_GLOBAL") {
      varid = NC_GLOBAL;
    } else {
      stat = nc_inq_varid(m_file_id, variable_name.c_str(), &varid); check(stat);
    }

    stat = nc_put_att_double(m_file_id, varid, att_name.c_str(),
                             pism_type_to_nc_type(nctype), data.size(), &data[0]); check(stat);
  }

  MPI_Barrier(m_com);
  MPI_Bcast(&stat, 1, MPI_INT, 0, m_com);

  return stat;
}



//! \brief Writes a text attribute.
/*!
 * Use "PISM_GLOBAL" as the "variable_name" to get the number of global attributes.
 */
int NC3File::put_att_text_impl(const std::string &variable_name, const std::string &att_name, const std::string &value) const {
  int stat = 0;

  stat = redef_impl(); check(stat);

  if (m_rank == 0) {
    int varid = -1;

    if (variable_name == "PISM_GLOBAL") {
      varid = NC_GLOBAL;
    } else {
      stat = nc_inq_varid(m_file_id, variable_name.c_str(), &varid); check(stat);
    }

    stat = nc_put_att_text(m_file_id, varid, att_name.c_str(), value.size(), value.c_str()); check(stat);
  }

  MPI_Barrier(m_com);
  MPI_Bcast(&stat, 1, MPI_INT, 0, m_com);

  return stat;
}

//! \brief Gets the name of a numbered attribute.
/*!
 * Use "PISM_GLOBAL" as the "variable_name" to get the number of global attributes.
 */
int NC3File::inq_attname_impl(const std::string &variable_name, unsigned int n, std::string &result) const {
  int stat;
  char name[NC_MAX_NAME];
  memset(name, 0, NC_MAX_NAME);

  if (m_rank == 0) {
    int varid = -1;

    if (variable_name == "PISM_GLOBAL") {
      varid = NC_GLOBAL;
    } else {
      stat = nc_inq_varid(m_file_id, variable_name.c_str(), &varid); check(stat);
    }

    stat = nc_inq_attname(m_file_id, varid, n, name); check(stat);
  }
  MPI_Barrier(m_com);
  MPI_Bcast(name, NC_MAX_NAME, MPI_CHAR, 0, m_com);
  MPI_Bcast(&stat, 1, MPI_INT, 0, m_com);

  result = name;

  return stat;
}

//! \brief Gets the type of an attribute.
/*!
 * Use "PISM_GLOBAL" as the "variable_name" to get the number of global attributes.
 */
int NC3File::inq_atttype_impl(const std::string &variable_name, const std::string &att_name, IO_Type &result) const {
  int stat, tmp;

  if (m_rank == 0) {
    int varid = -1;

    if (variable_name == "PISM_GLOBAL") {
      varid = NC_GLOBAL;
    } else {
      stat = nc_inq_varid(m_file_id, variable_name.c_str(), &varid); check(stat);
    }

    // In NetCDF 3.6.x nc_type is an enum; in 4.x it is 'typedef int'.
    nc_type nctype = NC_NAT;
    stat = nc_inq_atttype(m_file_id, varid, att_name.c_str(), &nctype);
    if (stat == NC_ENOTATT) {
      tmp = NC_NAT;
    } else {
      tmp = static_cast<int>(nctype);
      check(stat);
    }
  }
  MPI_Barrier(m_com);
  MPI_Bcast(&tmp, 1, MPI_INT, 0, m_com);

  result = nc_type_to_pism_type(tmp);

  return 0;
}


//! \brief Sets the fill mode.
int NC3File::set_fill_impl(int fillmode, int &old_modep) const {
  int stat;

  if (m_rank == 0) {
    stat = nc_set_fill(m_file_id, fillmode, &old_modep); check(stat);
  }

  MPI_Barrier(m_com);
  MPI_Bcast(&old_modep, 1, MPI_INT, 0, m_com);
  MPI_Bcast(&stat, 1, MPI_INT, 0, m_com);

  return stat;
}

std::string NC3File::get_format_impl() const {
  int format;

  if (m_rank == 0) {
    int stat = nc_inq_format(m_file_id, &format); check(stat);
  }
  MPI_Barrier(m_com);
  MPI_Bcast(&format, 1, MPI_INT, 0, m_com);

  switch(format) {
  case NC_FORMAT_CLASSIC:
  case NC_FORMAT_64BIT:
    return "netcdf3";
  case NC_FORMAT_NETCDF4:
  case NC_FORMAT_NETCDF4_CLASSIC:
  default:
    return "netcdf4";
  }
}

} // end of namespace pism
