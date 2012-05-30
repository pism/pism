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

#include "PISMNC3File.hh"

// The following is a stupid kludge necessary to make NetCDF 4.x work in
// serial mode in an MPI program:
#ifndef MPI_INCLUDED
#define MPI_INCLUDED 1
#endif
#include <netcdf.h>

#include "pism_type_conversion.hh" // This has to be included *after* netcdf.h.

#include <cstring>              // memset
#include <cstdio>		// stderr, fprintf

PISMNC3File::PISMNC3File(MPI_Comm c, int r)
  : PISMNCFile(c, r) {
}

PISMNC3File::~PISMNC3File() {
  if (ncid >= 0) {
    if (rank == 0) {
      nc_close(ncid);
      fprintf(stderr, "PISMNC3File::~PISMNC3File: NetCDF file %s is still open\n",
              filename.c_str());
    }
    ncid = -1;
  }
}

// open/create/close
int PISMNC3File::open(string fname, int mode) {
  int stat;

  filename = fname;

  if (rank == 0) {
    stat = nc_open(filename.c_str(), mode, &ncid);
  }

  MPI_Barrier(com);
  MPI_Bcast(&ncid, 1, MPI_INT, 0, com);
  MPI_Bcast(&stat, 1, MPI_INT, 0, com);

  define_mode = false;

  return stat;
}

//! \brief Create a NetCDF file.
int PISMNC3File::create(string fname) {
  int stat;

  filename = fname;

  if (rank == 0) {
    stat = nc_create(filename.c_str(), NC_CLOBBER|NC_64BIT_OFFSET, &ncid);
  }

  MPI_Barrier(com);
  MPI_Bcast(&ncid, 1, MPI_INT, 0, com);
  MPI_Bcast(&stat, 1, MPI_INT, 0, com);

  define_mode = true;

  return stat;
}

//! \brief Close a NetCDF file.
int PISMNC3File::close() {
  int stat;

  if (rank == 0) {
    stat = nc_close(ncid);
    ncid = -1;
  }

  MPI_Barrier(com);
  MPI_Bcast(&ncid, 1, MPI_INT, 0, com);
  MPI_Bcast(&stat, 1, MPI_INT, 0, com);

  filename.clear();

  return stat;
}


//! \brief Exit define mode.
int PISMNC3File::enddef() const {
  int stat;

  if (define_mode == false)
    return 0;

  if (rank == 0) {
    //! 50000 (below) means that we allocate ~50Kb for metadata in NetCDF files
    //! created by PISM.
    stat = nc__enddef(ncid, 50000, 4, 0, 4); check(stat);
  }

  MPI_Barrier(com);
  MPI_Bcast(&stat, 1, MPI_INT, 0, com);

  define_mode = false;

  return stat;
}

//! \brief Enter define mode.
int PISMNC3File::redef() const {
  int stat;

  if (define_mode == true)
    return 0;

  if (rank == 0) {
    stat = nc_redef(ncid);
  }

  MPI_Barrier(com);
  MPI_Bcast(&stat, 1, MPI_INT, 0, com);

  define_mode = true;

  return stat;
}


//! \brief Define a dimension.
int PISMNC3File::def_dim(string name, size_t length) const {
  int stat;

  if (rank == 0) {
    int dimid;
    stat = nc_def_dim(ncid, name.c_str(), length, &dimid); check(stat);
  }

  MPI_Barrier(com);
  MPI_Bcast(&stat, 1, MPI_INT, 0, com);

  return stat;
}

int PISMNC3File::inq_dimid(string dimension_name, bool &exists) const {
  int stat, flag = -1;

  if (rank == 0) {
    stat = nc_inq_dimid(ncid, dimension_name.c_str(), &flag);

    if (stat == NC_NOERR)
      flag = 1;
    else
      flag = 0;

  }
  MPI_Barrier(com);
  MPI_Bcast(&flag, 1, MPI_INT, 0, com);

  exists = (flag == 1);

  return 0;
}


//! \brief Get a dimension length.
int PISMNC3File::inq_dimlen(string dimension_name, unsigned int &result) const {
  int stat;

  if (rank == 0) {
    int dimid;
    size_t length;

    stat = nc_inq_dimid(ncid, dimension_name.c_str(), &dimid); check(stat);
    if (stat == NC_NOERR) {
      stat = nc_inq_dimlen(ncid, dimid, &length); check(stat);
      result = static_cast<unsigned int>(length);
    }
  }

  MPI_Barrier(com);
  MPI_Bcast(&result, 1, MPI_UNSIGNED, 0, com);
  MPI_Bcast(&stat,   1, MPI_INT,      0, com);

  return stat;
}

//! \brief Get an unlimited dimension.
int PISMNC3File::inq_unlimdim(string &result) const {
  int stat;
  char dimname[NC_MAX_NAME];
  memset(dimname, 0, NC_MAX_NAME);

  if (rank == 0) {
    int dimid;
    stat = nc_inq_unlimdim(ncid, &dimid); check(stat);

    if (dimid != -1) {
      stat = nc_inq_dimname(ncid, dimid, dimname); check(stat);
    }
  }

  MPI_Barrier(com);

  MPI_Bcast(&stat,   1, MPI_INT, 0, com);
  MPI_Bcast(dimname, NC_MAX_NAME, MPI_CHAR, 0, com);

  result = dimname;

  return stat;
}


//! \brief Define a variable.
int PISMNC3File::def_var(string name, PISM_IO_Type nctype, vector<string> dims) const {
  int stat;

  if (rank == 0) {
    vector<int> dimids;
    int varid;

    vector<string>::iterator j;
    for (j = dims.begin(); j != dims.end(); ++j) {
      int dimid;
      stat = nc_inq_dimid(ncid, j->c_str(), &dimid); check(stat);
      dimids.push_back(dimid);
    }

    stat = nc_def_var(ncid, name.c_str(), pism_type_to_nc_type(nctype),
		      static_cast<int>(dims.size()), &dimids[0], &varid); check(stat);
  }

  MPI_Barrier(com);
  MPI_Bcast(&stat,   1, MPI_INT, 0, com);

  return stat;
}

int PISMNC3File::get_varm_double(string variable_name,
                                 vector<unsigned int> start,
                                 vector<unsigned int> count,
                                 vector<unsigned int> imap, double *op) const {
  return this->get_var_double(variable_name,
                              start, count, imap, op, true);
}

int PISMNC3File::get_vara_double(string variable_name,
                                 vector<unsigned int> start,
                                 vector<unsigned int> count,
                                 double *op) const {
  vector<unsigned int> dummy;
  return this->get_var_double(variable_name,
                              start, count, dummy, op, false);
}

//! \brief Get variable data.
int PISMNC3File::get_var_double(string variable_name,
                                vector<unsigned int> start,
                                vector<unsigned int> count,
                                vector<unsigned int> imap, double *ip,
                                bool mapped) const {
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

  if (mapped == false)
    imap.resize(ndims);

  // get the size of the communicator
  MPI_Comm_size(com, &com_size);

  // compute the size of a local chunk
  for (int k = 0; k < ndims; ++k)
    local_chunk_size *= count[k];

  // compute the maximum and send it to processor 0; this is the size of the
  // buffer processor 0 will need
  MPI_Reduce(&local_chunk_size, &processor_0_chunk_size, 1, MPI_UNSIGNED, MPI_MAX, 0, com);

  // now we need to send start, count and imap data to processor 0 and receive data
  if (rank == 0) {
    // Note: this could be optimized: if processor_0_chunk_size <=
    // max(local_chunk_size) we can avoid allocating this buffer. The inner for
    // loop will have to be re-ordered, though.
    processor_0_buffer = new double[processor_0_chunk_size];

    // MPI calls below require C datatypes (so that we don't have to worry
    // about sizes of size_t and ptrdiff_t), so we make local copies of start,
    // count, and imap to use in the nc_get_varm_double() call.
    vector<size_t> nc_start(ndims), nc_count(ndims);
    vector<ptrdiff_t> nc_imap(ndims), nc_stride(ndims);
    int varid;

    stat = nc_inq_varid(ncid, variable_name.c_str(), &varid); check(stat);

    for (int r = 0; r < com_size; ++r) {

      if (r != 0) {
        // Note: start, count, imap, and local_chunk_size on processor zero are
        // used *before* they get overwritten by these calls
        MPI_Recv(&start[0],         ndims, MPI_UNSIGNED, r, start_tag,      com, &mpi_stat);
        MPI_Recv(&count[0],         ndims, MPI_UNSIGNED, r, count_tag,      com, &mpi_stat);
        MPI_Recv(&imap[0],          ndims, MPI_UNSIGNED, r, imap_tag,       com, &mpi_stat);
        MPI_Recv(&local_chunk_size, 1,     MPI_UNSIGNED, r, chunk_size_tag, com, &mpi_stat);
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
        stat = nc_get_varm_double(ncid, varid, &nc_start[0], &nc_count[0], &nc_stride[0], &nc_imap[0],
                                  processor_0_buffer); check(stat);
      } else {
        stat = nc_get_vara_double(ncid, varid, &nc_start[0], &nc_count[0],
                                  processor_0_buffer); check(stat);
      }

      if (r != 0) {
        MPI_Send(processor_0_buffer, local_chunk_size, MPI_DOUBLE, r, data_tag, com);
      } else {
        for (unsigned int k = 0; k < local_chunk_size; ++k)
          ip[k] = processor_0_buffer[k];
      }

    } // end of the for loop

    delete[] processor_0_buffer;
  } else {
    MPI_Send(&start[0],          ndims, MPI_UNSIGNED, 0, start_tag,      com);
    MPI_Send(&count[0],          ndims, MPI_UNSIGNED, 0, count_tag,      com);
    MPI_Send(&imap[0],           ndims, MPI_UNSIGNED, 0, imap_tag,       com);
    MPI_Send(&local_chunk_size,  1,     MPI_UNSIGNED, 0, chunk_size_tag, com);

    MPI_Recv(ip, local_chunk_size, MPI_DOUBLE, 0, data_tag, com, &mpi_stat);
  }

  return stat;
}

int PISMNC3File::put_varm_double(string variable_name,
                                 vector<unsigned int> start,
                                 vector<unsigned int> count,
                                 vector<unsigned int> imap, const double *op) const {
  return this->put_var_double(variable_name,
                              start, count, imap, op, true);
}

int PISMNC3File::put_vara_double(string variable_name,
                                 vector<unsigned int> start,
                                 vector<unsigned int> count,
                                 const double *op) const {
  vector<unsigned int> dummy;
  return this->put_var_double(variable_name,
                              start, count, dummy, op, false);
}


//! \brief Put variable data (mapped).
int PISMNC3File::put_var_double(string variable_name,
				vector<unsigned int> start,
				vector<unsigned int> count,
				vector<unsigned int> imap, const double *op,
                                bool mapped) const {
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

  if (mapped == false)
    imap.resize(ndims);

  // get the size of the communicator
  MPI_Comm_size(com, &com_size);

  // compute the size of a local chunk
  for (int k = 0; k < ndims; ++k)
    local_chunk_size *= count[k];

  // compute the maximum and send it to processor 0; this is the size of the
  // buffer processor 0 will need
  MPI_Reduce(&local_chunk_size, &processor_0_chunk_size, 1, MPI_UNSIGNED, MPI_MAX, 0, com);

  // now we need to send start, count and imap data to processor 0 and receive data
  if (rank == 0) {
    processor_0_buffer = new double[processor_0_chunk_size];

    // MPI calls below require C datatypes (so that we don't have to worry
    // about sizes of size_t and ptrdiff_t), so we make local copies of start,
    // count, and imap to use in the nc_get_varm_double() call.
    vector<size_t> nc_start(ndims), nc_count(ndims);
    vector<ptrdiff_t> nc_imap(ndims), nc_stride(ndims);
    int varid;

    stat = nc_inq_varid(ncid, variable_name.c_str(), &varid); check(stat);

    for (int r = 0; r < com_size; ++r) {

      if (r != 0) {
        // Note: start, count, imap, and local_chunk_size on processor zero are
        // used *before* they get overwritten by these calls
        MPI_Recv(&start[0],         ndims, MPI_UNSIGNED, r, start_tag,      com, &mpi_stat);
        MPI_Recv(&count[0],         ndims, MPI_UNSIGNED, r, count_tag,      com, &mpi_stat);
        MPI_Recv(&imap[0],          ndims, MPI_UNSIGNED, r, imap_tag,       com, &mpi_stat);
        MPI_Recv(&local_chunk_size, 1,     MPI_UNSIGNED, r, chunk_size_tag, com, &mpi_stat);

        MPI_Recv(processor_0_buffer, local_chunk_size, MPI_DOUBLE, r, data_tag, com, &mpi_stat);
      } else {
        for (unsigned int k = 0; k < local_chunk_size; ++k)
          processor_0_buffer[k] = op[k];
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
        stat = nc_put_varm_double(ncid, varid, &nc_start[0], &nc_count[0], &nc_stride[0], &nc_imap[0],
                                  processor_0_buffer); check(stat);
      } else {
        stat = nc_put_vara_double(ncid, varid, &nc_start[0], &nc_count[0],
                                  processor_0_buffer); check(stat);
      }

      if (stat != NC_NOERR) {
        fprintf(stderr, "NetCDF call nc_put_var?_double failed with return code %d, '%s'\n",
                stat, nc_strerror(stat));
        fprintf(stderr, "while writing '%s' to '%s'\n",
                variable_name.c_str(), filename.c_str());

        for (int k = 0; k < ndims; ++k)
          fprintf(stderr, "start[%d] = %d\n", k, start[k]);

        for (int k = 0; k < ndims; ++k)
          fprintf(stderr, "count[%d] = %d\n", k, count[k]);

        for (int k = 0; k < ndims; ++k)
          fprintf(stderr, "imap[%d] = %d\n", k, imap[k]);
      }

    } // end of the for loop

    delete[] processor_0_buffer;
  } else {
    MPI_Send(&start[0],          ndims, MPI_UNSIGNED, 0, start_tag,      com);
    MPI_Send(&count[0],          ndims, MPI_UNSIGNED, 0, count_tag,      com);
    MPI_Send(&imap[0],           ndims, MPI_UNSIGNED, 0, imap_tag,       com);
    MPI_Send(&local_chunk_size,  1,     MPI_UNSIGNED, 0, chunk_size_tag, com);

    MPI_Send(const_cast<double*>(op), local_chunk_size, MPI_DOUBLE, 0, data_tag, com);
  }

  return stat;
}

//! \brief Get the number of variables.
int PISMNC3File::inq_nvars(int &result) const {
  int stat;

  if (rank == 0) {
    stat = nc_inq_nvars(ncid, &result); check(stat);
  }
  MPI_Barrier(com);
  MPI_Bcast(&result, 1, MPI_INT, 0, com);

  return 0;
}

//! \brief Get dimensions a variable depends on.
int PISMNC3File::inq_vardimid(string variable_name, vector<string> &result) const {
  int stat, ndims, varid = -1;
  vector<int> dimids;

  if (rank == 0) {
    stat = nc_inq_varid(ncid, variable_name.c_str(), &varid); check(stat);

    stat = nc_inq_varndims(ncid, varid, &ndims); check(stat);
  }
  MPI_Bcast(&ndims, 1, MPI_INT, 0, com);

  if (ndims == 0) {
    result.clear();
    return 0;
  }

  result.resize(ndims);
  dimids.resize(ndims);

  if (rank == 0) {
    stat = nc_inq_vardimid(ncid, varid, &dimids[0]); check(stat);
  }

  MPI_Barrier(com);

  for (int k = 0; k < ndims; ++k) {
    char name[NC_MAX_NAME];
    memset(name, 0, NC_MAX_NAME);

    if (rank == 0) {
      stat = nc_inq_dimname(ncid, dimids[k], name); check(stat);
    }

    MPI_Barrier(com);
    MPI_Bcast(name, NC_MAX_NAME, MPI_CHAR, 0, com);

    result[k] = name;
  }

  return 0;
}

//! \brief Get the number of attributes of a variable.
/*!
 * Use "PISM_GLOBAL" as the "variable_name" to get the number of global attributes.
 */
int PISMNC3File::inq_varnatts(string variable_name, int &result) const {
  int stat;

  if (rank == 0) {
    int varid = -1;

    if (variable_name == "PISM_GLOBAL") {
      varid = NC_GLOBAL;
    } else {
      stat = nc_inq_varid(ncid, variable_name.c_str(), &varid); check(stat);
    }

    stat = nc_inq_varnatts(ncid, varid, &result); check(stat);
  }
  MPI_Barrier(com);
  MPI_Bcast(&result, 1, MPI_INT, 0, com);

  return 0;
}

//! \brief Finds a variable and sets the "exists" flag.
int PISMNC3File::inq_varid(string variable_name, bool &exists) const {
  int stat, flag = -1;

  if (rank == 0) {
    stat = nc_inq_varid(ncid, variable_name.c_str(), &flag);

    if (stat == NC_NOERR)
      flag = 1;
    else
      flag = 0;

  }
  MPI_Barrier(com);
  MPI_Bcast(&flag, 1, MPI_INT, 0, com);

  exists = (flag == 1);

  return 0;
}

int PISMNC3File::inq_varname(unsigned int j, string &result) const {
  int stat;
  char varname[NC_MAX_NAME];
  memset(varname, 0, NC_MAX_NAME);

  if (rank == 0) {
    stat = nc_inq_varname(ncid, j, varname); check(stat);
  }

  MPI_Barrier(com);

  MPI_Bcast(&stat,   1, MPI_INT, 0, com);
  MPI_Bcast(varname, NC_MAX_NAME, MPI_CHAR, 0, com);

  result = varname;

  return stat;
}



//! \brief Gets a double attribute.
/*!
 * Use "PISM_GLOBAL" as the "variable_name" to get the number of global attributes.
 */
int PISMNC3File::get_att_double(string variable_name, string att_name, vector<double> &result) const {
  int stat, len, varid = -1;

  // Read and broadcast the attribute length:
  if (rank == 0) {
    size_t attlen;

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
  }
  MPI_Bcast(&len, 1, MPI_INT, 0, com);

  if (len == 0) {
    result.clear();
    return 0;
  }

  result.resize(len);

  // Now read data and broadcast stat to see if we succeeded:
  if (rank == 0) {
    stat = nc_get_att_double(ncid, varid, att_name.c_str(), &result[0]); check(stat);
  }
  MPI_Bcast(&stat, 1, MPI_INT, 0, com);

  // On success, broadcast the data. On error, stop.
  if (stat == NC_NOERR) {
    MPI_Bcast(&result[0], len, MPI_DOUBLE, 0, com);
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
int PISMNC3File::get_att_text(string variable_name, string att_name, string &result) const {
  char *str = NULL;
  int stat, len, varid = -1;

  // Read and broadcast the attribute length:
  if (rank == 0) {
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
  }
  MPI_Bcast(&len, 1, MPI_INT, 0, com);

  // Allocate some memory or set result to NULL and return:
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
  if (rank == 0) {
    stat = nc_get_att_text(ncid, varid, att_name.c_str(), str);
  }
  MPI_Bcast(&stat, 1, MPI_INT, 0, com);

  // On success, broadcast the string. On error, set str to "".
  if (stat == NC_NOERR) {
    MPI_Bcast(str, len + 1, MPI_CHAR, 0, com);
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
int PISMNC3File::put_att_double(string variable_name, string att_name,
                               PISM_IO_Type nctype, vector<double> &data) const {

  int stat = 0;

  stat = redef(); check(stat);

  if (rank == 0) {
    int varid = -1;

    if (variable_name == "PISM_GLOBAL") {
      varid = NC_GLOBAL;
    } else {
      stat = nc_inq_varid(ncid, variable_name.c_str(), &varid); check(stat);
    }

    stat = nc_put_att_double(ncid, varid, att_name.c_str(),
                             pism_type_to_nc_type(nctype), data.size(), &data[0]); check(stat);
  }

  MPI_Barrier(com);
  MPI_Bcast(&stat, 1, MPI_INT, 0, com);

  return stat;
}



//! \brief Writes a text attribute.
/*!
 * Use "PISM_GLOBAL" as the "variable_name" to get the number of global attributes.
 */
int PISMNC3File::put_att_text(string variable_name, string att_name, string value) const {
  int stat = 0;

  stat = redef(); check(stat);

  if (rank == 0) {
    int varid = -1;

    if (variable_name == "PISM_GLOBAL") {
      varid = NC_GLOBAL;
    } else {
      stat = nc_inq_varid(ncid, variable_name.c_str(), &varid); check(stat);
    }

    stat = nc_put_att_text(ncid, varid, att_name.c_str(), value.size(), value.c_str()); check(stat);
  }

  MPI_Barrier(com);
  MPI_Bcast(&stat, 1, MPI_INT, 0, com);

  return stat;
}

//! \brief Gets the name of a numbered attribute.
/*!
 * Use "PISM_GLOBAL" as the "variable_name" to get the number of global attributes.
 */
int PISMNC3File::inq_attname(string variable_name, unsigned int n, string &result) const {
  int stat;
  char name[NC_MAX_NAME];
  memset(name, 0, NC_MAX_NAME);

  if (rank == 0) {
    int varid = -1;

    if (variable_name == "PISM_GLOBAL") {
      varid = NC_GLOBAL;
    } else {
      stat = nc_inq_varid(ncid, variable_name.c_str(), &varid); check(stat);
    }

    stat = nc_inq_attname(ncid, varid, n, name); check(stat);
  }
  MPI_Barrier(com);
  MPI_Bcast(name, NC_MAX_NAME, MPI_CHAR, 0, com);
  MPI_Bcast(&stat, 1, MPI_INT, 0, com);

  result = name;

  return stat;
}

//! \brief Gets the type of an attribute.
/*!
 * Use "PISM_GLOBAL" as the "variable_name" to get the number of global attributes.
 */
int PISMNC3File::inq_atttype(string variable_name, string att_name, PISM_IO_Type &result) const {
  int stat, tmp;

  if (rank == 0) {
    int varid = -1;

    if (variable_name == "PISM_GLOBAL") {
      varid = NC_GLOBAL;
    } else {
      stat = nc_inq_varid(ncid, variable_name.c_str(), &varid); check(stat);
    }

    // In NetCDF 3.6.x nc_type is an enum; in 4.x it is 'typedef int'.
    nc_type nctype = NC_NAT;
    stat = nc_inq_atttype(ncid, varid, att_name.c_str(), &nctype);
    if (stat == NC_ENOTATT) {
      tmp = NC_NAT;
    } else {
      tmp = static_cast<int>(nctype);
      check(stat);
    }
  }
  MPI_Barrier(com);
  MPI_Bcast(&tmp, 1, MPI_INT, 0, com);

  result = nc_type_to_pism_type(static_cast<nc_type>(tmp));

  return 0;
}


//! \brief Sets the fill mode.
int PISMNC3File::set_fill(int fillmode, int &old_modep) const {
  int stat;

  if (rank == 0) {
    stat = nc_set_fill(ncid, fillmode, &old_modep); check(stat);
  }

  MPI_Barrier(com);
  MPI_Bcast(&old_modep, 1, MPI_INT, 0, com);
  MPI_Bcast(&stat, 1, MPI_INT, 0, com);

  return stat;
}
