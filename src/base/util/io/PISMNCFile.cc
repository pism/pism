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
  : com(c) {
  ncid = -1;
  define_mode = false;
  m_xs = m_xm = m_ys = m_ym = -1;
}

NCFile::~NCFile() {
  // empty
}

std::string NCFile::get_filename() const {
  return m_filename;
}

int NCFile::put_att_double(const std::string &variable_name, const std::string &att_name, IO_Type nctype, double value) const {
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

void NCFile::set_local_extent(unsigned int xs, unsigned int xm,
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
int NCFile::move_if_exists(const std::string &file_to_move, int rank_to_use) {
  int stat, rank = 0;
  MPI_Comm_rank(com, &rank);

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
int NCFile::remove_if_exists(const std::string &file_to_remove, int rank_to_use) {
  int stat, rank = 0;
  MPI_Comm_rank(com, &rank);

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

} // end of namespace pism
