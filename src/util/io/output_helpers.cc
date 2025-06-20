/* Copyright (C) 2025 PISM Authors
 *
 * This file is part of PISM.
 *
 * PISM is free software; you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation; either version 3 of the License, or (at your option) any later
 * version.
 *
 * PISM is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License
 * along with PISM; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 */
#include "pism/util/error_handling.hh"
#include "pism/util/io/io_helpers.hh"

namespace pism {
namespace io {

//! \brief Moves the file aside (file.nc -> file.nc~).
/*!
 * Note: only one processor does the renaming.
 */
void move_if_exists(MPI_Comm com, const std::string &file_to_move, int rank_to_use) {
  int stat = 0, rank = 0;
  MPI_Comm_rank(com, &rank);
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
    }
  } // end of "if (rank == rank_to_use)"

  int global_stat = 0;
  MPI_Allreduce(&stat, &global_stat, 1, MPI_INT, MPI_SUM, com);

  if (global_stat != 0) {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION, "PISM ERROR: can't move '%s' to '%s'",
                                  file_to_move.c_str(), backup_filename.c_str());
  }
}

//! \brief Check if a file is present are remove it.
/*!
 * Note: only processor 0 does the job.
 */
void remove_if_exists(MPI_Comm com, const std::string &file_to_remove, int rank_to_use) {
  int stat = 0, rank = 0;
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
    }
  } // end of "if (rank == rank_to_use)"

  int global_stat = 0;
  MPI_Allreduce(&stat, &global_stat, 1, MPI_INT, MPI_SUM, com);

  if (global_stat != 0) {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION, "PISM ERROR: can't remove '%s'",
                                  file_to_remove.c_str());
  }
}

} // namespace io
} // namespace pism
