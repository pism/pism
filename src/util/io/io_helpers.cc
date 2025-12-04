/* Copyright (C) 2015--2025 PISM Authors
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

#include <mpi.h>
#include <string>
#include <cstdio>

namespace pism {
namespace io {

bool file_exists(MPI_Comm com, const std::string &filename) {
  int file_exists_flag = 0, rank = 0;
  MPI_Comm_rank(com, &rank);

  if (rank == 0) {
    // Check if the file exists:
    if (FILE *f = fopen(filename.c_str(), "r")) {
      file_exists_flag = 1;
      fclose(f);
    } else {
      file_exists_flag = 0;
    }
  }
  MPI_Bcast(&file_exists_flag, 1, MPI_INT, 0, com);

  return file_exists_flag == 1;
}

} // end of namespace io
} // end of namespace pism
