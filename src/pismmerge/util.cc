// Copyright (C) 2013 PISM Authors
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

#include "pismmerge.hh"

//! Checks if a NetCDF call succeeded.
void check(int return_code) {
  if (return_code != NC_NOERR) {
    fprintf(stderr, "NC_ERR: %s\n", nc_strerror(return_code));
  }
  assert(return_code == NC_NOERR);
}

//! \brief Computes the file name corresponding to a patch written by mpi_rank.
string patch_filename(string input, int mpi_rank) {
  char tmp[TEMPORARY_STRING_LENGTH];

  snprintf(tmp, TEMPORARY_STRING_LENGTH, "%04d", mpi_rank);

  string::size_type n = input.find("RANK");
  if (n != string::npos) {
    snprintf(tmp, TEMPORARY_STRING_LENGTH, "%04d", mpi_rank);
    input.replace(n, 4, tmp);
  } else {
    snprintf(tmp, TEMPORARY_STRING_LENGTH, "-rank%04d", mpi_rank);
    input = pism_filename_add_suffix(input, tmp, "");
  }

  return input;
}

//! \brief Computes the output file name (if not given using -o).
string output_filename(string input, string var_name) {
  string::size_type n = input.find("RANK");
  if (n != string::npos) {
    input.replace(n, 4, var_name);
  } else {
    input = pism_filename_add_suffix(input, string("-") + var_name, "");
  }

  return input;
}

//! \brief Gets the total number of patches.
int get_quilt_size(PISMNC4_Serial &input, int &mpi_size) {
  int stat;

  vector<double> tmp;
  stat = input.get_att_double("x_patch", "mpi_size", tmp);
  if (stat != 0 || tmp.size() != 1) {
    printf("ERROR: x_patch:mpi_size does not exist or has the wrong length.\n");
    PISMEnd();
  }
  mpi_size = static_cast<int>(tmp[0]);

  return 0;
}

//! \brief Checks if all input files are present. (We do this before creating
//! the output file to make sure we don't end up bailing in the middle of it.)
int check_input_files(string filename) {
  PISMNC4_Serial nc(MPI_COMM_SELF, 0, 0);
  int stat;

  stat = nc.open(patch_filename(filename, 0), PISM_NOWRITE);
  if (stat != 0) {
    printf("ERROR: Cannot open %s!\n", patch_filename(filename, 0).c_str());
    PISMEnd();
  }

  int mpi_size;
  stat = get_quilt_size(nc, mpi_size); check(stat);

  nc.close();

  for (int j = 1; j < mpi_size; ++j) {
    stat = nc.open(patch_filename(filename, j), PISM_NOWRITE);

    if (stat != 0) {
      printf("ERROR: Cannot open %s!\n", patch_filename(filename, j).c_str());
      PISMEnd();
    }

    nc.close();
  }

  return 0;
}

//! \brief Reads the size of the local patch and its location within the
//! dataset from an input file.
int patch_geometry(PISMNC4_Serial &input, int &xs, int &ys,
                   unsigned int &xm, unsigned int &ym) {
  int stat;
  vector<double> tmp;

  stat = input.get_att_double("x_patch", "patch_offset", tmp);
  if (stat != 0 || tmp.size() != 1) {
    printf("ERROR: x_patch:patch_offset does not exist or has the wrong length.\n");
    PISMEnd();
  }
  xs = (int)tmp[0];

  stat = input.get_att_double("y_patch", "patch_offset", tmp);
  if (stat != 0 || tmp.size() != 1) {
    printf("ERROR: y_patch:patch_offset does not exist or has the wrong length.\n");
    PISMEnd();
  }
  ys = (int)tmp[0];

  stat = input.inq_dimlen("x_patch", xm); check(stat);
  stat = input.inq_dimlen("y_patch", ym); check(stat);

  return 0;
}
