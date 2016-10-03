// Copyright (C) 2013, 2014, 2015, 2016 PISM Authors
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

#include "pismmerge.hh"

#include "base/util/error_handling.hh"

using pism::io::NC4_Serial;

//! \brief Computes the file name corresponding to a patch written by mpi_rank.
std::string patch_filename(const std::string &input, int mpi_rank) {
  char tmp[pism::TEMPORARY_STRING_LENGTH];
  std::string result = input;
  snprintf(tmp, pism::TEMPORARY_STRING_LENGTH, "%04d", mpi_rank);

  std::string::size_type n = result.find("RANK");
  if (n != std::string::npos) {
    snprintf(tmp, pism::TEMPORARY_STRING_LENGTH, "%04d", mpi_rank);
    result.replace(n, 4, tmp);
  } else {
    snprintf(tmp, pism::TEMPORARY_STRING_LENGTH, "-rank%04d", mpi_rank);
    result = pism::pism_filename_add_suffix(result, tmp, "");
  }

  return result;
}

//! \brief Computes the output file name (if not given using -o).
std::string output_filename(const std::string &input, const std::string &var_name) {
  std::string result = input;
  std::string::size_type n = result.find("RANK");
  if (n != std::string::npos) {
    result.replace(n, 4, var_name);
  } else {
    result = pism::pism_filename_add_suffix(result, std::string("-") + var_name, "");
  }

  return result;
}

//! \brief Gets the total number of patches.
int get_quilt_size(const NC4_Serial &input) {
  std::vector<double> tmp;
  input.get_att_double("x_patch", "mpi_size", tmp);
  if (tmp.size() != 1) {
    throw pism::RuntimeError(PISM_ERROR_LOCATION, "x_patch:mpi_size does not exist or has the wrong length.");
  }
  return static_cast<int>(tmp[0]);
}

//! \brief Checks if all input files are present. (We do this before creating
//! the output file to make sure we don't end up bailing in the middle of it.)
void check_input_files(const std::string &filename) {
  NC4_Serial nc(MPI_COMM_SELF, 0);

  nc.open(patch_filename(filename, 0), pism::PISM_READONLY);

  int mpi_size = get_quilt_size(nc);

  nc.close();

  for (int j = 1; j < mpi_size; ++j) {
    nc.open(patch_filename(filename, j), pism::PISM_READONLY);
    nc.close();
  }
}

//! \brief Reads the size of the local patch and its location within the
//! dataset from an input file.
void patch_geometry(const NC4_Serial &input, int &xs, int &ys,
                    unsigned int &xm, unsigned int &ym) {
  std::vector<double> tmp;

  input.get_att_double("x_patch", "patch_offset", tmp);
  if (tmp.size() != 1) {
    throw pism::RuntimeError(PISM_ERROR_LOCATION, "x_patch:patch_offset does not exist or has the wrong length.");
  }
  xs = (int)tmp[0];

  input.get_att_double("y_patch", "patch_offset", tmp);
  if (tmp.size() != 1) {
    throw pism::RuntimeError(PISM_ERROR_LOCATION, "y_patch:patch_offset does not exist or has the wrong length.");
  }
  ys = (int)tmp[0];

  input.inq_dimlen("x_patch", xm);
  input.inq_dimlen("y_patch", ym);
}
