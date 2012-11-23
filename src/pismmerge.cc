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

#include "pism_const.hh"
#include "pism_options.hh"
#include "PISMNC4_Serial.hh"
#include <string>
#include <cstdio>
#include <stdlib.h>
#include <map>
#include <assert.h>

// The following is a stupid kludge necessary to make NetCDF 4.x work in
// serial mode in an MPI program:
#ifndef MPI_INCLUDED
#define MPI_INCLUDED 1
#endif
#include <netcdf.h>

using namespace std;

static char help[] =
  "Tool for merging PISM output files produced using '-o_format quilt'.\n";

//! Checks if a NetCDF call succeeded.
void check(int return_code) {
  assert(return_code == NC_NOERR);
  if (return_code != NC_NOERR) {
    fprintf(stderr, "NC_ERR: %s\n", nc_strerror(return_code));
  }
}

//! \brief Computes the file name corresponding to a patch written by mpi_rank.
static string patch_filename(string input, int mpi_rank) {
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
string output_filename(string input) {
  string::size_type n = input.find("RANK");
  if (n != string::npos) {
    input.replace(n, 4, "ALL");
  } else {
    input = pism_filename_add_suffix(input, "-ALL", "");
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
static int check_input_files(string filename) {
  PISMNC4_Serial nc(MPI_COMM_SELF, 0);
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

//! \brief Copies variable attributes.
int copy_attributes(PISMNC4_Serial &input, PISMNC4_Serial &output, string var_name) {
  int stat;
  int n_attrs;

  stat = input.inq_varnatts(var_name, n_attrs); check(stat);

  for (int j = 0; j < n_attrs; ++j) {
    string att_name;
    PISM_IO_Type att_type;

    stat = input.inq_attname(var_name, j, att_name); check(stat);

    stat = input.inq_atttype(var_name, att_name, att_type); check(stat);

    if (att_type == PISM_CHAR) {
      string tmp;

      stat = input.get_att_text(var_name, att_name, tmp); check(stat);

      stat = output.put_att_text(var_name, att_name, tmp); check(stat);
    } else {
      vector<double> tmp;

      stat = input.get_att_double(var_name, att_name, tmp); check(stat);

      stat = output.put_att_double(var_name, att_name, att_type, tmp); check(stat);
    }
  }

  return 0;
}

//! \brief Defines dimensions in the output file.
int define_dimensions(PISMNC4_Serial &input, PISMNC4_Serial &output) {
  int stat, n_dims;

  stat = input.inq_ndims(n_dims); check(stat);

  for (int j = 0; j < n_dims; ++j) {
    string dim_name;
    unsigned int dim_len;

    stat = input.inq_dimname(j, dim_name); check(stat);

    if (dim_name == "x_patch" || dim_name == "y_patch")
      continue;

    if (dim_name == "time") {
      stat = output.def_dim("time", PISM_UNLIMITED); check(stat);
    } else {
      stat = input.inq_dimlen(dim_name, dim_len); check(stat);
      stat = output.def_dim(dim_name, dim_len); check(stat);
    }
  }

  return 0;
}

//! \brief Defines variables in the output file.
int define_variables(PISMNC4_Serial &input, PISMNC4_Serial &output) {
  int stat, n_vars;

  stat = input.inq_nvars(n_vars); check(stat);

  for (int j = 0; j < n_vars; ++j) {
    vector<string> dimensions;
    string var_name;

    stat = input.inq_varname(j, var_name); check(stat);

    if (var_name == "x_patch" || var_name == "y_patch")
      continue;

    stat = input.inq_vardimid(var_name, dimensions); check(stat);

    for (unsigned int k = 0; k < dimensions.size(); ++k) {
      if (dimensions[k] == "x_patch")
        dimensions[k] = "x";

      if (dimensions[k] == "y_patch")
        dimensions[k] = "y";
    }

    PISM_IO_Type var_type;

    stat = input.inq_vartype(var_name, var_type); check(stat);

    stat = output.def_var(var_name, var_type, dimensions); check(stat);

    copy_attributes(input, output, var_name);
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

//! \brief Copies 1D variables (mostly coordinates).
/*!
 * Uses rank-0 data.
 */
int copy_coordinate_data(string filename, string var_name, PISMNC4_Serial &output) {
  int stat;
  unsigned int dim_len;
  vector<unsigned int> start, count;
  PISMNC4_Serial input(MPI_COMM_SELF, 0);

  stat = input.open(patch_filename(filename, 0), PISM_NOWRITE); check(stat);

  vector<string> dims;
  stat = input.inq_vardimid(var_name, dims); check(stat);

  assert(dims.size() == 1);
  stat = input.inq_dimlen(dims[0], dim_len); check(stat);

  start.push_back(0);
  count.push_back(dim_len);

  double *data = new double[dim_len];

  stat = input.get_vara_double(var_name, start, count, data); check(stat);
  stat = output.put_vara_double(var_name, start, count, data); check(stat);

  delete[] data;
  input.close();

  return 0;
}

//! \brief Copies 2D and 3D variables.
/*!
 * This is where most of the time is spent.
 *
 * Variables are copied one at a time (all patches of the first variable, then
 * all patches of the second, etc). This means that we end up opening and
 * closing input files many, many times. It might be better to process all the
 * variables in the first file, then all variables in the second file, etc, but
 * it is not clear it this access pattern of the output file is better or not.
 *
 * Note that patches may have different sizes. This code allocates and frees
 * the buffer for every read/write cycle. It might be better to compute the
 * maximum buffer size and allocate once, although it is also not clear if this
 * would give any performance benefit.
 */
int copy_spatial_data(string filename, string var_name,
                      vector<string> dims,
                      PISMNC4_Serial &output) {
  map<string, int> dim_lengths;
  PISMNC4_Serial input(MPI_COMM_SELF, 0);
  int stat;

  for (unsigned int d = 0; d < dims.size(); ++d) {
    unsigned int tmp;
    stat = output.inq_dimlen(dims[d], tmp); check(stat);
    dim_lengths[dims[d]] = tmp;
  }

  int mpi_size;
  stat = input.open(patch_filename(filename, 0), PISM_NOWRITE); check(stat);
  stat = get_quilt_size(input, mpi_size); check(stat);
  input.close();

  for (int r = 0; r < mpi_size; ++r) { // for each patch...
    int xs, ys;
    unsigned int xm, ym;

    stat = input.open(patch_filename(filename, r), PISM_NOWRITE); check(stat);

    stat = patch_geometry(input, xs, ys, xm, ym); check(stat);

    // for each time record...

    int max_time_start = dim_lengths["time"] > 0 ? dim_lengths["time"] : 1;

    for (int time_start = 0; time_start < max_time_start; ++time_start) {
      vector<unsigned int> in_start, out_start, count;

      // prepare start and count
      for (unsigned int d = 0; d < dims.size(); ++d) {
        // start
        if (dims[d] == "time") {
          in_start.push_back(time_start);
          out_start.push_back(time_start);
        } else if (dims[d] == "x") {
          in_start.push_back(0);
          out_start.push_back(xs);
        } else if (dims[d] == "y") {
          in_start.push_back(0);
          out_start.push_back(ys);
        } else {
          in_start.push_back(0);
          out_start.push_back(0);
        }

        // count
        if (dims[d] == "time") {
          count.push_back(1);
        } else if (dims[d] == "x") {
          count.push_back(xm);
        } else if (dims[d] == "y") {
          count.push_back(ym);
        } else {
          count.push_back(dim_lengths[dims[d]]);
        }
      }

      // allocate a buffer...
      unsigned long long int buffer_size = 1;
      for (unsigned int k = 0; k < count.size(); ++k)
        buffer_size *= count[k];

      double *data = (double*)malloc(sizeof(double) * buffer_size);
      if (data == NULL) {
        printf("ERROR: memory allocation failed while processing %s (variable %s)! Ending...\n",
               input.get_filename().c_str(), var_name.c_str());
        PISMEnd();
      }

      stat = input.get_vara_double(var_name, in_start, count, data); check(stat);

      stat = output.put_vara_double(var_name, out_start, count, data); check(stat);

      free(data);
    }

    input.close();
  }

  return 0;
}

//! \brief Copies variables.
int copy_data(string filename, PISMNC4_Serial &output) {
  int n_vars;
  int stat;
  PISMNC4_Serial input(MPI_COMM_SELF, 0);

  stat = output.inq_nvars(n_vars); check(stat);

  stat = input.open(patch_filename(filename, 0), PISM_NOWRITE); check(stat);
  input.close();

  for (int j = 0; j < n_vars; ++j) {
    vector<string> dimensions;
    string var_name;

    stat = output.inq_varname(j, var_name); check(stat);
    stat = output.inq_vardimid(var_name, dimensions); check(stat);

    // copy coordinate variables from the rank 0 file:
    if (dimensions.size() == 1) {
      stat = copy_coordinate_data(filename, var_name, output); check(stat);
    } else {
      // 2D or 3D variables
      stat = copy_spatial_data(filename, var_name, dimensions, output); check(stat);
    }
  }

  return 0;
}

int main(int argc, char *argv[])
{
  PetscErrorCode  ierr;

  MPI_Comm    com;
  PetscMPIInt rank, size;

  ierr = PetscInitialize(&argc, &argv, PETSC_NULL, help); CHKERRQ(ierr);

  com = PETSC_COMM_WORLD;
  ierr = MPI_Comm_rank(com, &rank); CHKERRQ(ierr);
  ierr = MPI_Comm_size(com, &size); CHKERRQ(ierr);

  /* This explicit scoping forces destructors to be called before PetscFinalize() */
  {
    ierr = verbosityLevelFromOptions(); CHKERRQ(ierr);

    ierr = verbPrintf(2,com, "PISM-MERGE %s (output file merging tool)\n",
		      PISM_Revision); CHKERRQ(ierr);
    ierr = stop_on_version_option(); CHKERRQ(ierr);

    bool i_set, o_set;
    string i_name, o_name;
    ierr = PISMOptionsString("-i", "Input file name",
                             i_name, i_set); CHKERRQ(ierr);
    ierr = PISMOptionsString("-o", "Output file name",
                             o_name, o_set); CHKERRQ(ierr);
    string usage =
      "  Merges output file created using '-o_format quilt'.\n\n"
      "  pismmerge {-i in.nc} [-o out.nc]\n"
      "where:\n"
      "  -i          in.nc is the name used with -extra_file or -o, e.g. ex-RANK.nc\n"
      "  -o          out.nc is the name of the output file that will contain merged data\n"
      "notes:\n"
      "  * -o is optional\n";

    vector<string> required;
    required.push_back("-i");
    ierr = show_usage_check_req_opts(com, "pismmerge", required, usage.c_str()); CHKERRQ(ierr);

    check_input_files(i_name);


    if (o_set == false) {
      o_name = output_filename(i_name);
    }

    if (rank == 0) {
      PISMNC4_Serial input(com, rank), output(com, rank);

      ierr = input.open(patch_filename(i_name, 0), PISM_NOWRITE); CHKERRQ(ierr);

      ierr = output.move_if_exists(o_name); CHKERRQ(ierr);
      ierr = output.create(o_name); CHKERRQ(ierr);

      ierr = define_dimensions(input, output); CHKERRQ(ierr);

      ierr = define_variables(input, output); CHKERRQ(ierr);

      copy_attributes(input, output, "PISM_GLOBAL");

      ierr = input.close(); CHKERRQ(ierr);

      ierr = copy_data(i_name, output); CHKERRQ(ierr);

      ierr = output.close(); CHKERRQ(ierr);
    }

    MPI_Barrier(com);
  }

  ierr = PetscFinalize(); CHKERRQ(ierr);

  return 0;
}
