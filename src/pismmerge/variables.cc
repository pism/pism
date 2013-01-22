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

//! \brief Copies 1D variables.
/*!
 * This function processes coordinate variables and time bounds.
 */
int copy_coordinate_variable(PISMNC4_Serial &input, string var_name, PISMNC4_Serial &output) {
  int stat;
  unsigned int dim1_len = 0, dim2_len = 0;
  vector<unsigned int> start, count;
  vector<string> dims;
  double *data = NULL;

  stat = input.inq_vardimid(var_name, dims); check(stat);

  if (dims.size() == 1) {

    stat = input.inq_dimlen(dims[0], dim1_len); check(stat);

    start.push_back(0);
    count.push_back(dim1_len);

    assert(dim1_len > 0);
    data = new double[dim1_len];

    stat = input.get_vara_double(var_name, start, count, data); check(stat);
    stat = output.put_vara_double(var_name, start, count, data); check(stat);

  } else if (dims.size() == 2) {

    stat = input.inq_dimlen(dims[0], dim1_len); check(stat);
    stat = input.inq_dimlen(dims[1], dim2_len); check(stat);

    start.push_back(0);
    start.push_back(0);
    count.push_back(dim1_len);
    count.push_back(dim2_len);

    assert(dim1_len*dim2_len > 0);
    data = new double[dim1_len*dim2_len];

    stat = input.get_vara_double(var_name, start, count, data); check(stat);
    stat = output.put_vara_double(var_name, start, count, data); check(stat);

  }

  delete[] data;

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
 * the buffer each time it goes to the next patch. It might be better to compute the
 * maximum buffer size and allocate once, although it is also not clear if this
 * would give any performance benefit.
 */
int copy_spatial_variable(string filename, string var_name, PISMNC4_Serial &output) {
  map<string, int> dim_lengths;
  PISMNC4_Serial input(MPI_COMM_SELF, 0, 0);
  int stat;
  vector<string> dims;
  vector<unsigned int> in_start, out_start, count;

  stat = output.inq_vardimid(var_name, dims); check(stat);

  for (unsigned int d = 0; d < dims.size(); ++d) {
    unsigned int tmp;
    stat = output.inq_dimlen(dims[d], tmp); check(stat);
    dim_lengths[dims[d]] = tmp;
  }

  int mpi_size;
  stat = input.open(patch_filename(filename, 0), PISM_NOWRITE); check(stat);
  stat = get_quilt_size(input, mpi_size); check(stat);
  stat = input.close(); check(stat);

  for (int r = 0; r < mpi_size; ++r) { // for each patch...
    int xs, ys;
    unsigned int xm, ym;

    stat = input.open(patch_filename(filename, r), PISM_NOWRITE); check(stat);

    stat = patch_geometry(input, xs, ys, xm, ym); check(stat);

    int max_time_start = dim_lengths["time"] > 0 ? dim_lengths["time"] : 1;

    in_start.clear();
    count.clear();
    out_start.clear();
    // prepare start and count
    int time_idx = -1;
    for (unsigned int d = 0; d < dims.size(); ++d) {
      // start
      if (dims[d] == "time") {
        time_idx = d;
        in_start.push_back(0);
        out_start.push_back(0);
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

    // for each time record...
    for (int time_start = 0; time_start < max_time_start; ++time_start) {

      if (time_idx >= 0) {
        in_start[time_idx] = time_start;
        out_start[time_idx] = time_start;
      }

      stat = input.get_vara_double(var_name, in_start, count, data); check(stat);

      stat = output.put_vara_double(var_name, out_start, count, data); check(stat);
    }

    free(data);

    input.close();
  } // end of "for each patch..."

  return 0;
}

//! \brief Copies all variables.
/*!
 * Loops over variables present in an output file. This allows us to process
 * both cases ("-v foo" and without "-v").
 */
int copy_all_variables(string filename, PISMNC4_Serial &output) {
  int n_vars, stat;
  PISMNC4_Serial input(MPI_COMM_SELF, 0, 0);
  vector<string> dimensions, spatial_vars;

  stat = input.open(patch_filename(filename, 0), PISM_NOWRITE); check(stat);

  stat = output.inq_nvars(n_vars); check(stat);

  for (int j = 0; j < n_vars; ++j) {
    string var_name;

    stat = output.inq_varname(j, var_name); check(stat);
    stat = output.inq_vardimid(var_name, dimensions); check(stat);

    // copy coordinate variables from the rank 0 file:
    if (dimensions.size() == 1 || var_name == "time_bounds") {
      stat = copy_coordinate_variable(input, var_name, output); check(stat);
    } else {
      spatial_vars.push_back(var_name);
    }
  }

  stat = input.close(); check(stat);

  for (unsigned int k = 0; k < spatial_vars.size(); ++k) {
    // 2D or 3D variables
    fprintf(stderr, "Copying %s... ", spatial_vars[k].c_str());
    stat = copy_spatial_variable(filename, spatial_vars[k], output); check(stat);
    fprintf(stderr, "done.\n");
  }

  return 0;
}
