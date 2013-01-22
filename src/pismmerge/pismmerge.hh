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

#ifndef _PISMMERGE_H_
#define _PISMMERGE_H_

#include "PISMNC4_Serial.hh"
#include "pism_const.hh"
#include <vector>
#include <string>
#include <map>
#include <assert.h>
#include <stdlib.h>
#include <cstdio>

using namespace std;

// The following is a stupid kludge necessary to make NetCDF 4.x work in
// serial mode in an MPI program:
#ifndef MPI_INCLUDED
#define MPI_INCLUDED 1
#endif
#include <netcdf.h>


// metadata.cc
int define_dimension(PISMNC4_Serial &input, PISMNC4_Serial &output, string dim_name);
int define_variable(PISMNC4_Serial &input, PISMNC4_Serial &output, string variable_name);
int copy_attributes(PISMNC4_Serial &input, PISMNC4_Serial &output, string var_name);

// variables.cc
int copy_coordinate_variable(PISMNC4_Serial &input, string var_name, PISMNC4_Serial &output);
int copy_spatial_variable(string filename, string var_name, PISMNC4_Serial &output);
int copy_all_variables(string filename, PISMNC4_Serial &output);

// util.cc
void check(int return_code);
string patch_filename(string input, int mpi_rank);
string output_filename(string input, string var_name);
int get_quilt_size(PISMNC4_Serial &input, int &mpi_size);
int check_input_files(string filename);
int patch_geometry(PISMNC4_Serial &input, int &xs, int &ys, unsigned int &xm, unsigned int &ym);

#endif /* _PISMMERGE_H_ */
