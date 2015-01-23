// Copyright (C) 2013, 2014, 2015 PISM Authors
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

#ifndef _PISMMERGE_H_
#define _PISMMERGE_H_

#include "PISMNC4_Serial.hh"
#include "pism_const.hh"
#include <vector>
#include <string>
#include <map>
#include <stdlib.h>
#include <cstdio>

// The following is a stupid kludge necessary to make NetCDF 4.x work in
// serial mode in an MPI program:
#ifndef MPI_INCLUDED
#define MPI_INCLUDED 1
#endif
#include <netcdf.h>


// metadata.cc
void define_dimension(const pism::NC4_Serial &input, const pism::NC4_Serial &output,
                      const std::string &dim_name);
void define_variable(const pism::NC4_Serial &input, const pism::NC4_Serial &output,
                     const std::string &variable_name);
void copy_attributes(const pism::NC4_Serial &input, const pism::NC4_Serial &output,
                     const std::string &var_name);

// variables.cc
void copy_coordinate_variable(const pism::NC4_Serial &input, const std::string &var_name,
                              const pism::NC4_Serial &output);
void copy_spatial_variable(const std::string &filename, const std::string &var_name,
                           const pism::NC4_Serial &output);
void copy_all_variables(const std::string &filename, const pism::NC4_Serial &output);

// util.cc
std::string patch_filename(const std::string &input, int mpi_rank);
std::string output_filename(const std::string &input, const std::string &var_name);
int get_quilt_size(const pism::NC4_Serial &input);
void check_input_files(const std::string &filename);
void patch_geometry(const pism::NC4_Serial &input, int &xs, int &ys,
                    unsigned int &xm, unsigned int &ym);

#endif /* _PISMMERGE_H_ */
