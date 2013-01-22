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

//! \brief Defines a dimension in an output file.
/*!
 * We assume that the time dimension is always unlimited (this is true about
 * PISM's output files.)
 *
 * Ignores dimensions that already exist in the output file or don't exist in
 * the input file.
 */
int define_dimension(PISMNC4_Serial &input, PISMNC4_Serial &output, string dim_name) {
  int stat;
  bool exists;

  stat = input.inq_dimid(dim_name, exists); check(stat);
  if (exists == false)
    return 0;

  stat = output.inq_dimid(dim_name, exists); check(stat);
  if (exists == true)
    return 0;

  if (dim_name == "time") {
    stat = output.def_dim("time", PISM_UNLIMITED); check(stat);
  } else {
    unsigned int dim_len;
    stat = input.inq_dimlen(dim_name, dim_len); check(stat);
    stat = output.def_dim(dim_name, dim_len); check(stat);
  }

  return 0;
}

//! \brief Defines a variable in an output file.
/*!
 * Also defines all the dimensions that this variable requires.
 *
 * The \c extra_vars output argument will contain names of coordinate variables
 * corresponding to dimensions used by this variable.
 */
int define_variable(PISMNC4_Serial &input, PISMNC4_Serial &output, string variable_name) {
  int stat;
  bool exists;
  vector<string> dimensions;

  stat = input.inq_varid(variable_name, exists); check(stat);
  if (exists == false)
    return 0;

  stat = output.inq_varid(variable_name, exists); check(stat);
  if (exists == true)
    return 0;

  stat = input.inq_vardimid(variable_name, dimensions); check(stat);

  for (unsigned int k = 0; k < dimensions.size(); ++k) {

    string dim_name = dimensions[k];

    if (dim_name == "x_patch")
      dim_name = "x";

    if (dim_name == "y_patch")
      dim_name = "y";

    stat = define_dimension(input, output, dim_name); check(stat);

    if (dim_name != variable_name) {
      stat = define_variable(input, output, dim_name); check(stat);
    }

    dimensions[k] = dim_name;
  }

  PISM_IO_Type var_type;
  stat = input.inq_vartype(variable_name, var_type); check(stat);

  stat = output.def_var(variable_name, var_type, dimensions); check(stat);

  stat = copy_attributes(input, output, variable_name); check(stat);

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
