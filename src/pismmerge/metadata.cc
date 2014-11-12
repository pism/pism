// Copyright (C) 2013, 2014 PISM Authors
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

//! \brief Defines a dimension in an output file.
/*!
 * We assume that the time dimension is always unlimited (this is true about
 * PISM's output files.)
 *
 * Ignores dimensions that already exist in the output file or don't exist in
 * the input file.
 */
int define_dimension(pism::NC4_Serial &input, pism::NC4_Serial &output, std::string dim_name) {
  bool exists;

  input.inq_dimid(dim_name, exists);
  if (exists == false) {
    return 0;
  }

  output.inq_dimid(dim_name, exists);
  if (exists == true) {
    return 0;
  }

  if (dim_name == "time") {
    output.def_dim("time", pism::PISM_UNLIMITED);
  } else {
    unsigned int dim_len;
    input.inq_dimlen(dim_name, dim_len);
    output.def_dim(dim_name, dim_len);
  }

  return 0;
}

//! \brief Defines a variable in an output file.
/*!
 * Also defines all the dimensions that this variable requires.
 *
 * The `extra_vars` output argument will contain names of coordinate variables
 * corresponding to dimensions used by this variable.
 */
int define_variable(pism::NC4_Serial &input, pism::NC4_Serial &output, std::string variable_name) {
  int stat;
  bool exists;
  std::vector<std::string> dimensions;

  input.inq_varid(variable_name, exists);
  if (exists == false) {
    return 0;
  }

  output.inq_varid(variable_name, exists);
  if (exists == true) {
    return 0;
  }

  input.inq_vardimid(variable_name, dimensions);

  for (unsigned int k = 0; k < dimensions.size(); ++k) {

    std::string dim_name = dimensions[k];

    if (dim_name == "x_patch") {
      dim_name = "x";
    }

    if (dim_name == "y_patch") {
      dim_name = "y";
    }

    stat = define_dimension(input, output, dim_name); check(stat);

    if (dim_name != variable_name) {
      stat = define_variable(input, output, dim_name); check(stat);
    }

    dimensions[k] = dim_name;
  }

  pism::IO_Type var_type;
  input.inq_vartype(variable_name, var_type);

  output.def_var(variable_name, var_type, dimensions);

  stat = copy_attributes(input, output, variable_name); check(stat);

  return 0;
}

//! \brief Copies variable attributes.
int copy_attributes(pism::NC4_Serial &input, pism::NC4_Serial &output, std::string var_name) {
  int n_attrs;

  input.inq_varnatts(var_name, n_attrs);

  for (int j = 0; j < n_attrs; ++j) {
    std::string att_name;
    pism::IO_Type att_type;

    input.inq_attname(var_name, j, att_name);

    input.inq_atttype(var_name, att_name, att_type);

    if (att_type == pism::PISM_CHAR) {
      std::string tmp;

      input.get_att_text(var_name, att_name, tmp);

      output.put_att_text(var_name, att_name, tmp);
    } else {
      std::vector<double> tmp;

      input.get_att_double(var_name, att_name, tmp);

      output.put_att_double(var_name, att_name, att_type, tmp);
    }
  }

  return 0;
}
