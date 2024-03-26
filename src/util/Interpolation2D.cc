/* Copyright (C) 2024 PISM Authors
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

#include "pism/util/Interpolation2D.hh"
#include "pism/util/io/IO_Flags.hh"

namespace pism {

/*!
 * Return the string that describes a 2D grid present in a NetCDF file.
 *
 * Here `variable_name` is the name of a 2D variable used to extract
 * grid information.
 *
 * We assume that a file may contain more than one grid, so the file
 * name alone is not sufficient.
 *
 * The output has the form "input_file.nc:y:x".
 */
std::string Interpolation2D::grid_name(const pism::File &file, const std::string &variable_name,
                                        pism::units::System::Ptr sys) {
  std::string result = file.name();
  for (const auto &d : file.dimensions(variable_name)) {
    auto type = file.dimension_type(d, sys);

    if (type == pism::X_AXIS or type == pism::Y_AXIS) {
      result += ":";
      result += d;
    }
  }
  return result;
}

double Interpolation2D::regrid(const pism::File &file, const SpatialVariableMetadata &metadata,
                               petsc::Vec &target) const {
  return regrid_impl(file, metadata, target);
}

Interpolation2D::Interpolation2D() {
  // empty
}

} // namespace pism
