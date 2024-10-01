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

#ifndef PISM_YACINTERPOLATION_H
#define PISM_YACINTERPOLATION_H

#include <memory>
#include <string>

#include "pism/util/InputInterpolation.hh"

namespace pism {
class Grid;
class File;
class SpatialVariableMetadata;

namespace array {
class Scalar;
}

namespace petsc {
class Vec;
}

/*!
 * Interpolation from a Cartesian projected grid in an `input_file`.
 *
 * An `input_file` has to contain projection information (any of the options supported by
 * MappingInfo::FromFile()).
 */
class InputInterpolationYAC : public InputInterpolation {
public:
  InputInterpolationYAC(const Grid &target_grid, const File &input_file,
                        const std::string &variable_name);
  virtual ~InputInterpolationYAC();

  void regrid(const File &file, array::Scalar &output) const;

private:
  double regrid_impl(const SpatialVariableMetadata &metadata, const pism::File &file,
                     int record_index, const Grid &target_grid, petsc::Vec &output) const;

  double interpolate(const array::Scalar &source, petsc::Vec &target) const;

  static int interpolation_coarse_to_fine(double missing_value);
  static int interpolation_fine_to_coarse(double missing_value);

  static int define_field(int component_id, const std::vector<double> &x,
                          const std::vector<double> &y, const std::string &proj_string,
                          const std::string &name);
  static int define_grid(const std::vector<double> &x, const std::vector<double> &y,
                         const std::string &grid_name, const std::string &projection);

  int m_instance_id;
  int m_source_field_id;
  int m_target_field_id;

  std::shared_ptr<array::Scalar> m_buffer;
};

} // namespace pism

#endif /* PISM_YACINTERPOLATION_H */
