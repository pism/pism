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

#include "pism/util/io/IO_Flags.hh"
#include "pism/util/Units.hh"

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

class YACInterpolation {
public:
  YACInterpolation(const Grid &target_grid, const File &input_file,
                   const std::string &variable_name);
  ~YACInterpolation();

  void regrid(const File &file, io::Default default_value, array::Scalar &target) const;

  double regrid(const pism::File &file, pism::io::Default default_value,
                const SpatialVariableMetadata &metadata, petsc::Vec &target) const;

  static std::string grid_name(const File &file, const std::string &variable_name,
                               units::System::Ptr sys);

private:
  double interpolate(const array::Scalar &source, petsc::Vec &target) const;

  static int interpolation_coarse_to_fine(double missing_value);
  static int interpolation_fine_to_coarse(double missing_value);

  static int define_field(int component_id, const Grid &pism_grid, const std::string &name);
  static int define_grid(const Grid &grid, const std::string &grid_name,
                         const std::string &projection);

  int m_instance_id;
  int m_source_field_id;
  int m_target_field_id;

  std::shared_ptr<array::Scalar> m_buffer;
};

} // namespace pism

#endif /* PISM_YACINTERPOLATION_H */
