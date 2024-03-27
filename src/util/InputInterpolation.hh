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

#ifndef PISM_INPUT_INTERPOLATION_H
#define PISM_INPUT_INTERPOLATION_H

#include <memory>
#include <vector>

#include "pism/util/Interpolation1D.hh" // InterpolationType
#include "pism/util/Units.hh"           // units::System::Ptr

namespace pism {

class File;
class SpatialVariableMetadata;
class Grid;
class LocalInterpCtx;

namespace petsc {
class Vec;
}

/*!
 * Interpolation from a 2D grid corresponding to a variable in an input file to a PISM's
 * internal 2D grid.
 */
class InputInterpolation {
public:
  virtual ~InputInterpolation() = default;

  double regrid(const pism::File &file, const SpatialVariableMetadata &metadata,
                petsc::Vec &output) const;

  static std::string grid_name(const File &file, const std::string &variable_name,
                               units::System::Ptr sys);

protected:
  InputInterpolation();
  virtual double regrid_impl(const pism::File &file, const SpatialVariableMetadata &metadata,
                             petsc::Vec &output) const = 0;
};

/*!
 * Legacy 2D and 3D interpolation code used to "regrid" (read with interpolation) inputs.
 *
 */
class InputInterpolation3D : public InputInterpolation {
public:
  InputInterpolation3D(std::shared_ptr<const Grid> target_grid, const std::vector<double> &levels,
                       const File &input_file, const std::string &variable_name,
                       InterpolationType type);

private:
  double regrid_impl(const pism::File &file, const SpatialVariableMetadata &metadata,
                     petsc::Vec &output) const;

  std::shared_ptr<const Grid> m_target_grid;
  std::shared_ptr<LocalInterpCtx> m_interp_context;
};

} // namespace pism

#endif /* PISM_INPUT_INTERPOLATION_H */
