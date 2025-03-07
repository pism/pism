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

namespace pism {

class File;
class SpatialVariableMetadata;
class Grid;
class LocalInterpCtx;

namespace petsc {
class Vec;
}

/*!
 * Interpolation from a grid corresponding to a variable in an input file to a PISM's
 * internal grid.
 */
class InputInterpolation {
public:
  virtual ~InputInterpolation() = default;

  static std::shared_ptr<InputInterpolation>
  create(const Grid &target_grid, const std::vector<double> &levels, const File &input_file,
         const std::string &variable_name, InterpolationType type);

  /*!
   * Read a record `record_index` of the variable corresponding to the provided `metadata`
   * from a `file` and interpolate onto the target grid. Store results in `output`.
   *
   * Set `record_index` to -1 to read the last record available in `file`.
   *
   */
  double regrid(const SpatialVariableMetadata &metadata, const pism::File &file, int record_index,
                const Grid &grid, petsc::Vec &output) const;

protected:
  InputInterpolation();
  virtual double regrid_impl(const SpatialVariableMetadata &metadata, const pism::File &file,
                             int record_index, const Grid &grid, petsc::Vec &output) const = 0;
};

/*!
 * Legacy 2D and 3D interpolation code used to "regrid" (read with interpolation) inputs.
 *
 */
class InputInterpolation3D : public InputInterpolation {
public:
  InputInterpolation3D(const Grid &target_grid, const std::vector<double> &levels,
                       const File &input_file, const std::string &variable_name,
                       InterpolationType type);

private:
  double regrid_impl(const SpatialVariableMetadata &metadata, const pism::File &file,
                     int record_index, const Grid &grid, petsc::Vec &output) const;

  std::shared_ptr<LocalInterpCtx> m_interp_context;
};

} // namespace pism

#endif /* PISM_INPUT_INTERPOLATION_H */
