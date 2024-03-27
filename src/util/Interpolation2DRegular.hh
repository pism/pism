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

#include "pism/util/Interpolation1D.hh"
#include "pism/util/io/LocalInterpCtx.hh"
#include "pism/util/Interpolation2D.hh"
#include <memory>
#include <vector>

namespace pism {

class Interpolation2DRegular : public Interpolation2D {
public:
  Interpolation2DRegular(std::shared_ptr<const Grid> target_grid,
                         const std::vector<double> &levels,
                         const File &input_file,
                         const std::string &variable_name, InterpolationType type);

private:
  double regrid_impl(const pism::File &file, const SpatialVariableMetadata &metadata,
                     petsc::Vec &output) const;

  std::shared_ptr<const Grid> m_target_grid;
  std::shared_ptr<LocalInterpCtx> m_interp_context;
};

} // namespace pism
