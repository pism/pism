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

#include "pism/util/VariableMetadata.hh"
#include "pism/util/petscwrappers/Vec.hh"
#include "pism/util/io/File.hh"

namespace pism {

/*!
 * Interpolation from a 2D grid corresponding to a variable in an input file to a PISM's
 * internal 2D grid.
 */
class Interpolation2D {
public:
  virtual ~Interpolation2D() = default;

  double regrid(const pism::File &file,
                const SpatialVariableMetadata &metadata, petsc::Vec &target) const;

  static std::string grid_name(const File &file, const std::string &variable_name,
                               units::System::Ptr sys);
protected:
  Interpolation2D();
  virtual double regrid_impl(const pism::File &file, const SpatialVariableMetadata &metadata,
                             petsc::Vec &target) const = 0;
};

}
