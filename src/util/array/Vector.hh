/* Copyright (C) 2020, 2021, 2022, 2023 PISM Authors
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
#ifndef PISM_ARRAY_VECTOR_HH
#define PISM_ARRAY_VECTOR_HH

#include "pism/util/array/Array2D.hh"
#include "pism/util/Vector2d.hh"

namespace pism {
namespace array {

class Scalar;

/** Class for storing and accessing 2D vector fields
*/
class Vector : public array::Array2D<pism::Vector2d> {
public:
  Vector(std::shared_ptr<const IceGrid> grid, const std::string &short_name);

  virtual ~Vector() = default;

  typedef std::shared_ptr<Vector> Ptr;
  typedef std::shared_ptr<const Vector> ConstPtr;

  std::shared_ptr<Vector> duplicate() const;
protected:
  Vector(std::shared_ptr<const IceGrid> grid, const std::string &name,
         unsigned int stencil_width);
};

class Vector1 : public Vector {
public:
  Vector1(std::shared_ptr<const IceGrid> grid, const std::string &name);

  typedef std::shared_ptr<Vector1> Ptr;
  typedef std::shared_ptr<const Vector1> ConstPtr;
  using Array2D<Vector2d>::star;
  using Array2D<Vector2d>::box;
protected:
  Vector1(std::shared_ptr<const IceGrid> grid, const std::string &name,
          unsigned int stencil_width);
};

class Vector2 : public Vector1 {
public:
  Vector2(std::shared_ptr<const IceGrid> grid, const std::string &name);

  typedef std::shared_ptr<Vector2> Ptr;
  typedef std::shared_ptr<const Vector2> ConstPtr;
};

void compute_magnitude(const array::Vector &input, array::Scalar &result);

} // end of namespace array
} // end of namespace pism

#endif /* PISM_ARRAY_VECTOR_HH */
