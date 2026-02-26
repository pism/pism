// Copyright (C) 2020--2026 PISM Authors
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

#ifndef PISM_BLATTERTESTCASE_H
#define PISM_BLATTERTESTCASE_H

#include "pism/stressbalance/blatter/Blatter.hh"
#include "pism/geometry/Geometry.hh"
#include "pism/util/array/Array3D.hh"
#include "pism/util/array/Scalar.hh"
#include <memory>

namespace pism {

class Context;
class Grid;

namespace stressbalance {

/*! A BlatterTestCase manages running a Blatter solver instance against a particular
  verification test. Subclasses must implement:

  1) initializeGeometry (to initialize ice geometry, basal yield stress, and enthalpy)
  2) exactSolution (to return the exact velocity at a given 3D point)

  A driver uses a BlatterTestCase by calling 1-3 below and 4,5 as desired:

  1) its constructor
  2) init (to set coefficients)
  3) run (to actually solve)
  4) report (to print error norms)
  5) write (to save the results to a file)
*/
class BlatterTestCase {
public:
  BlatterTestCase(std::shared_ptr<Blatter> solver);

  virtual ~BlatterTestCase() = default;

  virtual void init();

  virtual void run();

  virtual void report(const std::string &testname);

  virtual void write(const std::string &filename);

protected:
  //! Set up geometry, yield stress, and enthalpy as appropriate for the test case.
  virtual void initializeGeometry() = 0;

  //! Return the exact solution at grid index (i,j) and sigma level k.
  //! The z coordinate is the physical elevation (bed + H * sigma).
  virtual Vector2d exactSolution(int i, int j, double x, double y,
                                 double z) const;

  std::shared_ptr<const pism::Grid> m_grid;

  const std::shared_ptr<const Context> m_ctx;
  const std::shared_ptr<const Config> m_config;

  const units::System::Ptr m_sys;

  Geometry m_geometry;
  array::Scalar m_tauc;
  array::Array3D m_ice_enthalpy;

  std::shared_ptr<Blatter> m_solver;
};

} // end of namespace stressbalance
} // end of namespace pism

#endif /* PISM_BLATTERTESTCASE_H */
