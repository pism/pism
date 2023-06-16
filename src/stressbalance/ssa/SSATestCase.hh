// Copyright (C) 2009--2017, 2021, 2022, 2023 Ed Bueler, Constantine Khroulev and David Maxwell
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

#ifndef _SSATESTCASE_H_
#define _SSATESTCASE_H_

#include "SSA.hh"
#include "pism/util/Context.hh"
#include "pism/util/EnthalpyConverter.hh"
#include "pism/basalstrength/basal_resistance.hh"
#include "pism/util/Vars.hh"
#include "pism/util/IceGrid.hh"
#include "pism/util/array/Scalar.hh"
#include "pism/util/array/Vector.hh"
#include "pism/util/array/Scalar.hh"
#include "pism/util/array/Array3D.hh"
#include "pism/geometry/Geometry.hh"

namespace pism {
namespace stressbalance {

/*! An SSATestCase manages running an SSA instance against a particular
  test.  Subclasses must implement the following abstract methods to define
  the input to an SSA for a test case:

  1) initializeSSACoefficients (to initialize the ssa coefficients, e.g. ice thickness)

  Additionally, a subclass can implement `report` to handle
  printing statistics after a run.  The default report method relies
  on subclasses implementing the exactSolution method for comparison.

  A driver uses an SSATestCase by calling 1-3 below and 4,5 as desired:

  1) its constructor
  2) init (to set coefficients)
  3) run (to actually solve the ssa)
  4) report
  5) write (to save the results of the computation to a file)
*/
class SSATestCase
{
public:
  SSATestCase(std::shared_ptr<Context> ctx, int Mx, int My,
              double Lx, double Ly,
              grid::Registration registration,
              grid::Periodicity periodicity);

  virtual ~SSATestCase();

  virtual void init();

  virtual void run();

  virtual void report(const std::string &testname);

  virtual void write(const std::string &filename);

protected:

  //! Set up the coefficient variables as appropriate for the test case.
  virtual void initializeSSACoefficients() = 0;

  //! Return the value of the exact solution at grid index (i,j) or equivalently
  //! at coordinates (x,y).
  virtual void exactSolution(int i, int j,
                             double x, double y, double *u, double *v);

  void report_netcdf(const std::string &testname,
                     double max_vector,
                     double rel_vector,
                     double max_u,
                     double max_v,
                     double avg_u,
                     double avg_v);

  MPI_Comm m_com;
  const std::shared_ptr<Context> m_ctx;
  const Config::Ptr m_config;
  IceGrid::Ptr m_grid;
  const units::System::Ptr m_sys;

  //! "wide" stencil width
  int m_stencil_width;

  // SSA model variables.
  EnthalpyConverter::Ptr m_enthalpyconverter;

  // SSA coefficient variables.
  array::Scalar1 m_tauc;
  array::Array3D m_ice_enthalpy;

  array::Vector2 m_bc_values;
  array::Scalar2 m_bc_mask;

  Geometry m_geometry;

  SSA *m_ssa;
};

} // end of namespace stressbalance
} // end of namespace pism

#endif /* _SSATESTCASE_H_ */
