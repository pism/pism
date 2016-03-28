// Copyright (C) 2009--2016 Ed Bueler, Constantine Khroulev and David Maxwell
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
#include "base/util/Context.hh"
#include "base/enthalpyConverter.hh"
#include "base/basalstrength/basal_resistance.hh"
#include "base/util/PISMVars.hh"
#include "base/util/IceGrid.hh"
#include "base/util/IceModelVec2CellType.hh"

namespace pism {
namespace stressbalance {

/*! An SSATestCase manages running an SSA instance against a particular
  test.  Subclasses must implement the following abstract methods to define
  the input to an SSA for a test case:

  1) initializeGrid (to build a grid of the specified size appropriate for the test)
  2) initializeSSAModel (to specify the laws used by the model, e.g. ice flow and basal sliding laws)
  3) initializeSSACoefficients (to initialize the ssa coefficients, e.g. ice thickness)

  The SSA itself is constructed between steps 2) and 3).

  Additionally, a subclass can implement `report` to handle
  printing statistics after a run.  The default report method relies
  on subclasses implementing the exactSolution method for comparision.

  A driver uses an SSATestCase by calling 1-3 below and 4,5 as desired:

  1) its constructor
  2) init (to specify the grid size and choice of SSA algorithm)
  3) run (to actually solve the ssa)
  4) report
  5) write (to save the results of the computation to a file)
*/
class SSATestCase
{
public:
  SSATestCase(Context::Ptr ctx);

  virtual ~SSATestCase();

  virtual void init(int Mx, int My, SSAFactory ssafactory);

  virtual void run();

  virtual void report(const std::string &testname);

  virtual void write(const std::string &filename);

protected:

  virtual void buildSSACoefficients();

  //! Initialize the member variable grid as appropriate for the test case.
  virtual void initializeGrid(int Mx,int My) = 0;

  //! Allocate the member variables basal, ice, and enthalpyconverter as
  //! appropriate for the test case.
  virtual void initializeSSAModel() = 0;

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
  const Config::Ptr m_config;
  IceGrid::Ptr m_grid;
  const Context::Ptr m_ctx;
  const units::System::Ptr m_sys;

  // SSA model variables.
  EnthalpyConverter::Ptr m_enthalpyconverter;

  // SSA coefficient variables.
  IceModelVec2S  m_surface, m_thickness, m_bed, m_tauc, m_melange_back_pressure;
  IceModelVec3 m_ice_enthalpy;
  IceModelVec2V m_bc_values;
  IceModelVec2CellType m_ice_mask;
  IceModelVec2Int m_bc_mask;

  SSA *m_ssa;
};

} // end of namespace stressbalance
} // end of namespace pism

#endif /* _SSATESTCASE_H_ */
