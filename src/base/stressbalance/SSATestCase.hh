// Copyright (C) 2009--2012 Ed Bueler, Constantine Khroulev and David Maxwell
//
// This file is part of PISM.
//
// PISM is free software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the Free Software
// Foundation; either version 2 of the License, or (at your option) any later
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
#include "enthalpyConverter.hh"
#include "basal_resistance.hh"
#include "PISMVars.hh"

//! Helper function for initializing a grid with the given dimensions and periodicity.
//! The grid is shallow (3 z-layers).
PetscErrorCode init_shallow_grid(IceGrid &grid, 
                                 PetscReal Lx, PetscReal Ly, 
                                 PetscInt Mx, PetscInt My, Periodicity p);



/*! An SSATestCase manages running an SSA instance against a particular
test.  Subclasses must implement the following abstract methods to define
the input to an SSA for a test case:

1) initializeGrid (to build a grid of the specified size appropriate for the test)
2) initializeSSAModel (to specify the laws used by the model, e.g. ice flow and basal sliding laws)
3) initializeSSACoefficients (to initialize the ssa coefficients, e.g. ice thickness)

The SSA itself is constructed between steps 2) and 3).

Additionally, a subclass can implement \c report to handle
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
  SSATestCase( MPI_Comm com, PetscMPIInt rank, 
               PetscMPIInt size, NCConfigVariable &c ): 
                  config(c), grid(com,rank,size,config), 
                  basal(0), enthalpyconverter(0), ssa(0),
                  report_velocity_scale(secpera)
  {  };

  virtual ~SSATestCase()
  {
    delete basal;
    delete enthalpyconverter;
    delete ssa;
  }

  virtual PetscErrorCode init(PetscInt Mx, PetscInt My,SSAFactory ssafactory);

  virtual PetscErrorCode run();

  virtual PetscErrorCode report(string testname);

  virtual PetscErrorCode write(const string &filename);

protected:

  virtual PetscErrorCode buildSSACoefficients();

  //! Initialize the member variable grid as appropriate for the test case.
  virtual PetscErrorCode initializeGrid(PetscInt Mx,PetscInt My) = 0;

  //! Allocate the member variables basal, ice, and enthalpyconverter as
  //! appropriate for the test case.
  virtual PetscErrorCode initializeSSAModel() = 0;

  //! Set up the coefficient variables as appropriate for the test case.
  virtual PetscErrorCode initializeSSACoefficients() = 0;

  //! Return the value of the exact solution at grid index (i,j) or equivalently
  //! at coordinates (x,y).
  virtual PetscErrorCode exactSolution(PetscInt i, PetscInt j,
    PetscReal x, PetscReal y, PetscReal *u, PetscReal *v );

  PetscErrorCode report_netcdf(string testname,
                               double max_vector,
                               double rel_vector,
                               double max_u,
                               double max_v,
                               double avg_u,
                               double avg_v);
  NCConfigVariable &config;
  IceGrid grid;

  // SSA model variables.
  IceBasalResistancePlasticLaw *basal;
  EnthalpyConverter *enthalpyconverter;

  // SSA coefficient variables.
  PISMVars vars;
  IceModelVec2S  surface, thickness, bed, tauc;
  IceModelVec3 enthalpy;
  IceModelVec2V vel_bc;
  IceModelVec2Int ice_mask, bc_mask;

  SSA *ssa;

  // Scale for converting velocities from their units during computation to human friendly units.
  PetscScalar report_velocity_scale;

};

#endif /* _SSATESTCASE_H_ */
