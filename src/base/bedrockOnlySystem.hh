// Copyright (C) 2010 Ed Bueler
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

#ifndef __bedrockOnlySystem_hh
#define __bedrockOnlySystem_hh

#include "columnSystem.hh"
#include "NCVariable.hh"

//! Tridiagonal linear system for conservation of energy in vertical column of bedrock only.
/*!
See the page documenting \ref bombproofenth.  We implement equations
\ref bedrockeqn and \ref geothermalbedeqn.  The top of the bedrock has a
Dirichlet condition and the bottom has a Neumann condition.  Heat flux at top of
bedrock can be extracted from the solution.
*/
class bedrockOnlySystemCtx : public columnSystemCtx {

public:
  bedrockOnlySystemCtx(const NCConfigVariable &config, int my_Mbz);
  ~bedrockOnlySystemCtx();

  PetscErrorCode initAllColumns(
      const PetscScalar my_dtTemp, const PetscScalar my_dzbEQ);

  PetscErrorCode setBoundaryValuesThisColumn(
      const PetscScalar my_Tbedtop, const PetscScalar my_Ghf);

  PetscErrorCode viewConstants(PetscViewer viewer, bool show_col_dependent);

  PetscErrorCode solveThisColumn(PetscScalar **x);
  PetscScalar    extractHeatFluxFromSoln(const PetscScalar *x);

public:
  // array must be filled before calling solveThisColumn():
  PetscScalar  *Tb;     // temperature in bedrock at prev step

private:
  PetscInt    Mbz;
  PetscScalar dtTemp, dzbEQ, bed_rho, bed_c, bed_k, bedK, bedR, Ghf, Tbedtop;
};

#endif  // ifndef __bedrockOnlySystem_hh


