// Copyright (C) 2009-2010 Andreas Aschwanden and Ed Bueler
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

#ifndef __iceenthOnlySystem_hh
#define __iceenthOnlySystem_hh

#include "columnSystem.hh"

//! Tridiagonal linear system for conservation of energy in vertical column of ice enthalpy.
/*!
See the page documenting \ref bombproofenth.  This class is used either
when the ice base is at the pressure melting temperature or when there is no
bedrock thermal layer of positive thickness (i.e. because Mbz==1).  The top of
the ice has a Dirichlet condition.  The boundary condition at the bottom of the
ice depends on various cases, and these cases are not decided here.  Instead, 
the user of this class sets the lowest-level (z=0) equation.
*/
class iceenthOnlySystemCtx : public columnSystemCtx {

public:
  iceenthOnlySystemCtx(
    const NCConfigVariable &config, IceModelVec3 &my_Enth3, int my_Mz);
  ~iceenthOnlySystemCtx();

  PetscErrorCode initAllColumns(
      const PetscScalar my_dx, const PetscScalar my_dy, 
      const PetscScalar my_dtTemp, const PetscScalar my_dzEQ);

  PetscErrorCode setSchemeParamsThisColumn(
      const bool my_ismarginal, const PetscScalar my_lambda);  
  PetscErrorCode setBoundaryValuesThisColumn(
      const PetscScalar my_Enth_surface);
  PetscErrorCode setLevel0EqnThisColumn(
      const PetscScalar my_a0, const PetscScalar my_a1, const PetscScalar my_b);

  PetscErrorCode viewConstants(PetscViewer viewer, bool show_col_dependent);

  PetscErrorCode solveThisColumn(PetscScalar **x);

public:
  // arrays must be filled before calling solveThisColumn():
  PetscScalar  *Enth,   // enthalpy in ice at prev time step
               *Enth_s, // enthalpy level for CTS; function only of pressure
               *u,
               *v,
               *w,
               *Sigma;

private:
  PetscInt     Mz;
  PetscScalar  ice_rho, ice_c, ice_k, ice_nu,
               dx, dy, dtTemp, dzEQ, nuEQ, iceK, iceRcold, iceRtemp;
  IceModelVec3 *Enth3;
  PetscScalar  lambda, Enth_ks, a0, a1, b;
  bool         ismarginal;
};

#endif   //  ifndef __iceenthOnlySystem_hh

