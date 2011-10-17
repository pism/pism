// Copyright (C) 2009-2011 Andreas Aschwanden and Ed Bueler
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

#ifndef __enthSystem_hh
#define __enthSystem_hh

#include <vector>

#include "columnSystem.hh"

class NCConfigVariable;
class IceModelVec3;

//! Tridiagonal linear system for conservation of energy in vertical column of ice enthalpy.
/*!
See the page documenting \ref bombproofenth.  The top of
the ice has a Dirichlet condition.  

FIXME:  IS THIS THE RIGHT STRUCTURE?:  The boundary condition at the bottom of the
ice depends on various cases, and these cases are not decided here.  Instead, 
the user of this class sets the lowest-level (z=0) equation.
*/
class enthSystemCtx : public columnSystemCtx {

public:
  enthSystemCtx(const NCConfigVariable &config, IceModelVec3 &my_Enth3, int my_Mz, string my_prefix);
  virtual ~enthSystemCtx();

  PetscErrorCode initAllColumns(PetscScalar my_dx, PetscScalar my_dy, 
                                PetscScalar my_dtTemp, PetscScalar my_dzEQ);

  PetscErrorCode initThisColumn(bool my_ismarginal,
                                PetscScalar my_lambda,
                                PetscReal ice_thickness);  

  PetscScalar k_from_T(PetscScalar T);

  PetscErrorCode setBoundaryValuesThisColumn(PetscScalar my_Enth_surface);
  PetscErrorCode setDirichletBasal(PetscScalar Y);
  PetscErrorCode setBasalHeatFlux(PetscScalar hf);

  PetscErrorCode viewConstants(PetscViewer viewer, bool show_col_dependent);
  PetscErrorCode viewSystem(PetscViewer viewer) const;

  PetscErrorCode solveThisColumn(PetscScalar **x, PetscErrorCode &pivoterrorindex);

public:
  // arrays must be filled before calling solveThisColumn():
  PetscScalar  *Enth,   // enthalpy in ice at prev time step
               *Enth_s, // enthalpy level for CTS; function only of pressure
               *u,
               *v,
               *w,
               *Sigma;
protected:
  PetscInt     Mz;
  PetscScalar  ice_rho, ice_c, ice_k, ice_K, ice_K0,
               dx, dy, dtTemp, dzEQ, nuEQ, iceK, iceRcold, iceRtemp;
  IceModelVec3 *Enth3;
  PetscScalar  lambda, Enth_ks, a0, a1, b;
  bool         ismarginal;
  vector<PetscScalar> R; // value of k \Delta t / (\rho c \Delta x^2) in
                         // column, using current values of enthalpy

  virtual PetscErrorCode assemble_R();
  PetscErrorCode checkReadyToSolve();
};

#endif   //  ifndef __enthSystem_hh

