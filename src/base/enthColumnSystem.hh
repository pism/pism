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

#ifndef __enthColumnSystem_hh
#define __enthColumnSystem_hh

#include "columnSystem.hh"

//! Tridiagonal linear system for vertical column of enthalpy-based conservation of energy.
class enthSystemCtx : public columnSystemCtx {

public:
  enthSystemCtx(int my_Mz, int my_Mbz);

  PetscErrorCode initAllColumns();
  PetscErrorCode viewConstants(PetscViewer viewer);

  PetscErrorCode setSchemeParamsThisColumn(
            const bool my_isfloating, const bool my_ismarginal,
            const PetscScalar my_lambda);  
  PetscErrorCode setSurfaceBoundaryValuesThisColumn(
            const PetscScalar my_Enth_surface);
  PetscErrorCode setBasalBoundaryValuesThisColumn(
            const PetscScalar my_Ghf, const PetscScalar my_Enth_shelfbase, 
            const PetscScalar my_Rb);

  PetscErrorCode solveThisColumn(PetscScalar **x);

public:
  // constants which should be set before calling initAllColumns()
  PetscScalar  dx,
               dy,
               dtTemp,
               dzEQ,
               dzbEQ,
               ice_rho,
               ice_c,
               ice_k,
               ice_nu,
               bed_rho,
               bed_c,
               bed_k;
  // pointers which should be set before calling initAllColumns()
  PetscScalar  *Enth,   // enthalpy in ice
               *Enth_s, // enthalpy level for CTS; function only of pressure
               *Tb,     // temperature in bedrock
               *u,
               *v,
               *w,
               *Sigma;
  IceModelVec3 *Enth3;

protected: // used internally
  PetscInt    Mz, Mbz;
  PetscScalar lambda, Enth_ks, Ghf, Enth_shelfbase, Rb;
  bool        isfloating, ismarginal;
  PetscScalar nuEQ,
              dzav,
              iceK,
              iceRcold,
              iceRtemp,
              bedK,
              bedR;

private:
  bool        initAllDone,
              schemeParamsValid,
              surfBCsValid,
              basalBCsValid;
};

#endif

