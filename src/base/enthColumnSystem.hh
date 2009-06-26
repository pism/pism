// Copyright (C) 2009 Andreas Aschwanden and Ed Bueler
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
  PetscErrorCode setIndicesThisColumn(PetscInt i, PetscInt j, PetscInt ks);  
  PetscErrorCode setSchemeParamsThisColumn(
                     PetscScalar my_mask, bool my_isMarginal, PetscScalar my_lambda);  
  PetscErrorCode setSurfaceBoundaryValuesThisColumn(PetscScalar my_Enth_surface);
  PetscErrorCode setBasalBoundaryValuesThisColumn(
                     PetscScalar my_Ghf, PetscScalar my_Enth_shelfbase, PetscScalar my_Rb);
  PetscErrorCode solveThisColumn(PetscScalar **x);
  
  PetscErrorCode view(MPI_Comm &com);

public:
  // constants which should be set before calling initForAllColumns()
  PetscScalar  dx,
               dy,
               dtTemp,
               dzEQ,
               dzbEQ,
               ice_rho,
               ice_c,
               ice_k,
               bed_thermal_rho,
               bed_thermal_c,
               bed_thermal_k;
  // pointers which should be set before calling initForAllColumns()
  PetscScalar  *Enth,   // enthalpy in ice
               *Enth_s, // enthalpy level for CTS; function only of pressure
               *Enth_b, // enthalpy in bedrock
               *u,
               *v,
               *w,
               *Sigma;
  IceModelVec3 *Enth3;

protected: // used internally
  PetscInt    Mz, Mbz, k0;
  PetscInt    i, j, ks;
  PetscScalar mask, lambda, Enth_ks, Ghf, Enth_shelfbase, Rb;
  bool        isMarginal;
  PetscScalar nuEQ,
              dzav,
              rho_c_I,
              rho_c_br,
              iceK,
              iceR,
              brK,
              brR;
  bool        initAllDone,
              indicesValid,
              schemeParamsValid,
              surfBCsValid,
              basalBCsValid;
};

#endif


