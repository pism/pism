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

#ifndef __combinedSystem_hh
#define __combinedSystem_hh

#include "columnSystem.hh"

//! Tridiagonal linear system for vertical column of enthalpy-based conservation of energy in ice and bedrock.
/*!
This derived class is used only in case of cold ice base \e and a bedrock 
thermal layer of positive thickness (Mbz > 1).

See the page documenting \ref bombproofenth.
*/
class combinedSystemCtx : public columnSystemCtx {

public:
  combinedSystemCtx(
    const NCConfigVariable &config, IceModelVec3 &my_Enth3,
    int my_Mz, int my_Mbz);
  ~combinedSystemCtx();

  PetscErrorCode initAllColumns(
      const PetscScalar my_dx, const PetscScalar my_dy, 
      const PetscScalar my_dtTemp,
      const PetscScalar my_dzEQ, const PetscScalar my_dzbEQ);

  PetscErrorCode setSchemeParamsThisColumn(
      const bool my_ismarginal, const PetscScalar my_lambda);  
  PetscErrorCode setBoundaryValuesThisColumn(
      const PetscScalar my_Enth_surface, const PetscScalar my_Ghf);

  PetscErrorCode viewConstants(PetscViewer viewer, bool show_col_dependent);

  PetscErrorCode solveThisColumn(PetscScalar **x);

public:
  // arrays must be filled before calling solveThisColumn():
  PetscScalar  *Enth,   // enthalpy in ice at prev time step
               *Enth_s, // enthalpy level for CTS; function only of pressure
               *u,
               *v,
               *w,
               *Sigma,
               *Tb;     // temperature in bedrock at prev step

private:
  PetscInt     Mz, Mbz;
  PetscScalar  ice_rho, ice_c, ice_k,
               bed_rho, bed_c, bed_k,
               dx, dy, dtTemp, dzEQ, dzbEQ,
               nuEQ, iceK, iceRcold, bedK, bedR;
  IceModelVec3 *Enth3;
  PetscScalar  lambda, Enth_ks, Ghf;
  bool         ismarginal;
};

#endif   //  ifndef __combinedSystem_hh

