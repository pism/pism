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

#include "bedrockOnlySystem.hh"


bedrockOnlySystemCtx::bedrockOnlySystemCtx(const NCConfigVariable &config, int my_Mbz)
      : columnSystemCtx(my_Mbz) {  // <- critical: sets size of sys
  Mbz = my_Mbz;

  // set values so we can check if init was called
  dtTemp  = -1.0;
  dzbEQ   = -1.0;
  bedR    = -1.0;
  Tbedtop = -1.0;
  Ghf     = -1.0;

  bed_rho = config.get("bedrock_thermal_density");                // kg m-3
  bed_c   = config.get("bedrock_thermal_specific_heat_capacity"); // J kg-1 K-1
  bed_k   = config.get("bedrock_thermal_conductivity");           // W m-1 K-1
  bedK = bed_k / (bed_rho * bed_c);
  
  Tb = new PetscScalar[Mbz];  // bedrock temps at prev step
}


bedrockOnlySystemCtx::~bedrockOnlySystemCtx() {
  delete [] Tb;
}


PetscErrorCode bedrockOnlySystemCtx::initAllColumns(
                  const PetscScalar my_dtTemp, const PetscScalar my_dzbEQ) {

  dtTemp = my_dtTemp;
  dzbEQ = my_dzbEQ;
  bedR = bedK * dtTemp / PetscSqr(dzbEQ);
  return 0;
}


PetscErrorCode bedrockOnlySystemCtx::setBoundaryValuesThisColumn(
                     const PetscScalar my_Tbedtop, const PetscScalar my_Ghf) {
  if ((dtTemp < 0.0) || (dzbEQ < 0.0) || (bedR < 0.0)) {
    SETERRQ(2,"setBoundaryValuesThisColumn() should only be called after\n"
              "  initAllColumns() in bedrockOnlySystemCtx"); }
  if ((Tbedtop >= 0.0) || (Ghf >= 0.0)) {  
    SETERRQ(3,"setBoundaryValuesThisColumn() called twice (?) in bedrockOnlySystemCtx"); }
  Tbedtop = my_Tbedtop;
  Ghf = my_Ghf;
  return 0;
}


PetscErrorCode bedrockOnlySystemCtx::viewConstants(
                     PetscViewer viewer, bool show_col_dependent) {
  PetscErrorCode ierr;

  if (!viewer) {
    ierr = PetscViewerASCIIGetStdout(PETSC_COMM_SELF,&viewer); CHKERRQ(ierr);
  }
  PetscTruth iascii;
  ierr = PetscTypeCompare((PetscObject)viewer,PETSC_VIEWER_ASCII,&iascii); CHKERRQ(ierr);
  if (!iascii) { SETERRQ(1,"Only ASCII viewer for bedrockOnlySystemCtx::viewConstants()\n"); }
  
  ierr = PetscViewerASCIIPrintf(viewer,
                  "\n<<VIEWING bedrockOnlySystemCtx:\n"); CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,
                     "for ALL columns:\n"
                     "  dtTemp,dzbEQ = %10.3e,%8.2f\n"
                     "  bed_rho,bed_c,bed_k = %10.3e,%10.3e,%10.3e\n"
                     "  bedK,bedR = %10.3e,%10.3e\n",
		     dtTemp,dzbEQ,bed_rho,bed_c,bed_k,bedK,bedR); CHKERRQ(ierr);
  if (show_col_dependent) {
    ierr = PetscViewerASCIIPrintf(viewer,
                     "for THIS column:\n"
                     "  Tbedtop,Ghf = %10.3f,%10.3e\n",
                     Tbedtop,Ghf); CHKERRQ(ierr);
  }
  ierr = PetscViewerASCIIPrintf(viewer,
                     ">>\n\n"); CHKERRQ(ierr);
  return 0;
}


/*! \brief Solve the tridiagonal system, in a single column, which determines
the bedrock temperature. */
PetscErrorCode bedrockOnlySystemCtx::solveThisColumn(PetscScalar **x) {

  if ((dtTemp < 0.0) || (dzbEQ < 0.0) || (bedR < 0.0)) {
    SETERRQ(2, "solveThisColumn() should only be called after\n"
               "  initAllColumns() in bedrockOnlySystemCtx"); }
  if ((Tbedtop < 0.0) || (Ghf < 0.0)) {  
    SETERRQ(3, "solveThisColumn() should only be called after\n"
               "  setBoundaryValuesThisColumn() in bedrockOnlySystemCtx"); }

  if (Mbz <= 1) {
    SETERRQ(4, "solveThisColumn() should only be called\n"
               "  if Mbz > 1 (for bedrockOnlySystemCtx)"); }

  // eqn:  - k_b (d Tb / d zb) = G + (heat equation); uses "add a point
  // past the end" trick (Morton & Mayers)
  // L[0] is not allocated
  D[0] = (1.0 + 2.0 * bedR);
  U[0] = - 2.0 * bedR;  
  rhs[0] = Tb[0] + 2.0 * dtTemp * Ghf / (bed_rho * bed_c * dzbEQ);

  // k=1:Mbz-2  bedrock only; pure vertical conduction problem
  for (PetscInt k = 1; k < Mbz-1; k++) {
    L[k] = -bedR;
    D[k] = 1.0 + 2.0 * bedR;
    U[k] = -bedR;
    rhs[k] = Tb[k];
  }
    
  // k=Mbz-1 equation says temperature at top of bedrock is known
  L[Mbz-1] = 0.0;
  D[Mbz-1] = 1.0;
  // U[Mbz-1] not allocated
  rhs[Mbz-1] = Tbedtop;

  PetscErrorCode retval = solveTridiagonalSystem(Mbz, x);  // solve it

  if (!retval) { // mark column as done by making b.c.s invalid
    Tbedtop = -1.0;
    Ghf     = -1.0;
  }
  return retval;
}


//! After solve: get flux at top of bedrock (e.g. to determine melt).
PetscScalar bedrockOnlySystemCtx::extractHeatFluxFromSoln(const PetscScalar *x) {
  return - bed_k * (x[Mbz-1] - x[Mbz-2]) / dzbEQ;
}

