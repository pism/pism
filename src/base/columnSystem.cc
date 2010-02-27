// Copyright (C) 2004-2010 Jed Brown, Ed Bueler and Constantine Khroulev
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

#include <petsc.h>
#include "pism_const.hh"   // e.g. MASK_FLOATING
#include "iceModelVec.hh"
#include "columnSystem.hh"


//! Allocate a tridiagonal system of maximum size my_nmax.
/*!
Let N = \c nmax.  Then allocated locations are like this:
\verbatim
D[0]   U[0]    0      0      0    ...
L[1]   D[1]   U[1]    0      0    ...
 0     L[2]   D[2]   U[2]    0    ...
 0      0     L[3]   D[3]   U[3]  ...
\endverbatim
with the last row
\verbatim
0       0     ...     0  L[N-1]  D[N-1]
\endverbatim
Thus the index into the arrays L, D, U is always the row number.

Note L[0] is not allocated and U[N-1] is not allocated.
 */
columnSystemCtx::columnSystemCtx(PetscInt my_nmax) : nmax(my_nmax) {
  if (nmax < 1) {
    PetscPrintf(PETSC_COMM_WORLD,
      "columnSystem ERROR: nmax of system too small\n");
    PetscEnd();
  }
  if (nmax > 1000000) {
    PetscPrintf(PETSC_COMM_WORLD,
      "columnSystem ERROR: nmax of system unreasonable (> 10^6)\n");
    PetscEnd();
  }
  Lp   = new PetscScalar[nmax-1];
  L    = Lp-1; // ptr arith.; note L[0]=Lp[-1] not allocated
  D    = new PetscScalar[nmax];
  U    = new PetscScalar[nmax-1];
  rhs  = new PetscScalar[nmax];
  work = new PetscScalar[nmax];

  resetColumn();

  indicesValid = false;
}


columnSystemCtx::~columnSystemCtx() {
  delete [] Lp;
  delete [] D;
  delete [] U;
  delete [] rhs;
  delete [] work;
}


//! Zero all entries.
PetscErrorCode columnSystemCtx::resetColumn() {
  // done this way for memory contiguity:
  for (PetscInt k = 0; k < nmax-1; k++) Lp[k]   = 0.0;
  for (PetscInt k = 0; k < nmax-1; k++) U[k]    = 0.0;
  for (PetscInt k = 0; k < nmax;   k++) D[k]    = 0.0;
  for (PetscInt k = 0; k < nmax;   k++) rhs[k]  = 0.0;
  for (PetscInt k = 0; k < nmax;   k++) work[k] = 0.0;
  return 0;
}


//! Compute 1-norm, which is max sum of absolute values of columns.
PetscScalar columnSystemCtx::norm1(const PetscInt n) const {
  if (n > nmax) {
    PetscPrintf(PETSC_COMM_WORLD,"PISM ERROR:  n > nmax in columnSystemCtx::norm1()\n");
    PetscEnd();
  }
  if (n == 1)  return fabs(D[0]);   // only 1x1 case is special
  PetscScalar z = fabs(D[0]) + fabs(L[1]);
  for (PetscInt k = 1; k < n; k++) {  // k is column index (zero-based)
    z = PetscMax(z, fabs(U[k-1])) + fabs(D[k]) + fabs(L[k+1]);
  }
  z = PetscMax(z, fabs(U[n-2]) + fabs(D[n-1]));
  return z;
}


//! Compute diagonal-dominance ratio.  If this is less than one then the matrix is strictly diagonally-dominant.
/*!
Let \f$A = (a_{ij})\f$ be the tridiagonal matrix
described by L,D,U for row indices 0 through \c n.  The computed ratio is 
  \f[ \max_{j=1,\dots,n} \frac{|a_{j,j-1}|+|a_{j,j+1}|}{|a_{jj}|}, \f]
where \f$a_{1,0}\f$ and \f$a_{n,n+1}\f$ are interpreted as zero.

If this is smaller than one then it is a theorem that the tridiagonal solve will
succeed.

We return -1.0 if the absolute value of any diagonal element is less than
1e-12 of the 1-norm of the matrix.
 */
PetscScalar columnSystemCtx::ddratio(const PetscInt n) const {
  if (n > nmax) {
    PetscPrintf(PETSC_COMM_WORLD,"PISM ERROR:  n > nmax in columnSystemCtx::ddratio()\n");
    PetscEnd();
  }
  const PetscScalar scale = norm1(n);

  if ( (fabs(D[0]) / scale) < 1.0e-12)  return -1.0;
  PetscScalar z = fabs(U[0]) / fabs(D[0]);

  for (PetscInt j = 1; j < n-1; j++) {  // j is row index (zero-based)
    if ( (fabs(D[j]) / scale) < 1.0e-12)  return -1.0;
    const PetscScalar s = fabs(L[j]) + fabs(U[j]);
    z = PetscMax(z, s / fabs(D[j]) );
  } 

  if ( (fabs(D[n-1]) / scale) < 1.0e-12)  return -1.0;
  z = PetscMax(z, fabs(L[n-1]) / fabs(D[n-1]) );

  return z;
}


PetscErrorCode columnSystemCtx::setIndicesAndClearThisColumn(
                  PetscInt my_i, PetscInt my_j, PetscInt my_ks) {
  if (indicesValid) {  SETERRQ(3,
     "setIndicesAndClearThisColumn() called twice in same column"); }
  i = my_i;
  j = my_j;
  ks = my_ks;

  resetColumn();
  
  indicesValid = true;
  return 0;
}


//! Utility for simple ascii view of column quantity.
/*!
Give first argument NULL to get standard out.  No binary viewer.

Give description string as \c info argument.
 */
PetscErrorCode columnSystemCtx::viewColumnValues(PetscViewer viewer,
                                                 PetscScalar *v, PetscInt m, const char* info) const {
  PetscErrorCode ierr;

  if (v==NULL) {
    SETERRQ1(2,"columnSystem ERROR: can't view '%s' by v=NULL pointer ... ending ...\n", info);
  }
  if ((m<1) || (m>1000)) {
    SETERRQ1(3,"columnSystem ERROR: can't view '%s' because m<1 or m>1000 column values ... ending ...\n",info);
  }

  PetscTruth iascii;
  if (!viewer) {
    ierr = PetscViewerASCIIGetStdout(PETSC_COMM_SELF,&viewer); CHKERRQ(ierr);
  }
  ierr = PetscTypeCompare((PetscObject)viewer,PETSC_VIEWER_ASCII,&iascii); CHKERRQ(ierr);
  if (!iascii) { SETERRQ(1,"Only ASCII viewer for ColumnSystem\n"); }

  ierr = PetscViewerASCIIPrintf(viewer,
     "\n<viewing ColumnSystem column object with description '%s':\n"
      "  k     value\n",info); CHKERRQ(ierr);
  for (PetscInt k=0; k<m; k++) {
    ierr = PetscViewerASCIIPrintf(viewer,
      "  %5d %12.5e\n",k,v[k]); CHKERRQ(ierr);
  }
  ierr = PetscViewerASCIIPrintf(viewer,
      "end viewing>\n"); CHKERRQ(ierr);
  return 0;
}


//! View the tridiagonal matrix as a full matrix.  Unwise if size exceeds 20 or so. 
/*! 
Give first argument NULL to get standard out.  No binary viewer.

Give description string as \c info argument.
 */
PetscErrorCode columnSystemCtx::viewMatrix(PetscViewer viewer, const char* info) const {
  PetscErrorCode ierr;
  PetscTruth iascii;
  if (!viewer) {
    ierr = PetscViewerASCIIGetStdout(PETSC_COMM_SELF,&viewer); CHKERRQ(ierr);
  }
  ierr = PetscTypeCompare((PetscObject)viewer,PETSC_VIEWER_ASCII,&iascii); CHKERRQ(ierr);
  if (!iascii) { SETERRQ(1,"Only ASCII viewer for ColumnSystem\n"); }

  if (L==NULL) {
    SETERRQ1(2,"columnSystem ERROR: can't matrix '%s' because L==NULL ... ending ...\n", info);
  }
  if (D==NULL) {
    SETERRQ1(3,"columnSystem ERROR: can't matrix '%s' because D==NULL ... ending ...\n", info);
  }
  if (U==NULL) {
    SETERRQ1(4,"columnSystem ERROR: can't matrix '%s' because U==NULL ... ending ...\n", info);
  }

  if (nmax < 2) {
    ierr = PetscViewerASCIIPrintf(viewer,
      "\n\n<nmax >= 2 required to view ColumnSystem tridiagonal matrix '%s' ... skipping view\n",info);
    CHKERRQ(ierr);
    return 0;
  }

//  if (nmax > 12) {
  if (nmax > 120) {
    ierr = PetscViewerASCIIPrintf(viewer,
      "\n\n<nmax > 120: matrix too big to display as full; viewing tridiagonal matrix '%s' by diagonals ...\n",info); CHKERRQ(ierr);
    char diag_info[PETSC_MAX_PATH_LEN];
    snprintf(diag_info,PETSC_MAX_PATH_LEN, "super-diagonal U for system '%s'", info);
    ierr = viewColumnValues(viewer,U,nmax-1,diag_info); CHKERRQ(ierr);
    snprintf(diag_info,PETSC_MAX_PATH_LEN, "main diagonal D for system '%s'", info);
    ierr = viewColumnValues(viewer,D,nmax,diag_info); CHKERRQ(ierr);
    snprintf(diag_info,PETSC_MAX_PATH_LEN, "sub-diagonal L for system '%s'", info);
    ierr = viewColumnValues(viewer,Lp,nmax-1,diag_info); CHKERRQ(ierr);
  } else {
    ierr = PetscViewerASCIIPrintf(viewer,
        "\n<viewing ColumnSystem tridiagonal matrix, '%s':\n",info); CHKERRQ(ierr);
    for (PetscInt k=0; k<nmax; k++) {    // k+1 is row  (while j+1 is column)
      if (k == 0) {              // viewing first row
        ierr = PetscViewerASCIIPrintf(viewer,"%10.7f %10.7f ",D[k],U[k]); CHKERRQ(ierr);
        for (PetscInt j=2; j<nmax; j++) {
          ierr = PetscViewerASCIIPrintf(viewer,"%10.7f ",0.0); CHKERRQ(ierr);
        }
      } else if (k < nmax-1) {   // viewing generic row
        for (PetscInt j=0; j<k-1; j++) {
          ierr = PetscViewerASCIIPrintf(viewer,"%10.7f ",0.0); CHKERRQ(ierr);
        }
        ierr = PetscViewerASCIIPrintf(viewer,"%10.7f %10.7f %10.7f ",L[k],D[k],U[k]); CHKERRQ(ierr);
        for (PetscInt j=k+2; j<nmax; j++) {
          ierr = PetscViewerASCIIPrintf(viewer,"%10.7f ",0.0); CHKERRQ(ierr);
        }
      } else {                   // viewing last row
        for (PetscInt j=0; j<k-1; j++) {
          ierr = PetscViewerASCIIPrintf(viewer,"%10.7f ",0.0); CHKERRQ(ierr);
        }
        ierr = PetscViewerASCIIPrintf(viewer,"%10.7f %10.7f ",L[k],D[k]); CHKERRQ(ierr);
      }
      ierr = PetscViewerASCIIPrintf(viewer,"\n"); CHKERRQ(ierr);  // end of row
    }
  }
  
  ierr = PetscViewerASCIIPrintf(viewer,
      "end viewing tridiagonal matrix>\n"); CHKERRQ(ierr);
  return 0;
}


//! View the tridiagonal system A x = b, both A as a full matrix and b as a vector.  Unwise if size exceeds 20 or so. 
/*! 
Calls viewMatrix() to view A and viewColumnValues() to view right-hand-side b.
 */
PetscErrorCode columnSystemCtx::viewSystem(PetscViewer viewer, const char* info) const {
  PetscErrorCode ierr;
  ierr = viewMatrix(viewer,info); CHKERRQ(ierr);
  char  rhs_info[PETSC_MAX_PATH_LEN];
  snprintf(rhs_info,PETSC_MAX_PATH_LEN, "right-hand side vector for system '%s'", info);
  ierr = viewColumnValues(viewer,rhs,nmax,rhs_info); CHKERRQ(ierr);
  return 0;
}


//! The actual code for solving a tridiagonal system.  Return code has diagnostic importance.
/*!
This is modified slightly from a Numerical Recipes version.

Input size n is size of instance.  Requires n <= columnSystemCtx::nmax.

Solution of system in x.

Success is return code zero.  Positive return code gives location of zero pivot.
Negative return code indicates a software problem.
 */
PetscErrorCode columnSystemCtx::solveTridiagonalSystem(
                  const PetscInt n, PetscScalar **x) {
  if (x == NULL) { SETERRQ(-999,"x is NULL in columnSystemCtx"); }
  if (*x == NULL) { SETERRQ(-998,"*x is NULL in columnSystemCtx"); }
  if (n < 1) { SETERRQ(-997,"instance size n < 1 in columnSystemCtx"); }
  if (n > nmax) { SETERRQ(-996,"instance size n too large in columnSystemCtx"); }

  if (!indicesValid) { SETERRQ(-995,"column indices not valid in columnSystemCtx"); }

  PetscScalar b;
  b = D[0];
  if (b == 0.0) { return 1; }
  (*x)[0] = rhs[0]/b;
  for (int i=1; i<n; ++i) {
    work[i] = U[i-1]/b;
    b = D[i] - L[i] * work[i];
    if (b == 0.0) { return i+1; }
    (*x)[i] = (rhs[i] - L[i] * (*x)[i-1]) / b;
  }
  for (int i=n-2; i>=0; --i) {
    (*x)[i] -= work[i+1] * (*x)[i+1];
  }

  indicesValid = false;
  return 0;
}


ageSystemCtx::ageSystemCtx(PetscInt my_Mz)
      : columnSystemCtx(my_Mz) { // size of system is Mz
  initAllDone = false;
  // set values so we can check if init was called on all
  dx = -1.0;
  dy = -1.0;
  dtAge = -1.0;
  dzEQ = -1.0;
  u = NULL;
  v = NULL;
  w = NULL;
  tau3 = NULL;
}


PetscErrorCode ageSystemCtx::initAllColumns() {
  // check whether each parameter & pointer got set
  if (dx <= 0.0) { SETERRQ(2,"un-initialized dx in ageSystemCtx"); }
  if (dy <= 0.0) { SETERRQ(3,"un-initialized dy in ageSystemCtx"); }
  if (dtAge <= 0.0) { SETERRQ(4,"un-initialized dtAge in ageSystemCtx"); }
  if (dzEQ <= 0.0) { SETERRQ(5,"un-initialized dzEQ in ageSystemCtx"); }
  if (u == NULL) { SETERRQ(6,"un-initialized pointer u in ageSystemCtx"); }
  if (v == NULL) { SETERRQ(7,"un-initialized pointer v in ageSystemCtx"); }
  if (w == NULL) { SETERRQ(8,"un-initialized pointer w in ageSystemCtx"); }
  if (tau3 == NULL) { SETERRQ(9,"un-initialized pointer tau3 in ageSystemCtx"); }
  nuEQ = dtAge / dzEQ; // derived constant
  initAllDone = true;
  return 0;
}


PetscErrorCode ageSystemCtx::solveThisColumn(PetscScalar **x) {
  PetscErrorCode ierr;
  if (!initAllDone) {  SETERRQ(2,
     "solveThisColumn() should only be called after initAllColumns() in ageSystemCtx"); }

  // set up system: 0 <= k < ks
  for (PetscInt k = 0; k < ks; k++) {
    planeStar ss;  // note ss.ij = tau[k]
    ierr = tau3->getPlaneStarZ(i,j,k * dzEQ,&ss); CHKERRQ(ierr);
    // do lowest-order upwinding, explicitly for horizontal
    rhs[k] =  (u[k] < 0) ? u[k] * (ss.ip1 -  ss.ij) / dx
                         : u[k] * (ss.ij  - ss.im1) / dx;
    rhs[k] += (v[k] < 0) ? v[k] * (ss.jp1 -  ss.ij) / dy
                         : v[k] * (ss.ij  - ss.jm1) / dy;
    // note it is the age eqn: dage/dt = 1.0 and we have moved the hor.
    //   advection terms over to right:
    rhs[k] = ss.ij + dtAge * (1.0 - rhs[k]);

    // do lowest-order upwinding, *implicitly* for vertical
    PetscScalar AA = nuEQ * w[k];
    if (k > 0) {
      if (AA >= 0) { // upward velocity
        L[k] = - AA;
        D[k] = 1.0 + AA;
        U[k] = 0.0;
      } else { // downward velocity; note  -AA >= 0
        L[k] = 0.0;
        D[k] = 1.0 - AA;
        U[k] = + AA;
      }
    } else { // k == 0 case
      // note L[0] not an allocated location
      if (AA > 0) { // if strictly upward velocity apply boundary condition:
                    // age = 0 because ice is being added to base
        D[0] = 1.0;
        U[0] = 0.0;
        rhs[0] = 0.0;
      } else { // downward velocity; note  -AA >= 0
        D[0] = 1.0 - AA;
        U[0] = + AA;
        // keep rhs[0] as is
      }
    }
  }  // done "set up system: 0 <= k < ks"
      
  // surface b.c. at ks
  if (ks>0) {
    L[ks] = 0;
    D[ks] = 1.0;   // ignore U[ks]
    rhs[ks] = 0.0;  // age zero at surface
  }

  // solve it
  return solveTridiagonalSystem(ks+1,x);
}


tempSystemCtx::tempSystemCtx(PetscInt my_Mz, PetscInt my_Mbz)
      : columnSystemCtx(my_Mz + my_Mbz - 1) {
  Mz = my_Mz;
  Mbz = my_Mbz;
  k0 = Mbz - 1; // max size nmax of system is Mz + k0 = Mz + Mbz - 1
  // set flags to indicate nothing yet set
  initAllDone = false;
  schemeParamsValid = false;
  surfBCsValid = false;
  basalBCsValid = false;
  // set values so we can check if init was called on all
  dx = -1;
  dy = -1;
  dtTemp = -1;
  dzEQ = -1;
  dzbEQ = -1;
  ice_rho = -1;
  ice_c_p = -1;
  ice_k = -1;
  bed_thermal_rho = -1;
  bed_thermal_c_p = -1;
  bed_thermal_k = -1;
  T = NULL;
  Tb = NULL;
  u = NULL;
  v = NULL;
  w = NULL;
  Sigma = NULL;
  T3 = NULL;
}


PetscErrorCode tempSystemCtx::initAllColumns() {
  // check whether each parameter & pointer got set
  if (dx <= 0.0) { SETERRQ(2,"un-initialized dx in tempSystemCtx"); }
  if (dy <= 0.0) { SETERRQ(3,"un-initialized dy in tempSystemCtx"); }
  if (dtTemp <= 0.0) { SETERRQ(4,"un-initialized dtTemp in tempSystemCtx"); }
  if (dzEQ <= 0.0) { SETERRQ(5,"un-initialized dzEQ in tempSystemCtx"); }
  if (dzbEQ <= 0.0) { SETERRQ(6,"un-initialized dzbEQ in tempSystemCtx"); }
  if (ice_rho <= 0.0) { SETERRQ(7,"un-initialized ice_rho in tempSystemCtx"); }
  if (ice_c_p <= 0.0) { SETERRQ(8,"un-initialized ice_c_p in tempSystemCtx"); }
  if (ice_k <= 0.0) { SETERRQ(9,"un-initialized ice_k in tempSystemCtx"); }
  if (bed_thermal_rho <= 0.0) { SETERRQ(10,"un-initialized bed_thermal_rho in tempSystemCtx"); }
  if (bed_thermal_c_p <= 0.0) { SETERRQ(11,"un-initialized bed_thermal_c_p in tempSystemCtx"); }
  if (bed_thermal_k <= 0.0) { SETERRQ(12,"un-initialized bed_thermal_k in tempSystemCtx"); }
  if (T == NULL) { SETERRQ(13,"un-initialized pointer T in tempSystemCtx"); }
  if (Tb == NULL) { SETERRQ(14,"un-initialized pointer Tb in tempSystemCtx"); }
  if (u == NULL) { SETERRQ(15,"un-initialized pointer u in tempSystemCtx"); }
  if (v == NULL) { SETERRQ(16,"un-initialized pointer v in tempSystemCtx"); }
  if (w == NULL) { SETERRQ(17,"un-initialized pointer w in tempSystemCtx"); }
  if (Sigma == NULL) { SETERRQ(18,"un-initialized pointer Sigma in tempSystemCtx"); }
  if (T3 == NULL) { SETERRQ(19,"un-initialized pointer T3 in tempSystemCtx"); }
  // set derived constants
  nuEQ = dtTemp / dzEQ;
  rho_c_I = ice_rho * ice_c_p;
  rho_c_br = bed_thermal_rho * bed_thermal_c_p;
  rho_c_av = (dzEQ * rho_c_I + dzbEQ * rho_c_br) / (dzEQ + dzbEQ);
  iceK = ice_k / rho_c_I;
  iceR = iceK * dtTemp / PetscSqr(dzEQ);
  brK = bed_thermal_k / rho_c_br;
  brR = brK * dtTemp / PetscSqr(dzbEQ);
  rho_c_ratio = rho_c_I / rho_c_av;
  dzav = 0.5 * (dzEQ + dzbEQ);
  iceReff = ice_k * dtTemp / (rho_c_av * dzEQ * dzEQ);
  brReff = bed_thermal_k * dtTemp / (rho_c_av * dzbEQ * dzbEQ);
  // done
  initAllDone = true;
  return 0;
}


PetscErrorCode tempSystemCtx::setSchemeParamsThisColumn(
                     PismMask my_mask, bool my_isMarginal, PetscScalar my_lambda) {
  if (!initAllDone) {  SETERRQ(2,
     "setSchemeParamsThisColumn() should only be called after initAllColumns() in tempSystemCtx"); }
  if (schemeParamsValid) {  SETERRQ(3,
     "setSchemeParamsThisColumn() called twice (?) in tempSystemCtx"); }
  mask = my_mask;
  isMarginal = my_isMarginal;
  lambda = my_lambda;
  schemeParamsValid = true;
  return 0;
}


PetscErrorCode tempSystemCtx::setSurfaceBoundaryValuesThisColumn(PetscScalar my_Ts) {
  if (!initAllDone) {  SETERRQ(2,
     "setSurfaceBoundaryValuesThisColumn() should only be called after initAllColumns() in tempSystemCtx"); }
  if (surfBCsValid) {  SETERRQ(3,
     "setSurfaceBoundaryValuesThisColumn() called twice (?) in tempSystemCtx"); }
  Ts = my_Ts;
  surfBCsValid = true;
  return 0;
}


PetscErrorCode tempSystemCtx::setBasalBoundaryValuesThisColumn(
                     PetscScalar my_Ghf, PetscScalar my_Tshelfbase, PetscScalar my_Rb) {
  if (!initAllDone) {  SETERRQ(2,
     "setIndicesThisColumn() should only be called after initAllColumns() in tempSystemCtx"); }
  if (basalBCsValid) {  SETERRQ(3,
     "setBasalBoundaryValuesThisColumn() called twice (?) in tempSystemCtx"); }
  Ghf = my_Ghf;
  Tshelfbase = my_Tshelfbase;
  Rb = my_Rb;
  basalBCsValid = true;
  return 0;
}


PetscErrorCode tempSystemCtx::solveThisColumn(PetscScalar **x) {
  PetscErrorCode ierr;
  if (!initAllDone) {  SETERRQ(2,
     "solveThisColumn() should only be called after initAllColumns() in tempSystemCtx"); }
  if (!schemeParamsValid) {  SETERRQ(3,
     "solveThisColumn() should only be called after setSchemeParamsThisColumn() in tempSystemCtx"); }
  if (!surfBCsValid) {  SETERRQ(3,
     "solveThisColumn() should only be called after setSurfaceBoundaryValuesThisColumn() in tempSystemCtx"); }
  if (!basalBCsValid) {  SETERRQ(3,
     "solveThisColumn() should only be called after setBasalBoundaryValuesThisColumn() in tempSystemCtx"); }

  if (Mbz > 1) { // bedrock present: build k=0:Mbz-2 eqns
    // gives O(\Delta t,\Delta z^2) convergence in Test K for equal spaced grid;
    // note L[0] not an allocated location:
    D[0] = (1.0 + 2.0 * brR);
    U[0] = - 2.0 * brR;  
    rhs[0] = Tb[0] + 2.0 * dtTemp * Ghf / (rho_c_br * dzbEQ);

    // bedrock only; pure vertical conduction problem
    for (PetscInt k=1; k < k0; k++) {
      L[k] = -brR;
      D[k] = 1.0 + 2.0 * brR;
      U[k] = -brR;
      rhs[k] = Tb[k];
    }
  }

  // bottom part of ice (and top of bedrock in some cases): k=k0=Mbz-1 eqn
  if (ks == 0) { // no ice; set T[0] to surface temp if grounded
    if (k0 > 0) { L[k0] = 0.0; } // note L[0] not allocated 
    D[k0] = 1.0;
    U[k0] = 0.0;
    // if floating and no ice then worry only about bedrock temps
    if (mask >= MASK_FLOATING) {
      // essentially no ice but floating ... ask PISMOceanCoupler
      rhs[k0] = Tshelfbase;
      // FIXME: split k0 into two grid points?
    } else { // top of bedrock sees atmosphere
      rhs[k0] = Ts; 
    }
  } else { // ks > 0; there is ice
    // for w, always difference *up* from base, but make it implicit
    if (mask >= MASK_FLOATING) {
      // just apply Dirichlet condition to base of column of ice in an ice shelf
      if (k0 > 0) { L[k0] = 0.0; } // note L[0] not allocated 
      D[k0] = 1.0;
      U[k0] = 0.0;
      rhs[k0] = Tshelfbase; // set by PISMOceanCoupler
    } else { 
      // there is *grounded* ice; ice/bedrock interface; from FV across interface
      rhs[k0] = T[0] + dtTemp * (Rb / (rho_c_av * dzav));
      if (!isMarginal) {
        rhs[k0] += dtTemp * rho_c_ratio * 0.5 * (Sigma[0] / rho_c_I);
        // WARNING: subtle consequences of finite volume argument across interface
        planeStar ss;
        ierr = T3->getPlaneStarZ(i,j,0.0,&ss);
        const PetscScalar UpTu = (u[0] < 0) ? u[0] * (ss.ip1 -  ss.ij) / dx :
                                              u[0] * (ss.ij  - ss.im1) / dx;
        const PetscScalar UpTv = (v[0] < 0) ? v[0] * (ss.jp1 -  ss.ij) / dy :
                                              v[0] * (ss.ij  - ss.jm1) / dy;
        rhs[k0] -= dtTemp * rho_c_ratio * (0.5 * (UpTu + UpTv));
      }
      const PetscScalar AA = dtTemp * rho_c_ratio * w[0] / (2.0 * dzEQ);
      if (Mbz > 1) { // there is bedrock; apply upwinding if w[0]<0,
                     // otherwise ignore advection; note 
                     // jump in diffusivity coefficient
        L[k0] = - brReff;
        if (w[0] >= 0.0) {  // velocity upward
          D[k0] = 1.0 + iceReff + brReff;
          U[k0] = - iceReff;
        } else { // velocity downward
          D[k0] = 1.0 + iceReff + brReff - AA;
          U[k0] = - iceReff + AA;
        }
      } else { // no bedrock; apply geothermal flux here
        // L[k0] = 0.0;  (note this is not an allocated location!) 
        if (w[0] >= 0.0) {  // velocity upward
          D[k0] = 1.0 + 2.0 * iceR;
          U[k0] = - 2.0 * iceR;
        } else { // velocity downward
          D[k0] = 1.0 + 2.0 * iceR - AA;
          U[k0] = - 2.0 * iceR + AA;
        }
        rhs[k0] += 2.0 * dtTemp * Ghf / (rho_c_I * dzEQ);
      }
    }
  }

  // generic ice segment: build k0+1:k0+ks-1 eqns
  for (PetscInt k = 1; k < ks; k++) {
    planeStar ss;
    ierr = T3->getPlaneStarZ(i,j,k * dzEQ,&ss);
    const PetscScalar UpTu = (u[k] < 0) ? u[k] * (ss.ip1 -  ss.ij) / dx :
                                          u[k] * (ss.ij  - ss.im1) / dx;
    const PetscScalar UpTv = (v[k] < 0) ? v[k] * (ss.jp1 -  ss.ij) / dy :
                                          v[k] * (ss.ij  - ss.jm1) / dy;
    const PetscScalar AA = nuEQ * w[k];      
    if (w[k] >= 0.0) {  // velocity upward
      L[k0+k] = - iceR - AA * (1.0 - lambda/2.0);
      D[k0+k] = 1.0 + 2.0 * iceR + AA * (1.0 - lambda);
      U[k0+k] = - iceR + AA * (lambda/2.0);
    } else {  // velocity downward
      L[k0+k] = - iceR - AA * (lambda/2.0);
      D[k0+k] = 1.0 + 2.0 * iceR - AA * (1.0 - lambda);
      U[k0+k] = - iceR + AA * (1.0 - lambda/2.0);
    }
    rhs[k0+k] = T[k];
    if (!isMarginal) {
      rhs[k0+k] += dtTemp * (Sigma[k] / rho_c_I - UpTu - UpTv);
    }
  }
      
  // surface b.c.
  if (ks>0) {
    L[k0+ks] = 0.0;
    D[k0+ks] = 1.0;
    // ignore U[k0+ks]
    rhs[k0+ks] = Ts;
  }

  // mark column as done
  schemeParamsValid = false;
  surfBCsValid = false;
  basalBCsValid = false;

  // solve it; note melting not addressed yet
  return solveTridiagonalSystem(k0+ks+1,x);
}

