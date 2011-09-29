// Copyright (C) 2011 Ed Bueler
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

#include "enthalpyConverter.hh"
#include "varkenthSystem.hh"
#include <gsl/gsl_math.h>


varkenthSystemCtx::varkenthSystemCtx(
    const NCConfigVariable &config, IceModelVec3 &my_Enth3, int my_Mz, string my_prefix)
      : enthSystemCtx(config, my_Enth3, my_Mz, my_prefix) {
      
  EC = new EnthalpyConverter(config);
}


varkenthSystemCtx::~varkenthSystemCtx() {
  delete EC;
}



PetscErrorCode varkenthSystemCtx::viewConstants(
                     PetscViewer viewer, bool show_col_dependent) {
  PetscErrorCode ierr;

  if (!viewer) {
    ierr = PetscViewerASCIIGetStdout(PETSC_COMM_SELF,&viewer); CHKERRQ(ierr);
  }

  PetscTruth iascii;
  ierr = PetscTypeCompare((PetscObject)viewer,PETSC_VIEWER_ASCII,&iascii); CHKERRQ(ierr);
  if (!iascii) { SETERRQ(1,"Only ASCII viewer for varkenthSystemCtx::viewConstants()\n"); }
  
  ierr = PetscViewerASCIIPrintf(viewer,
                   "\n<<varkenthSystemCtx IS A MODIFICATION OF enthSystemCtx>>\n"); CHKERRQ(ierr);

  ierr = enthSystemCtx::viewConstants(viewer, show_col_dependent); CHKERRQ(ierr);
  return 0;
}


// FIXME:  Neumann boundary condition should be fixed to use getvark() also,
//         but the base of the ice is where the k(T) value is close to k_i anyway

/*!  Equation (4.37) in \ref GreveBlatter2009 is
  \f[ k(T ) = 9.828 e^{−0.0057 T} \f]
where \f$T\f$ is in Kelvin and the resulting conductivity is in units W m−1 K−1.
 */
PetscScalar varkenthSystemCtx::getvark(PetscScalar T) {
  return 9.828 * exp(-0.0057 * T);
}


PetscErrorCode varkenthSystemCtx::solveThisColumn(PetscScalar **x, PetscErrorCode &pivoterrorindex) {
  PetscErrorCode ierr;
#ifdef PISM_DEBUG
  ierr = checkReadyToSolve(); CHKERRQ(ierr);
  if ((gsl_isnan(a0)) || (gsl_isnan(a1)) || (gsl_isnan(b))) {
    SETERRQ(1, "solveThisColumn() should only be called after\n"
               "  setting basal boundary condition in enthSystemCtx"); }
#endif
  // k=0 equation is already established
  // L[0] = 0.0;  // not allocated
  D[0] = a0;
  U[0] = a1;
  rhs[0] = b;

  const PetscScalar Rfactor = dtTemp / (PetscSqr(dzEQ) * ice_rho * ice_c);

  // generic ice segment in k location (if any; only runs if ks >= 2)
  for (PetscInt k = 1; k < ks; k++) {
    // FIXME:  below commits O(dz) error in depth because ks * dzEQ is not exactly the thickness
    // FIXME:  below is inefficient because k(T) at grid points are (expensively)
    //         recomputed instead of being stored in a column
    PetscScalar Tl, Tc, Tr;
    ierr = EC->getAbsTemp(Enth[k-1], EC->getPressureFromDepth((ks - (k-1))*dzEQ), Tl); CHKERRQ(ierr); 
    ierr = EC->getAbsTemp(Enth[k],   EC->getPressureFromDepth((ks - ( k ))*dzEQ), Tc); CHKERRQ(ierr); 
    ierr = EC->getAbsTemp(Enth[k+1], EC->getPressureFromDepth((ks - (k+1))*dzEQ), Tr); CHKERRQ(ierr); 
    const PetscScalar
        Rl = (Enth[k-1] < Enth_s[k-1]) ? (getvark(Tl) * Rfactor) : iceRtemp,
        Rc = (Enth[k]   < Enth_s[k]  ) ? (getvark(Tc) * Rfactor) : iceRtemp,
        Rr = (Enth[k+1] < Enth_s[k+1]) ? (getvark(Tr) * Rfactor) : iceRtemp,
        Rminus = 0.5 * (Rl + Rc),
        Rplus  = 0.5 * (Rc + Rr);
    L[k] = - Rminus;
    D[k] = 1.0 + Rminus + Rplus;
    U[k] = - Rplus;
    const PetscScalar AA = nuEQ * w[k];
    if (w[k] >= 0.0) {  // velocity upward
      L[k] -= AA * (1.0 - lambda/2.0);
      D[k] += AA * (1.0 - lambda);
      U[k] += AA * (lambda/2.0);
    } else {            // velocity downward
      L[k] -= AA * (lambda/2.0);
      D[k] -= AA * (1.0 - lambda);
      U[k] += AA * (1.0 - lambda/2.0);
    }
    rhs[k] = Enth[k];
    if (!ismarginal) {
      planeStar<PetscScalar> ss;
      ierr = Enth3->getPlaneStar_fine(i,j,k,&ss); CHKERRQ(ierr);
      const PetscScalar UpEnthu = (u[k] < 0) ? u[k] * (ss.e -  ss.ij) / dx :
                                               u[k] * (ss.ij  - ss.w) / dx;
      const PetscScalar UpEnthv = (v[k] < 0) ? v[k] * (ss.n -  ss.ij) / dy :
                                               v[k] * (ss.ij  - ss.s) / dy;
      rhs[k] += dtTemp * ((Sigma[k] / ice_rho) - UpEnthu - UpEnthv);
    }
  }

  // set Dirichlet boundary condition at top
  if (ks > 0) L[ks] = 0.0;
  D[ks] = 1.0;
  if (ks < Mz-1) U[ks] = 0.0;
  rhs[ks] = Enth_ks;

  // solve it; note drainage is not addressed yet and post-processing may occur
  pivoterrorindex = solveTridiagonalSystem(ks+1, x);

  // air above
  for (PetscInt k = ks+1; k < Mz; k++) {
    (*x)[k] = Enth_ks;
  }

#ifdef PISM_DEBUG
  if (pivoterrorindex == 0) {
    // if success, mark column as done by making scheme params and b.c. coeffs invalid
    lambda  = -1.0;
    a0 = GSL_NAN;
    a1 = GSL_NAN;
    b  = GSL_NAN;
  }
#endif
  return 0;
}

