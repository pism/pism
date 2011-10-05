// Copyright (C) 2011 The PISM Authors
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
#include "varcEnthalpyConverter.hh"
#include "varkenthSystem.hh"
#include <gsl/gsl_math.h>


varkenthSystemCtx::varkenthSystemCtx(const NCConfigVariable &config,
                                     IceModelVec3 &my_Enth3, int my_Mz,
                                     string my_prefix, EnthalpyConverter *my_EC)
  : enthSystemCtx(config, my_Enth3, my_Mz, my_prefix), EC(my_EC) {

  if (dynamic_cast<varcEnthalpyConverter*>(EC) != NULL)
    use_variable_c = true;
  else
    use_variable_c = false;
}

PetscErrorCode varkenthSystemCtx::initAllColumns(PetscScalar my_dx, PetscScalar my_dy,
                                                 PetscScalar my_dtTemp, PetscScalar my_dzEQ) {
  PetscErrorCode ierr = enthSystemCtx::initAllColumns(my_dx, my_dy, my_dtTemp, my_dzEQ); CHKERRQ(ierr);

  for (PetscInt k = 0; k < Mz; k++)
    R[k] = iceRcold;  // fill with cold constant value for safety

  return 0;
}


PetscErrorCode varkenthSystemCtx::viewConstants(PetscViewer viewer,
                                                bool show_col_dependent) {
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


/*!  Equation (4.37) in \ref GreveBlatter2009 is
  \f[ k(T ) = 9.828 e^{−0.0057 T} \f]
  where \f$T\f$ is in Kelvin and the resulting conductivity is in units W m−1 K−1.
 */
PetscScalar varkenthSystemCtx::getvark(PetscScalar T) {
  return 9.828 * exp(-0.0057 * T);
}

/*!
  Equation (4.39) in [\ref GreveBlatter2009] is
  \f$C(T) = c_i + 7.253 (T - T_r)\f$, with a reference temperature
  \f$T_r = 256.82\f$ K.

  FIXME: not sure if this method should be here or in varcEnthalpyConverter...
 */
PetscScalar varkenthSystemCtx::getvarc(PetscScalar T) {
  return ice_c + 7.253 * (T - 256.82);
}

//! \brief Assemble the R array.
PetscErrorCode varkenthSystemCtx::assemble_R() {
  PetscErrorCode ierr;

  const PetscScalar Rfactor = dtTemp / (PetscSqr(dzEQ) * ice_rho);

  for (PetscInt k = 0; k <= ks; k++) {
    if (Enth[k] < Enth_s[k]) {
      // cold case
      const PetscScalar depth = (ks - k) * dzEQ;
      // FIXME: commits O(dz) error in "depth" because ks * dzEQ is not
      // exactly the thickness
      PetscScalar T;
      ierr = EC->getAbsTemp(Enth[k], EC->getPressureFromDepth(depth), // FIXME: task #7297
                            T); CHKERRQ(ierr);

      if (use_variable_c)
        R[k] = getvark(T) / getvarc(T) * Rfactor;
      else
        R[k] = getvark(T) / ice_c * Rfactor;
    } else {
      // temperate case
      R[k] = iceRtemp;
    }
  }

  // R[k] for k > ks are never used
#ifdef PISM_DEBUG
  for (int k = ks + 1; k < Mz; ++k)
    R[k] = GSL_NAN;
#endif

  return 0;
}

