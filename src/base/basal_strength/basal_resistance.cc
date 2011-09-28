// Copyright (C) 2004-2011 Jed Brown, Ed Bueler, and Constantine Khroulev
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

#include "basal_resistance.hh"
#include "pism_const.hh"
#include "enthalpyConverter.hh"

IceBasalResistancePlasticLaw::IceBasalResistancePlasticLaw(
             PetscScalar regularizationConstant, bool pseudoPlastic,
             PetscScalar pseudoExponent, PetscScalar pseudoUThreshold) {
  plastic_regularize = regularizationConstant;
  pseudo_plastic = pseudoPlastic;
  pseudo_q = pseudoExponent;
  pseudo_u_threshold = pseudoUThreshold;
}

PetscErrorCode IceBasalResistancePlasticLaw::printInfo(int verbthresh, MPI_Comm com) {
  PetscErrorCode ierr;
  if (pseudo_plastic == PETSC_TRUE) {
    if (pseudo_q == 1.0) {
      ierr = verbPrintf(verbthresh, com, 
        "Using linearly viscous till with u_threshold = %.2f m/a.\n", 
                        convert(pseudo_u_threshold, "m/s", "m/year")); CHKERRQ(ierr);
    } else {
      ierr = verbPrintf(verbthresh, com, 
        "Using pseudo-plastic till with eps = %10.5e m/a, q = %.4f,"
        " and u_threshold = %.2f m/a.\n", 
                        convert(plastic_regularize, "m/s", "m/year"),
                        pseudo_q,
                        convert(pseudo_u_threshold, "m/s", "m/year")); 
        CHKERRQ(ierr);
    }
  } else {
    ierr = verbPrintf(verbthresh, com, 
      "Using purely plastic till with eps = %10.5e m/a.\n",
                      convert(plastic_regularize, "m/s", "m/year")); CHKERRQ(ierr);
  }
  return 0;
}


//! Compute the drag coefficient for the basal shear stress.
/*!
The basal shear stress term \f$\tau_b\f$ in the SSA stress balance for ice
is minus the return value here times (vx,vy).  Thus this method computes the
basal shear stress as
    \f[ \tau_b = - \frac{\tau_c}{|\mathbf{U}|^{1-q} U_{\mathtt{th}}^q} \mathbf{U} \f]
where \f$\tau_b=(\tau_{(b)x},\tau_{(b)y})\f$, \f$U=(u,v)\f$,
\f$q=\f$ <tt>pseudo_q</tt>, and \f$U_{\mathtt{th}}=\f$ <tt>pseudo_u_threshold</tt>.
Typical values for the constants are \f$q=0.25\f$ and \f$U_{\mathtt{th}} = 100\f$
m/a.

The linearly-viscous till case pseudo_q = 1.0 is allowed, in which case 
\f$\beta = \tau_c/U_{\mathtt{th}}\f$.  The purely-plastic till case pseudo_q = 0.0
is also allowed; note that there is still a regularization with data member
plastic_regularize.
 */
PetscScalar IceBasalResistancePlasticLaw::drag(PetscScalar tauc,
                                   PetscScalar vx, PetscScalar vy) {
  const PetscScalar magreg2 = PetscSqr(plastic_regularize) + PetscSqr(vx) + PetscSqr(vy);
  if (pseudo_plastic == PETSC_TRUE) {
    return tauc * pow(magreg2, 0.5*(pseudo_q - 1)) * pow(pseudo_u_threshold, -pseudo_q);
  } else { // pure plastic, but regularized
    return tauc / sqrt(magreg2);
  }
}

// Derivative of drag with respect to \f$ \alpha = \frac 1 2 (u_x^2 + u_y^2) \f$
void IceBasalResistancePlasticLaw::dragWithDerivative(PetscReal tauc, PetscScalar vx, PetscScalar vy,
						      PetscScalar *d, PetscScalar *dd) const
{
  const PetscScalar magreg2 = PetscSqr(plastic_regularize) + PetscSqr(vx) + PetscSqr(vy);
  if (pseudo_plastic == PETSC_TRUE) {
    *d = tauc * pow(magreg2, 0.5*(pseudo_q - 1)) * pow(pseudo_u_threshold, -pseudo_q);
    if (dd) *dd = (pseudo_q - 1) * *d / magreg2;
  } else { // pure plastic, but regularized
    *d = tauc / sqrt(magreg2);
    if (dd) *dd = -1 * *d / magreg2;
  }
}
