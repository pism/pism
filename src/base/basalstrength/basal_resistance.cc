// Copyright (C) 2004-2012 Jed Brown, Ed Bueler, and Constantine Khroulev
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

/* Purely plastic */

IceBasalResistancePlasticLaw::IceBasalResistancePlasticLaw(const NCConfigVariable &config) {
  plastic_regularize = config.get("plastic_regularization", "m/year", "m/second");
}

PetscErrorCode IceBasalResistancePlasticLaw::printInfo(int verbthresh, MPI_Comm com) {
  PetscErrorCode ierr;
  ierr = verbPrintf(verbthresh, com, 
                    "Using purely plastic till with eps = %10.5e m/a.\n",
                    convert(plastic_regularize, "m/s", "m/year")); CHKERRQ(ierr);

  return 0;
}


//! Compute the drag coefficient for the basal shear stress.
PetscScalar IceBasalResistancePlasticLaw::drag(PetscScalar tauc,
                                               PetscScalar vx, PetscScalar vy) {
  const PetscScalar magreg2 = PetscSqr(plastic_regularize) + PetscSqr(vx) + PetscSqr(vy);

  return tauc / sqrt(magreg2);
}

//! Compute the drag coefficient and its derivative with respect to \f$ \alpha = \frac 1 2 (u_x^2 + u_y^2) \f$
void IceBasalResistancePlasticLaw::dragWithDerivative(PetscReal tauc, PetscScalar vx, PetscScalar vy,
						      PetscScalar *d, PetscScalar *dd) const
{
  const PetscScalar magreg2 = PetscSqr(plastic_regularize) + PetscSqr(vx) + PetscSqr(vy);

  *d = tauc / sqrt(magreg2);

  if (dd)
    *dd = -1 * (*d) / magreg2;

}

/* Pseudo-plastic */

IceBasalResistancePseudoPlasticLaw::IceBasalResistancePseudoPlasticLaw(const NCConfigVariable &config)
  : IceBasalResistancePlasticLaw(config) {
  pseudo_q = config.get("pseudo_plastic_q");
  pseudo_u_threshold = config.get("pseudo_plastic_uthreshold", "m/year", "m/second");
  sliding_scale = config.get("sliding_scale_factor_reduces_tauc");
}

PetscErrorCode IceBasalResistancePseudoPlasticLaw::printInfo(int verbthresh, MPI_Comm com) {
  PetscErrorCode ierr;

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
/*! One can scale tauc if desired:

A scale factor of \f$A\f$ is intended to increase basal sliding rate by
\f$A\f$.  It would have exactly this effect \e if the driving stress were
\e hypothetically completely held by the basal resistance.  Thus this scale factor
is used to reduce (if \c -sliding_scale \f$A\f$ with \f$A > 1\f$) or increase
(if \f$A < 1\f$) the value of the (pseudo-) yield stress \c tauc.  The concept
behind this is described at
http://websrv.cs.umt.edu/isis/index.php/Category_1:_Whole_Ice_Sheet#Initial_Experiment_-_E1_-_Increased_Basal_Lubrication.

Specifically, the concept behind this mechanism is to suppose equality of driving
and basal shear stresses,
    \f[ \rho g H \nabla h = \frac{\tau_c}{|\mathbf{U}|^{1-q} U_{\mathtt{th}}^q} \mathbf{U}. \f]
(<i>For emphasis:</i> The membrane stress held by the ice itself is missing from
this incomplete stress balance.)  Thus the pseudo yield stress
\f$\tau_c\f$ would be related to the sliding speed \f$|\mathbf{U}|\f$ by
  \f[ |\mathbf{U}| = \frac{C}{\tau_c^{1/q}} \f]
for some (geometry-dependent) constant \f$C\f$.  Multiplying \f$|\mathbf{U}|\f$
by \f$A\f$ in this equation corresponds to dividing \f$\tau_c\f$ by \f$A^q\f$.

  Note that the mechanism has no effect whatsoever if
\f$q=0\f$, which is the purely plastic case. In that case there is \e no direct
relation between the yield stress and the sliding velocity, and the difference
between the driving stress and the yield stress is entirely held by the membrane
stresses.  (There is also no singular mathematical operation as \f$A^q = A^0 = 1\f$.)
*/
PetscScalar IceBasalResistancePseudoPlasticLaw::drag(PetscScalar tauc,
                                                     PetscScalar vx, PetscScalar vy) {
  const PetscScalar magreg2 = PetscSqr(plastic_regularize) + PetscSqr(vx) + PetscSqr(vy);

  if (sliding_scale > 0.0) {
    double Aq = pow(sliding_scale, pseudo_q);
    return (tauc / Aq) * pow(magreg2, 0.5*(pseudo_q - 1)) * pow(pseudo_u_threshold, -pseudo_q);
  } else {
    return tauc * pow(magreg2, 0.5*(pseudo_q - 1)) * pow(pseudo_u_threshold, -pseudo_q);
  }
}


//! Compute the drag coefficient and its derivative with respect to \f$ \alpha = \frac 1 2 (u_x^2 + u_y^2) \f$
void IceBasalResistancePseudoPlasticLaw::dragWithDerivative(PetscReal tauc, PetscScalar vx, PetscScalar vy,
                                                            PetscScalar *d, PetscScalar *dd) const
{
  const PetscScalar magreg2 = PetscSqr(plastic_regularize) + PetscSqr(vx) + PetscSqr(vy);

  if (sliding_scale > 0.0) {
    double Aq = pow(sliding_scale, pseudo_q);
    *d = (tauc / Aq) * pow(magreg2, 0.5*(pseudo_q - 1)) * pow(pseudo_u_threshold, -pseudo_q);
  } else {
    *d =  tauc * pow(magreg2, 0.5*(pseudo_q - 1)) * pow(pseudo_u_threshold, -pseudo_q);
  }

  if (dd)
    *dd = (pseudo_q - 1) * (*d) / magreg2;

}
