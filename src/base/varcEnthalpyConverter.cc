// Copyright (C) 2011, 2012 The PISM Authors
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


#include "pism_const.hh"
#include "varcEnthalpyConverter.hh"


/*!
A calculation only used in the cold case.

Let \f$c_i=2009.0\,\, \text{J}\,\text{kg}^{-1}\,\text{K}^{-1}\f$, the default
constant value.  Note equation (4.39) in [\ref GreveBlatter2009] says 
\f$C(T) = c_i + 7.253 (T - T_r)\f$ using a reference temperature
\f$T_r = 256.81786846822\f$ K.  Thus the calculation of enthalpy from cold ice
temperature is: 
\f{align*}{
  E(T) &= \int_{T_0}^T C(T')\,dT \\
       &= c_i (T-T_0) + \frac{7.253}{2} \left((T-T_r)^2 - (T_0-T_r)^2\right) \\
       &= \left(c_i + 7.253 \left(\frac{T+T_0}{2} - T_r\right)\right) (T- T_0).
\f}
 */
double varcEnthalpyConverter::EfromT(double T) const {
  const double
     Trefdiff = 0.5 * (T + T_0) - T_r; 
  return (c_i + c_gradient * Trefdiff) * (T - T_0);
}


/*!
A calculation only used in the cold case.

From the documentation for EfromT(), in the cold ice case we must solve the equation
  \f[ E = \left(c_i + 7.253 \left(\frac{T+T_0}{2} - T_r\right)\right) (T- T_0) \f]
for \f$T\f$.  This equation is quadratic in \f$T\f$.  If we write it in terms
of \f$\Delta T = T-T_0\f$ and \f$\alpha = 2/7.253\f$ then it says
  \f[ E = c_i \Delta T + \alpha^{-1} \left(\Delta T + 2(T_0-T_r)\right) \Delta T.\f]
Rearranging as a standard form quadratic in unknown \f$\Delta T\f$ gives
  \f[ 0 = \Delta T^2 + \left[\alpha c_i + 2(T_0-T_r)\right] \Delta T - \alpha E.\f]
Define
  \f[ \beta = \alpha c_i + 2(T_0-T_r) = 486.64, \f]
which has units K; the value comes from \f$c_i=2009\f$, \f$T_0=223.15\f$, and
\f$T_r=256.81786846822\f$.  Then the solution of the quadratic is the one which makes 
\f$\Delta T \ge 0\f$ assuming \f$E\ge 0\f$; we stop otherwise.  With the usual
rewriting to avoid cancellation we have
  \f[ \Delta T = \frac{-\beta + \sqrt{\beta^2 + 4 \alpha E}}{2} = \frac{2 \alpha E}{\sqrt{\beta^2 + 4 \alpha E} + \beta}.\f]
Of course, \f$T=T_0 + \Delta T\f$.
 */
double varcEnthalpyConverter::TfromE(double E) const {
  if (E < 0.0) {
    PetscPrintf(PETSC_COMM_WORLD,"\n\nE < 0 in varcEnthalpyConverter is not allowed.  FIXME.\n\n");
    PISMEnd();
  }
  const double
    ALPHA = 2.0 / c_gradient,
    BETA  = ALPHA * c_i + 2.0 * (T_0 - T_r),
    tmp   = 2.0 * ALPHA * E,
    dT    = tmp / (sqrt(BETA*BETA + 2.0*tmp) + BETA);
  return T_0 + dT;
}


PetscErrorCode varcEnthalpyConverter::viewConstants(PetscViewer viewer) const {
  PetscErrorCode ierr;

  PetscBool iascii;
  if (!viewer) {
    ierr = PetscViewerASCIIGetStdout(PETSC_COMM_WORLD,&viewer); CHKERRQ(ierr);
  }
  ierr = PetscObjectTypeCompare((PetscObject)viewer,PETSCVIEWERASCII,&iascii); CHKERRQ(ierr);
  if (!iascii) { SETERRQ(PETSC_COMM_SELF, 1,"Only ASCII viewer for EnthalpyConverter\n"); }

  ierr = PetscViewerASCIIPrintf(viewer,
    "\n<class varcEnthalpyConverter has two additional constants, so as to implement\n"
       "  equation (4.39) from Greve & Blatter (2009):"); CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,
      "   T_r   = %12.5f (K)\n",         T_r); CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,
      "   c_gradient = %12.8f  >\n",    c_gradient); CHKERRQ(ierr);

  ierr = EnthalpyConverter::viewConstants(viewer); CHKERRQ(ierr);
  return 0;
}

//! Redefined from EnthalpyConverter version, for use when specific heat capacity depends on temperature.
/*!
Calls EfromT().
 */
double varcEnthalpyConverter::getEnthalpyCTS(double p) const {
  return EfromT(getMeltingTemp(p));
}


//! Redefined from EnthalpyConverter version, for use when specific heat capacity depends on temperature.
/*!
Calls TfromE().
 */
PetscErrorCode varcEnthalpyConverter::getAbsTemp(double E, double p, double &T) const {
  double E_s, E_l;
  PetscErrorCode ierr = getEnthalpyInterval(p, E_s, E_l); CHKERRQ(ierr);
  if (E >= E_l) { // enthalpy equals or exceeds that of liquid water
    T = getMeltingTemp(p);
    return 1;
  }
  if (E < E_s) {
    T = TfromE(E);
  } else {
    T = getMeltingTemp(p);
  }
  return 0;
}


//! Redefined from EnthalpyConverter version, for use when specific heat capacity depends on temperature.
/*!
Calls EfromT().
 */
PetscErrorCode varcEnthalpyConverter::getEnth(
                  double T, double omega, double p, double &E) const {
  const double T_m = getMeltingTemp(p);
  if (T <= 0.0) {
    SETERRQ1(PETSC_COMM_SELF, 1,"\n\nT = %f <= 0 is not a valid absolute temperature\n\n",T);
  }
  if ((omega < 0.0 - 1.0e-6) || (1.0 + 1.0e-6 < omega)) {
    SETERRQ1(PETSC_COMM_SELF, 2,"\n\nwater fraction omega=%f not in range [0,1]\n\n",omega);
  }
  if (T > T_m + 1.0e-6) {
    SETERRQ2(PETSC_COMM_SELF, 3,"T=%f exceeds T_m=%f; not allowed\n\n",T,T_m);
  }
  if ((T < T_m - 1.0e-6) && (omega > 0.0 + 1.0e-6)) {
    SETERRQ3(PETSC_COMM_SELF, 4,"T < T_m AND omega > 0 is contradictory\n\n",T,T_m,omega);
  }
  if (T < T_m) {
    E = EfromT(T);
  } else {
    E = getEnthalpyCTS(p) + omega * L;
  }
  return 0;
}

