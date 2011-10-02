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

#include "varcEnthalpyConverter.hh"

double EnthalpyConverter::getEnthalpyCTS(double p) const {
  const double
     T        = getMeltingTemp(p),
     T_r      = 256.82,
     Trefdiff = 0.5 * (T + T_0) - T_r; 
  return c_i * (T - T_0) + 7.253 * Trefdiff * (T - T_0);
}


PetscErrorCode EnthalpyConverter::getAbsTemp(double E, double p, double &T) const {
  double E_s, E_l;
  PetscErrorCode ierr = getEnthalpyInterval(p, E_s, E_l); CHKERRQ(ierr);
  if (E >= E_l) { // enthalpy equals or exceeds that of liquid water
    T = getMeltingTemp(p);
    return 1;
  }
  if (E < E_s) {
    T = (E / c_i) + T_0;  FIXME
  } else {
    T = getMeltingTemp(p);
  }
  return 0;
}


PetscErrorCode varcEnthalpyConverter::getEnth(
                  double T, double omega, double p, double &E) const {
  const double T_m = getMeltingTemp(p);
  if (T <= 0.0) {
    SETERRQ1(1,"\n\nT = %f <= 0 is not a valid absolute temperature\n\n",T);
  }
  if ((omega < 0.0 - 1.0e-6) || (1.0 + 1.0e-6 < omega)) {
    SETERRQ1(2,"\n\nwater fraction omega=%f not in range [0,1]\n\n",omega);
  }
  if (T > T_m + 1.0e-6) {
    SETERRQ2(3,"T=%f exceeds T_m=%f; not allowed\n\n",T,T_m);
  }
  if ((T < T_m - 1.0e-6) && (omega > 0.0 + 1.0e-6)) {
    SETERRQ3(4,"T < T_m AND omega > 0 is contradictory\n\n",T,T_m,omega);
  }
  if (T < T_m) {
    const double
       T_r      = 256.82,
       Trefdiff = 0.5 * (T + T_0) - T_r; 
    E = c_i * (T - T_0) + 7.253 * Trefdiff * (T - T_0);
  } else {
    E = getEnthalpyCTS(p) + omega * L;
  }
  return 0;
}

