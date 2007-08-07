// Copyright (C) 2004-2007 Jed Brown and Ed Bueler
//
// This file is part of Pism.
//
// Pism is free software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the Free Software
// Foundation; either version 2 of the License, or (at your option) any later
// version.
//
// Pism is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
// details.
//
// You should have received a copy of the GNU General Public License
// along with Pism; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

#include <petsc.h>
#include "base/materials.hh"

int main(int argc, char *argv[]) {
  PetscScalar     T0=273.15, dT=10, p=2e7;
  PetscScalar     sigma[] = {1e4, 5e4, 1e5, 1.5e5};
  ThermoGlenIce   tglen;
  HybridIce       hyb;
  GKparts         pt;

  printf("flowtable:  [pressure = %10.2e in whole table]\n",p);
  printf("  (stress)   (temp)   =   (flow)\n");
  printf("THERMOGLENICE:\n");
  for (int i=0; i<4; ++i) {
    for (int j=0; j<5; ++j) {
      PetscScalar T = T0 - j*dT;
      printf("%10.2e %10.3f = %10.2e\n", sigma[i], T, tglen.flow(sigma[i], T, p));
    }
  }

  printf("HYBRIDICE:  [after (flow) are four parts: (diff, basal, gbs, disl)]\n");
  for (int i=0; i<4; ++i) {
    for (int j=0; j<5; ++j) {
      PetscScalar T = T0 - j*dT;
      pt=hyb.flowParts(sigma[i], T, p);
      printf("%10.2e %10.3f = %10.2e = (%8.2e, %8.2e, %8.2e, %8.2e)\n",
             sigma[i], T, pt.eps_total,pt.eps_diff,pt.eps_basal,pt.eps_gbs,pt.eps_disl);
    }
  }

  return 0;
}
