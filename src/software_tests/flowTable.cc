// Copyright (C) 2004-2007, 2010, 2011 Jed Brown and Ed Bueler
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
#include "pism_const.hh"
#include "flowlaws.hh"
#include "NCVariable.hh"
#include "enthalpyConverter.hh"
#include "pism_options.hh"


static char help[] =
  "Show a table of flow results from the Goldsby-Kohlstedt law (HYBRIDICE),\n"
  "  compared to results from Paterson-Budd (THERMOGLENICE).\n";

int main(int argc, char *argv[]) {
  PetscErrorCode  ierr;

  MPI_Comm    com;
  PetscMPIInt rank, size;

  ierr = PetscInitialize(&argc, &argv, PETSC_NULL, help); CHKERRQ(ierr);

  com = PETSC_COMM_WORLD;
  ierr = MPI_Comm_rank(com, &rank); CHKERRQ(ierr);
  ierr = MPI_Comm_size(com, &size); CHKERRQ(ierr);

  /* This explicit scoping forces destructors to be called before PetscFinalize() */
  {
  NCConfigVariable config, overrides;
  ierr = init_config(com, rank, config, overrides); CHKERRQ(ierr);

  ierr = verbosityLevelFromOptions(); CHKERRQ(ierr);
  ierr = verbPrintf(2,com, 
     "FLOWTABLE %s (test program for ThermoGlenIce,HybridIce,EnthalpyConverter classes)\n",
     PISM_Revision); CHKERRQ(ierr);
  ierr = stop_on_version_option(); CHKERRQ(ierr);

  EnthalpyConverter EC(config);

  double     T0=273.15, dT=10, p=2e7;
  double     sigma[] = {1e4, 5e4, 1e5, 1.5e5};
  // Ice constructors need a communicator and a prefix, since this is a serial code, the communicator should not be used
  // anyway so we should be okay just passing 0 for the communicator (note that the ABI for MPI_Comm is not defined, it
  // is typedef'd to int in MPICH and to an opaque pointer in OpenMPI, the constant 0 is implicitly cast to a pointer if
  // need be.
  ThermoGlenIce   tglen(0, NULL, config, &EC);
  HybridIce       hyb(0, NULL, config, &EC);
  GKparts         pt;

  printf("flowtable:  [pressure = %10.2e in whole table]\n",p);
  printf("  (stress)   (temp)   =   (flow)\n");
  printf("THERMOGLENICE:\n");
  for (int i=0; i<4; ++i) {
    for (int j=0; j<5; ++j) {
      double T = T0 - j*dT;
      printf("%10.2e %10.3f = %10.2e\n", sigma[i], T, tglen.flow_from_temp(sigma[i], T, p, 0));
    }
  }

  printf("HYBRIDICE:  [after (flow) are four parts: (diff, basal, gbs, disl)]\n");
  for (int i=0; i<4; ++i) {
    for (int j=0; j<5; ++j) {
      double T = T0 - j*dT;
      pt=hyb.flowParts(sigma[i], T, p);
      printf("%10.2e %10.3f = %10.2e = (%8.2e, %8.2e, %8.2e, %8.2e)\n",
             sigma[i], T, pt.eps_total,pt.eps_diff,pt.eps_basal,pt.eps_gbs,pt.eps_disl);
    }
  }

  } // end explicit scope

  ierr = PetscFinalize(); CHKERRQ(ierr);
  return 0;
}
