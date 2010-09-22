// Copyright (C) 2004-2010, 2010 Jed Brown, Ed Bueler and Constantine Khroulev
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
#include "../base/pism_const.hh"
#include "../base/materials.hh"
#include "../base/NCVariable.hh"

static char help[] =
  "Calls IceFlowLaw::flow...() with various values of arguments and prints results (for software tests).\n";

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

    EnthalpyConverter EC(config);

    IceFlowLaw *ice = NULL;
    IceFlowLawFactory ice_factory(com, NULL, config);

    string ice_type = ICE_GPBLD;
    ice_factory.setType(ICE_GPBLD); // set the default type

    ierr = ice_factory.setFromOptions(); CHKERRQ(ierr);
    ice_factory.create(&ice);

    bool dummy;
    ierr = PISMOptionsString("-ice_type", "Selects the flow law", ice_type, dummy); CHKERRQ(ierr);

    double     E0=560000, dE=10000, p=2e7, T = 0, gs = 1e-3;
    double     sigma[] = {1e4, 5e4, 1e5, 1.5e5};

    printf("flow law:   \"%s\"\n", ice_type.c_str());
    printf("flowtable:  [pressure = %10.2e throughout]\n",p);
    printf("  (stress)   (enthalpy)    (temp)     =   (flow)\n");

    PolyThermalGPBLDIce *poly_ice;
    if (ice_type == "gpbld") {
      poly_ice = dynamic_cast<PolyThermalGPBLDIce*>(ice);
      if (poly_ice == NULL) {
        PetscPrintf(com,
          "ERROR: poly_ice == NULL after dynamic cast, when given '-ice_type gpbld'\n");
        PetscEnd();
      }
    } else
      poly_ice = NULL;

    for (int i=0; i<4; ++i) {
      for (int j=0; j<5; ++j) {
        double E = E0 - j*dE;
        EC.getAbsTemp(E, p, T);
        if (ice_type == "gpbld") {
          printf("%10.2e   %10.3f    %10.6f = %10.6e\n",
                 sigma[i], E, T, poly_ice->flow_from_enth(sigma[i], E, p, gs));
        } else {
          printf("%10.2e   %10.3f    %10.6f = %10.6e\n",
                 sigma[i], E, T, ice->flow_from_temp(sigma[i], T, p, gs));
        }
      }
    }

    delete ice;
  } // end explicit scope

  ierr = PetscFinalize(); CHKERRQ(ierr);
  return 0;
}
