// Copyright (C) 2004-2011 Jed Brown, Ed Bueler and Constantine Khroulev
//
// This file is part of Pism.
//
// Pism is free software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the Free Software
// Foundation; either version 3 of the License, or (at your option) any later
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
#include "flowlaw_factory.hh"
#include "NCVariable.hh"
#include "enthalpyConverter.hh"
#include "pism_options.hh"

static char help[] =
  "Calls IceFlowLaw with various values of arguments and prints results.\n"
  "Used for software tests.  Tests the flow() method but prints\n"
  "temperature and liquid fraction as inputs and flow coefficient as output.\n"
  "Thus also tests methods getPressureFromDepth(), getMeltingTemp(), and\n"
  "getEnth() methods of EnthalpyConverter.  Nonetheless a change to the\n"
  "enthalpy normalization only should not affect the outcome.  Only physically-\n"
  "meaningful inputs and output appear at stdout.\n";

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

    IceFlowLaw *flow_law = NULL;
    IceFlowLawFactory ice_factory(com, NULL, config, &EC);

    string flow_law_name = ICE_GPBLD;
    ice_factory.setType(ICE_GPBLD); // set the default type

    ierr = ice_factory.setFromOptions(); CHKERRQ(ierr);
    ice_factory.create(&flow_law);

    bool dummy;
    ierr = PISMOptionsString("-flow_law", "Selects the flow law",
                             flow_law_name, dummy); CHKERRQ(ierr);

    double     TpaC[]  = {-30.0, -5.0, 0.0, 0.0},  // pressure-adjusted, deg C
               depth   = 2000.0,
               gs      = 1.0e-3, // some laws use grain size; fixed
               omega0  = 0.005,  // some laws use liquid fraction; used w TpaC[3]
               sigma[] = {1e4, 5e4, 1e5, 1.5e5};

    double     p       = EC.getPressureFromDepth(depth),
               Tm      = EC.getMeltingTemp(p);

    printf("flow law:   \"%s\"\n", flow_law_name.c_str());
    printf("pressure = %9.3e Pa = (hydrostatic at depth %7.2f m)\n",
           p,depth);
    printf("flowtable:\n");
    printf("  (dev stress)   (abs temp) (liq frac) =   (flow)\n");

    for (int i=0; i<4; ++i) {
      for (int j=0; j<4; ++j) {

        double T     = Tm + TpaC[j],
               omega = (j == 3) ? omega0 : 0.0;

        double E, flowcoeff;
        EC.getEnth(T, omega, p, E);
        flowcoeff = flow_law->flow(sigma[i], E, p, gs);

        printf("    %10.2e   %10.3f  %9.3f = %10.6e\n",
               sigma[i], T, omega, flowcoeff);

      }
    }

    delete flow_law;
  } // end explicit scope

  ierr = PetscFinalize(); CHKERRQ(ierr);
  return 0;
}
