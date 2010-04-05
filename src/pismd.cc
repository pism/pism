// Copyright (C) 2007--2010 Ed Bueler and Constantine Khroulev
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

#include <petsc.h>
#include "base/grid.hh"
#include "base/iceModel.hh"
#include "eismint/iceROSSModel.hh"

#include "coupler/PCFactory.hh"
#include "coupler/PISMAtmosphere.hh"
#include "coupler/PISMSurface.hh"
#include "coupler/PISMOcean.hh"

static char help[] =
  "Driver for EISMINT-Ross diagnostic velocity computation in ice shelf.)\n";

int main(int argc, char *argv[]) {
  PetscErrorCode  ierr;

  MPI_Comm    com;
  PetscMPIInt rank, size;
  ierr = PetscInitialize(&argc, &argv, PETSC_NULL, help); CHKERRQ(ierr);
  com = PETSC_COMM_WORLD;
  ierr = MPI_Comm_rank(com, &rank); CHKERRQ(ierr);
  ierr = MPI_Comm_size(com, &size); CHKERRQ(ierr);

  { /* This explicit scoping forces destructors to be called before PetscFinalize() */
    ierr = verbosityLevelFromOptions(); CHKERRQ(ierr);

    ierr = verbPrintf(2,com, "PISMROSS %s (EISMINT-Ross diagnostic velocity computation mode)\n",
		      PISM_Revision); CHKERRQ(ierr);
    ierr = stop_on_version_option(); CHKERRQ(ierr);

    bool iset, bfset;
    ierr = PISMOptionsIsSet("-i", iset); CHKERRQ(ierr);
    ierr = PISMOptionsIsSet("-boot_from", bfset); CHKERRQ(ierr);
    string usage =
      "  pismross {-i IN.nc|-boot_from IN.nc} [OTHER PISM & PETSc OPTIONS]\n\n"
      "where:\n"
      "  -i          input file in NetCDF format: contains PISM-written model state\n"
      "  -boot_from  input file in NetCDF format: contains a few fields, from which\n"
      "              heuristics will build initial model state\n"
      "  -ssaBC      read SSA boundary conditions from a file\n"
      "  -riggs      read RIGGS data from a file\n"
      "notes:\n"
      "  * one of -i or -boot_from is required\n"
      "  * if -boot_from is used then in fact '-Mx A -My B -Mz C -Lz D' is also required\n";
    if ((iset == PETSC_FALSE) && (bfset == PETSC_FALSE)) {
      ierr = PetscPrintf(com,
         "PISM ERROR: one of options -i,-boot_from is required\n\n"); CHKERRQ(ierr);
      ierr = show_usage_and_quit(com, "pismross", usage.c_str()); CHKERRQ(ierr);
    } else {
      vector<string> required;  required.clear();
      ierr = show_usage_check_req_opts(com, "pismross", required, usage.c_str()); CHKERRQ(ierr);
    }


    NCConfigVariable config, overrides;
    ierr = init_config(com, rank, config, overrides); CHKERRQ(ierr);

    config.set_flag("use_ssa_velocity", true);
    config.set_flag("use_ssa_when_grounded", false);
    config.set_flag("use_constant_nuh_for_ssa", false);
    config.set_flag("do_mass_conserve", false);
    config.set_flag("do_temp", false);
    config.set_flag("do_cold_ice_methods", true);
    config.set("epsilon_ssa", 0.0);  // don't use this lower bound on effective viscosity
    config.set("run_length_years", 0);
    config.set_flag("force_full_diagnostics", true);

    IceGrid    g(com, rank, size, config);

    // Initialize boundary models:
    PAFactory pa(g, config);
    PISMAtmosphereModel *atmosphere;

    PSFactory ps(g, config);
    PISMSurfaceModel *surface;
    ierr = ps.set_default("constant"); CHKERRQ(ierr);

    POFactory po(g, config);
    PISMOceanModel *ocean;

    pa.create(atmosphere);
    ps.create(surface);
    po.create(ocean);

    surface->attach_atmosphere_model(atmosphere);

    IceROSSModel m(g, config, overrides);

    m.attach_surface_model(surface);
    m.attach_ocean_model(ocean);

    ierr = m.setExecName("pismross"); CHKERRQ(ierr);

    ierr = m.init(); CHKERRQ(ierr);

    ierr = verbPrintf(2,com, "computing velocity field (diagnostically) ...\n"); CHKERRQ(ierr);
    ierr = m.run(); CHKERRQ(ierr);
    ierr = verbPrintf(2,com, "... done\n"); CHKERRQ(ierr);

    ierr = m.finishROSS(); CHKERRQ(ierr);

    ierr = m.writeFiles("unnamed_diag.nc"); CHKERRQ(ierr);  // default filename if no -o
    
  }
  ierr = PetscFinalize(); CHKERRQ(ierr);
  return 0;
}
