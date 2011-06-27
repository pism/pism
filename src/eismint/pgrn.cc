// Copyright (C) 2007-2011 Ed Bueler and Nathan Shemonski and Constantine Khroulev
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

//! \file pgrn.cc Implements EISMINT-Greenland experiments.
/*! \file pgrn.cc
This driver adds only the minimum functionality needed to implement the
choices stated in [\ref RitzEISMINT], the EISMINT-Greenland specification:
- A PDD melt model is used by default.
- Non-standard PDD parameters are used.
- PA_EISMINT_Greenland atmosphere model implements EISMINT-Greenland air temperature parameterization.
- An enhancement factor of 3.0 is used.
- -ocean_kill is used by default.

Which experiment to do is chosen by one of options -ssl2,-ccl3,-gwl3.  (Experiment
SSL3 is not implemented; see User's Manual.)
 */

static char help[] = 
  "Driver for PISM runs of EISMINT-Greenland intercomparison.\n";

#include <petsc.h>
#include "iceModel.hh"

#include "PCFactory.hh"
#include "PAEismintGreenland.hh"
#include "PISMSurface.hh"
#include "PISMOcean.hh"


static PetscErrorCode set_eismint_greenland_params(MPI_Comm com,
						   NCConfigVariable &config) {
  PetscErrorCode ierr;

  bool ssl2Set,ssl3Set;
  // we never care about default "-ssl2", but this avoids "Option left" message:
  ierr = PISMOptionsIsSet("-ssl2", ssl2Set); CHKERRQ(ierr); 
  ierr = PISMOptionsIsSet("-ssl3", ssl3Set); CHKERRQ(ierr);
  if (ssl3Set) {
    PetscPrintf(com,
       "PISM ERROR:  experiment SSL3 (-ssl3) is not implemented ... ENDING\n"
       "  (implement by choosing config parameters and runtime options to your taste)\n");
    PISMEnd();
  }

  ierr = verbPrintf(2, com, 
     "  setting flags equivalent to '-e 3 -ocean_kill'; user options may override ...\n");
     CHKERRQ(ierr);
  config.set("enhancement_factor", 3.0);
  config.set_flag("ocean_kill", true);
  // use the EISMINT-Greenland value if no value in -boot_file file
  config.set("bootstrapping_geothermal_flux_value_no_var", 0.050);

  bool ccl3Set, gwl3Set;
  ierr = PISMOptionsIsSet("-ccl3", ccl3Set); CHKERRQ(ierr);
  ierr = PISMOptionsIsSet("-gwl3", gwl3Set); CHKERRQ(ierr);
  if (ccl3Set || gwl3Set) {
    // use Lingle-Clark bed deformation model for CCL3 and GWL3 but not SSL2
    ierr = verbPrintf(2, com, 
      "  setting flags equivalent to: '-bed_def lc'; user options may override ...\n");
      CHKERRQ(ierr);
    config.set_string("bed_deformation_model", "lc");
  }
  bool gwl3_start_set;
  ierr = PISMOptionsIsSet("-gwl3_start_year", gwl3_start_set); CHKERRQ(ierr);
  if (gwl3_start_set && !gwl3Set) {
    PetscPrintf(com,
       "PISM ERROR: option -gwl3_start_year is allowed only if -gwl3 is set.\n");
    PISMEnd();
  }

  // degree-day factors in \ref RitzEISMINT are water-equivalent thickness
  // per degree day; ice-equivalent thickness melted per degree day is
  // slightly larger; for example, iwfactor = 1000/910
  const PetscReal iwfactor = config.get("fresh_water_density") / config.get("ice_density");
  PetscReal pdd_factor_ice = 0.008 * iwfactor, // convert from water- to ice-equivalent
    pdd_factor_snow        = 0.003 * iwfactor, // ditto
    pdd_refreeze           = 0.6,
    pdd_std_dev            = 5.0;

  ierr = verbPrintf(2, com,
		    "Setting PDD parameters to EISMINT-Greenland values...\n"
		    "  pdd_factor_ice  = %3.5f m (ice equivalent) K-1 day-1\n"
		    "  pdd_factor_snow = %3.5f m (ice equivalent) K-1 day-1\n"
		    "  pdd_refreeze    = %3.5f 1\n"
		    "  pdd_std_dev     = %3.5f K\n",
		    pdd_factor_ice, pdd_factor_snow, pdd_refreeze, pdd_std_dev);
  config.set("pdd_factor_ice",  pdd_factor_ice);
  config.set("pdd_factor_snow", pdd_factor_snow);
  config.set("pdd_refreeze",    pdd_refreeze);
  config.set("pdd_std_dev",     pdd_std_dev);

  return 0;
}

int main(int argc, char *argv[]){
  PetscErrorCode ierr;

  MPI_Comm com;
  PetscMPIInt rank, size;
  
  ierr = PetscInitialize(&argc, &argv, PETSC_NULL, help); CHKERRQ(ierr);

  com = PETSC_COMM_WORLD;
  ierr = MPI_Comm_rank(com, &rank); CHKERRQ(ierr);
  ierr = MPI_Comm_size(com, &size); CHKERRQ(ierr);
  
  { // explicit scoping does destructors before PetscFinalize() 
    ierr = verbosityLevelFromOptions(); CHKERRQ(ierr);

    ierr = verbPrintf(2, com, "PGRN %s (PISM EISMINT-Greenland mode)\n",
		      PISM_Revision); CHKERRQ(ierr);
    ierr = stop_on_version_option(); CHKERRQ(ierr);

    bool iset, bfset;
    ierr = PISMOptionsIsSet("-i", iset); CHKERRQ(ierr);
    ierr = PISMOptionsIsSet("-boot_file", bfset); CHKERRQ(ierr);
    string usage =
      "  pgrn {-i IN.nc|-boot_file IN.nc} [OTHER PISM & PETSc OPTIONS]\n"
      "where:\n"
      "  -i          IN.nc is input file in NetCDF format: contains PISM-written model state\n"
      "  -boot_file  IN.nc is input file in NetCDF format: contains a few fields, from which\n"
      "              heuristics will build initial model state\n"
      "and one of these EISMINT-Greenland experiments [-ssl2 is default]:\n"
      "  -ssl2, -ccl3, -gwl3\n"
      "notes:\n"
      "  * pgrn is a special executable for EISMINT-Greenland\n"
      "  * one of -i or -boot_file is required\n"
      "  * if -boot_file is used then '-Mx A -My B -Mz C -Lz D' is also required\n"
      "  * generally behaves like pismr after initialization\n";
    if ((iset == PETSC_FALSE) && (bfset == PETSC_FALSE)) {
      ierr = PetscPrintf(com,
         "PISM ERROR: one of options -i,-boot_file is required\n\n"); CHKERRQ(ierr);
      ierr = show_usage_and_quit(com, "pgrn", usage.c_str()); CHKERRQ(ierr);
    } else {
      vector<string> required;  required.clear();
      ierr = show_usage_check_req_opts(com, "pgrn", required, usage.c_str()); CHKERRQ(ierr);
    }

    NCConfigVariable config, overrides;
    ierr = init_config(com, rank, config, overrides); CHKERRQ(ierr);

    ierr = set_eismint_greenland_params(com, config); CHKERRQ(ierr);

    IceGrid  g(com, rank, size, config);
    IceModel m(g, config, overrides);
    ierr = m.setExecName("pgrn"); CHKERRQ(ierr);

    // boundary models:
    PAFactory pa(g, config);
    PISMAtmosphereModel *atmosphere;
    ierr = pa.set_default("eismint_greenland"); CHKERRQ(ierr);

    PSFactory ps(g, config);
    PISMSurfaceModel *surface;
    ierr = ps.set_default("pdd"); CHKERRQ(ierr);

    POFactory po(g, config);
    PISMOceanModel *ocean;

    ierr = PetscOptionsBegin(com, "", "Options choosing PISM boundary models", ""); CHKERRQ(ierr);
    pa.create(atmosphere);
    ps.create(surface);
    po.create(ocean);
    ierr = PetscOptionsEnd(); CHKERRQ(ierr);

    surface->attach_atmosphere_model(atmosphere);

    m.attach_surface_model(surface);
    m.attach_ocean_model(ocean);
 
    ierr = m.init(); CHKERRQ(ierr);
 
    ierr = m.run(); CHKERRQ(ierr);
    ierr = verbPrintf(2, com, "done with run ... \n"); CHKERRQ(ierr);

    ierr = m.writeFiles("grn_exper.nc"); CHKERRQ(ierr);
  }

  ierr = PetscFinalize(); CHKERRQ(ierr);
  return 0;
}

