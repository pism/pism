// Copyright (C) 2010 Ed Bueler and Constantine Khroulev
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

// TO DO:  make std out from velocitySSA() immediate by using PetscPrintf
// (and then run on greenland example again)

/*

Example uses, from a PISM directory:

$ pismv -test I -My 25 -Mx 3 -o testI.nc
$ ssa_test -i testI.nc

$ pismv -test J -Mx 40 -My 40 -o testJ.nc
$ ssa_test -i testJ.nc

$ cd example/pst/
$ ./pst.sh 8 >> out.pst &  # takes a while; produces P1.nc; can be shortened
$ ssa_test -i P1.nc


*/

static char help[] =
  "\nSSA_TEST\n"
  "  Testing program for SSA, time-independent calculations separate from\n"
  "  IceModel.  Also may be used in a PISM software (regression) test.\n\n";

#include <cmath>
#include <cstdio>
#include <string>
#include <petscksp.h>
#include "pism_const.hh"
#include "grid.hh"
#include "iceModelVec.hh"
#include "flowlaw_factory.hh" // IceFlowLawFactory, IceFlowLaw
#include "materials.hh" // SSAStrengthExtension, IceBasalResistancePlasticLaw
#include "nc_util.hh"
#include "PISMIO.hh"
#include "NCVariable.hh"
#include "SSAFD.hh"

int main(int argc, char *argv[]) {
  PetscErrorCode  ierr;

  MPI_Comm    com;  // won't be used except for rank,size
  PetscMPIInt rank, size;

  ierr = PetscInitialize(&argc, &argv, PETSC_NULL, help); CHKERRQ(ierr);

  com = PETSC_COMM_WORLD;
  ierr = MPI_Comm_rank(com, &rank); CHKERRQ(ierr);
  ierr = MPI_Comm_size(com, &size); CHKERRQ(ierr);
  
  /* This explicit scoping forces destructors to be called before PetscFinalize() */
  {  
    NCConfigVariable config, overrides;
    ierr = init_config(com, rank, config, overrides); CHKERRQ(ierr);
    //    config.set_flag("write_ssa_system_to_matlab", true);

    PetscTruth i_set, usage_set, help_set;
    ierr = PetscOptionsHasName(PETSC_NULL, "-usage", &usage_set); CHKERRQ(ierr);
    ierr = PetscOptionsHasName(PETSC_NULL, "-help", &help_set); CHKERRQ(ierr);
    if ((usage_set==PETSC_TRUE) || (help_set==PETSC_TRUE)) {
      PetscPrintf(com,
        "\nusage of SSA_TEST:\n"
          "  1) create a PISM output file foo.nc with variables thk, usurf, bed, tauc\n"
          "  2) do 'ssa_test -i foo.nc' or\n"
          "  3) or do 'mpiexec -n 2 ssa_test -i foo.nc -display :0'\n\n");
    }
    
    // get input file name and open
    char filename[PETSC_MAX_PATH_LEN];
    ierr = PetscOptionsGetString(PETSC_NULL, "-i", filename, 
                                 PETSC_MAX_PATH_LEN, &i_set); CHKERRQ(ierr);
    if (i_set==PETSC_FALSE) {
      PetscPrintf(com,
        "\nSSAFD_TEST ERROR:  requires PISM-written NetCDF file as input ... ending!\n\n");
      PetscEnd();
    }

    IceGrid grid(com, rank, size, config);
    ierr = PetscPrintf(grid.com, 
        "SSAFD_TEST\n  initializing from NetCDF file '%s'...\n", filename); CHKERRQ(ierr);
    PISMIO pio(&grid);
    ierr = pio.get_grid(filename); CHKERRQ(ierr); // fails if filename not present
    grid.start_year = grid.year;

    grid.compute_nprocs();
    grid.compute_ownership_ranges();
    ierr = grid.compute_horizontal_spacing(); CHKERRQ(ierr);
    ierr = grid.createDA(); CHKERRQ(ierr);

    ierr = setVerbosityLevel(5); CHKERRQ(ierr);
    ierr = grid.printInfo(1); CHKERRQ(ierr);
    //ierr = grid.printVertLevels(1); CHKERRQ(ierr); 

    IceFlowLaw *ice = NULL;
    IceFlowLawFactory ice_factory(com, NULL, config);
    string ice_type = ICE_GPBLD;
    ice_factory.setType(ICE_GPBLD); // set the default type
    ierr = ice_factory.setFromOptions(); CHKERRQ(ierr);
    ice_factory.create(&ice);
    //ierr = ice->printInfo(1); CHKERRQ(ierr);

    IceBasalResistancePlasticLaw basal(
           config.get("plastic_regularization") / secpera, 
           config.get_flag("do_pseudo_plastic_till"),
           config.get("pseudo_plastic_q"),
           config.get("pseudo_plastic_uthreshold") / secpera);
    //ierr = basal.printInfo(1,grid.com); CHKERRQ(ierr);

    EnthalpyConverter EC(config);

    IceModelVec2S vh, vH, vbed, vtauc;
    IceModelVec2Mask vMask;
    IceModelVec3 enthalpy;
    const PetscInt WIDE_STENCIL = 2;

    PISMVars vars;

    // ice upper surface elevation
    ierr = vh.create(grid, "usurf", true, WIDE_STENCIL); CHKERRQ(ierr);
    ierr = vh.set_attrs("diagnostic", "ice upper surface elevation",
          "m", "surface_altitude"); CHKERRQ(ierr);
    ierr = vars.add(vh); CHKERRQ(ierr);

    // land ice thickness
    ierr = vH.create(grid, "thk", true, WIDE_STENCIL); CHKERRQ(ierr);
    ierr = vH.set_attrs("model_state", "land ice thickness",
          "m", "land_ice_thickness"); CHKERRQ(ierr);
    ierr = vH.set_attr("valid_min", 0.0); CHKERRQ(ierr);
    ierr = vars.add(vH); CHKERRQ(ierr);

    // bedrock surface elevation
    ierr = vbed.create(grid, "topg", true, WIDE_STENCIL); CHKERRQ(ierr);
    ierr = vbed.set_attrs("model_state", "bedrock surface elevation",
          "m", "bedrock_altitude"); CHKERRQ(ierr);
    ierr = vars.add(vbed); CHKERRQ(ierr);

    // yield stress for basal till (plastic or pseudo-plastic model)
    ierr = vtauc.create(grid, "tauc", false); CHKERRQ(ierr);
    ierr = vtauc.set_attrs("diagnostic", 
          "yield stress for basal till (plastic or pseudo-plastic model)",
          "Pa", ""); CHKERRQ(ierr);
    ierr = vars.add(vtauc); CHKERRQ(ierr);

    ierr = enthalpy.create(grid, "enthalpy", true); CHKERRQ(ierr);
    ierr = enthalpy.set_attrs("model_state",
                              "ice enthalpy (includes sensible heat, latent heat, pressure)",
                              "J kg-1", ""); CHKERRQ(ierr);
    ierr = vars.add(enthalpy); CHKERRQ(ierr);

    // grounded_dragging_floating integer mask
    ierr = vMask.create(grid, "mask", true, WIDE_STENCIL); CHKERRQ(ierr);
    ierr = vMask.set_attrs("model_state", "grounded_dragging_floating integer mask",
			 "", ""); CHKERRQ(ierr);
    vector<double> mask_values(6);
    mask_values[0] = MASK_ICE_FREE_BEDROCK;
    mask_values[1] = MASK_SHEET;
    mask_values[2] = MASK_DRAGGING_SHEET;
    mask_values[3] = MASK_FLOATING;
    mask_values[4] = MASK_ICE_FREE_OCEAN;
    mask_values[5] = MASK_OCEAN_AT_TIME_0;
    ierr = vMask.set_attr("flag_values", mask_values); CHKERRQ(ierr);
    ierr = vMask.set_attr("flag_meanings",
			"ice_free_bedrock sheet dragging_sheet floating ice_free_ocean ocean_at_time_zero");
		  CHKERRQ(ierr);
    vMask.output_data_type = NC_BYTE;
    ierr = vars.add(vMask); CHKERRQ(ierr);
  
    // read fields
    NCTool nc(grid.com, grid.rank);
    ierr = nc.open_for_reading(filename); CHKERRQ(ierr);
    int last_record;
    ierr = nc.get_dim_length("t", &last_record); CHKERRQ(ierr);
    ierr = nc.close(); CHKERRQ(ierr);

    last_record -= 1;
    ierr = vh.read(filename, last_record); CHKERRQ(ierr);
    ierr = vH.read(filename, last_record); CHKERRQ(ierr);
    ierr = vbed.read(filename, last_record); CHKERRQ(ierr);
    ierr = vMask.read(filename, last_record); CHKERRQ(ierr);
    ierr = vtauc.read(filename, last_record); CHKERRQ(ierr);
    ierr = enthalpy.read(filename, last_record); CHKERRQ(ierr);

    bool show = true;
    const PetscInt  window = 400;
    if (show) {
      ierr = vh.view(window);  CHKERRQ(ierr);
      ierr = vH.view(window);  CHKERRQ(ierr);
      ierr = vbed.view(window);  CHKERRQ(ierr);
      // the following breaks for an unknown reason:
      //      ierr = vtauc.view(window);  CHKERRQ(ierr);
      ierr = vMask.view(window);  CHKERRQ(ierr);
      PetscPrintf(grid.com,"  before SSA: showing fields in X windows for 5 seconds ...\n");
      ierr = PetscSleep(5); CHKERRQ(ierr);
    }

    SSAFD ssa(grid, basal, *ice, EC, config);

    ierr = ssa.init(vars); CHKERRQ(ierr); 
    
    ierr = ssa.update(false); CHKERRQ(ierr); 

    IceModelVec2V *vel_ssa;
    ierr = ssa.get_advective_2d_velocity(vel_ssa); CHKERRQ(ierr); 

    if (show) {
      ierr = vel_ssa->view(window);  CHKERRQ(ierr);
      PetscPrintf(grid.com,"[after SSA: showing components of velocity solution in X windows for 5 seconds ...]\n");
      ierr = PetscSleep(5); CHKERRQ(ierr);
    }

    ierr = vel_ssa->dump("vel_ssa_new.nc"); CHKERRQ(ierr); 

    if (ice != NULL)  delete ice;
  }
  ierr = PetscFinalize(); CHKERRQ(ierr);
  return 0;
}
