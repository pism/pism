// Copyright (C) 2010--2011 Ed Bueler and Constantine Khroulev
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

static char help[] =
  "\nSSAFD_TEST\n"
  "  Testing program for the finite difference implementation of the SSA.\n"
  "  Does a time-independent calculation.  Does not run IceModel or a derived\n"
  "  class thereof. Uses verification test J. Also may be used in a PISM\n"
  "  software (regression) test.\n\n";

#include "pism_const.hh"
#include "iceModelVec.hh"
#include "flowlaws.hh" // IceFlowLaw
#include "materials.hh" // IceBasalResistancePlasticLaw
#include "PISMIO.hh"
#include "NCVariable.hh"
#include "SSAFD.hh"
#include "exactTestsIJ.h"

//! \brief Report errors relative to the test J exact solution.
PetscErrorCode reportErrors(IceGrid &grid, NCConfigVariable &config,
                            IceModelVec2V &vel_ssa) {
  PetscErrorCode  ierr;
  PetscScalar exactmaxu, maxvecerr = 0.0, avvecerr = 0.0, 
    avuerr = 0.0, avverr = 0.0, maxuerr = 0.0, maxverr = 0.0;
  PetscScalar gmaxvecerr = 0.0, gavvecerr = 0.0, gavuerr = 0.0, gavverr = 0.0,
    gmaxuerr = 0.0, gmaxverr = 0.0;

  if (config.get_flag("do_pseudo_plastic_till")) {
    ierr = verbPrintf(1,grid.com, 
                      "WARNING: numerical errors not valid for pseudo-plastic till\n"); CHKERRQ(ierr);
  }
  ierr = verbPrintf(1,grid.com, 
                    "NUMERICAL ERRORS in velocity relative to exact solution:\n"); CHKERRQ(ierr);

  ierr = vel_ssa.begin_access(); CHKERRQ(ierr);

  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; i++) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; j++) {
      PetscScalar junk1, junk2, uexact, vexact;
      PetscScalar myr,myx,myy;
      grid.mapcoords(i,j,myx,myy,myr);
      // eval exact solution
      const PetscInt ifrom0 = i - (grid.Mx)/2,
        jfrom0 = j - (grid.My)/2;
      myx = grid.dx * ifrom0;
      myy = grid.dy * jfrom0;
      exactJ(myx, myy, &junk1, &junk2, &uexact, &vexact);

      // compute maximum errors
      const PetscScalar uerr = PetscAbsReal(vel_ssa(i,j).u - uexact);
      const PetscScalar verr = PetscAbsReal(vel_ssa(i,j).v - vexact);
      avuerr = avuerr + uerr;      
      avverr = avverr + verr;      
      maxuerr = PetscMax(maxuerr,uerr);
      maxverr = PetscMax(maxverr,verr);
      const PetscScalar vecerr = sqrt(uerr * uerr + verr * verr);
      maxvecerr = PetscMax(maxvecerr,vecerr);
      avvecerr = avvecerr + vecerr;
    }
  }
  ierr = vel_ssa.end_access(); CHKERRQ(ierr);
   
  ierr = PetscGlobalMax(&maxuerr, &gmaxuerr, grid.com); CHKERRQ(ierr);
  ierr = PetscGlobalMax(&maxverr, &gmaxverr, grid.com); CHKERRQ(ierr);
  ierr = PetscGlobalSum(&avuerr, &gavuerr, grid.com); CHKERRQ(ierr);
  gavuerr = gavuerr/(grid.Mx*grid.My);
  ierr = PetscGlobalSum(&avverr, &gavverr, grid.com); CHKERRQ(ierr);
  gavverr = gavverr/(grid.Mx*grid.My);
  ierr = PetscGlobalMax(&maxvecerr, &gmaxvecerr, grid.com); CHKERRQ(ierr);
  ierr = PetscGlobalSum(&avvecerr, &gavvecerr, grid.com); CHKERRQ(ierr);
  gavvecerr = gavvecerr/(grid.Mx*grid.My);

  // following from "pismv -test J -Mx 601 -My 601 -Mz 3 -verbose -eo"
  exactmaxu = 181.366 / secpera;

  ierr = verbPrintf(1,grid.com, 
                    "velocity  :  maxvector   prcntavvec      maxu      maxv       avu       avv\n");
  CHKERRQ(ierr);
  ierr = verbPrintf(1,grid.com, 
                    "           %11.4f%13.5f%10.4f%10.4f%10.4f%10.4f\n", 
                    gmaxvecerr*secpera, (gavvecerr/exactmaxu)*100.0,
                    gmaxuerr*secpera, gmaxverr*secpera, gavuerr*secpera, 
                    gavverr*secpera); CHKERRQ(ierr);

  ierr = verbPrintf(1,grid.com, "NUM ERRORS DONE\n");  CHKERRQ(ierr);

  return 0;
}

//! \brief Set the test J initial state.
PetscErrorCode setInitStateJ(NCConfigVariable &config,
                             IceGrid &grid,
                             SSAFD &ssa,
                             IceFlowLaw &ice,
                             IceModelVec2S *bed,
                             IceModelVec2Mask *mask,
                             IceModelVec2S *surface,
                             IceModelVec2S *thickness,
                             IceModelVec2V *vel_bc,
                             IceModelVec2S *tauc,
                             IceModelVec3 *enthalpy) {
  PetscErrorCode ierr;


  ierr = tauc->set(0.0); CHKERRQ(ierr);    // irrelevant for test J
  ierr = bed->set(-5000.0); CHKERRQ(ierr); // assures shelf is floating
  ierr = mask->set(MASK_FLOATING); CHKERRQ(ierr);

  ierr = enthalpy->set(528668.35); CHKERRQ(ierr); // arbitrary; corresponds to
                                                  // 263.15 Kelvin at depth=0.

  double ocean_rho = config.get("sea_water_density");

  /* use Ritz et al (2001) value of 30 MPa yr for typical vertically-averaged viscosity */
  const PetscScalar nu0 = 30.0 * 1.0e6 * secpera; /* = 9.45e14 Pa s */
  const PetscScalar H0 = 500.0;       /* 500 m typical thickness */

  // Use the same nuH factor everywhere (maximum ice thickness is 770 m).
  ssa.strength_extension.set_notional_strength(nu0 * H0);
  ssa.strength_extension.set_min_thickness(800);

  ierr = thickness->begin_access(); CHKERRQ(ierr);
  ierr = surface->begin_access(); CHKERRQ(ierr);
  ierr = mask->begin_access(); CHKERRQ(ierr);
  ierr = vel_bc->begin_access(); CHKERRQ(ierr);

  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      PetscScalar junk1, myu, myv, H;
      const PetscInt ifrom0 = i - (grid.Mx)/2,
                     jfrom0 = j - (grid.My)/2;
      const PetscScalar myx = grid.dx * ifrom0,
        myy = grid.dy * jfrom0;

      // set H,h on regular grid
      ierr = exactJ(myx, myy, &H, &junk1, &myu, &myv); CHKERRQ(ierr);

      (*thickness)(i,j) = H;
      (*surface)(i,j) = (1.0 - ice.rho / ocean_rho) * H; // FIXME task #7297

      // special case at center point: here we set vel_bc at (i,j) by marking
      // this grid point as SHEET and setting vel_bc approriately
      if ( (i == (grid.Mx)/2) && (j == (grid.My)/2) ) {
        (*mask)(i,j) = MASK_SHEET; // FIXME: replace with MASK_BC
        (*vel_bc)(i,j).u = myu;
        (*vel_bc)(i,j).v = myv;
      }
    }
  }  

  ierr = surface->end_access(); CHKERRQ(ierr);    
  ierr = thickness->end_access(); CHKERRQ(ierr);
  ierr = mask->end_access(); CHKERRQ(ierr);
  ierr = vel_bc->end_access(); CHKERRQ(ierr);

  // communicate what we have set
  ierr = surface->beginGhostComm(); CHKERRQ(ierr);
  ierr = surface->endGhostComm(); CHKERRQ(ierr);

  ierr = thickness->beginGhostComm(); CHKERRQ(ierr);
  ierr = thickness->endGhostComm(); CHKERRQ(ierr);

  ierr = mask->beginGhostComm(); CHKERRQ(ierr);
  ierr = mask->endGhostComm(); CHKERRQ(ierr);

  return 0;
}

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

    PetscTruth usage_set, help_set;
    ierr = PetscOptionsHasName(PETSC_NULL, "-usage", &usage_set); CHKERRQ(ierr);
    ierr = PetscOptionsHasName(PETSC_NULL, "-help", &help_set); CHKERRQ(ierr);
    if ((usage_set==PETSC_TRUE) || (help_set==PETSC_TRUE)) {
      PetscPrintf(com,
                  "\n"
                  "usage of SSA_TEST:\n"
                  "  run ssafd_test -Mx <number> -My <number>\n"
                  "\n");
    }

    IceGrid grid(com, rank, size, config);

    PetscReal LforJ = 300.0e3; // 300 km half-width

    grid.Lx = grid.Ly = LforJ;
    // grid is truly periodic both in x and in y directions
    grid.periodicity = XY_PERIODIC;
    grid.start_year = grid.year = 0.0;
    grid.Mx = grid.My = 61;
    grid.Mz = 3;
    
    string output_file = "ssafd_test_J.nc";
    ierr = PetscOptionsBegin(grid.com, "", "SSAFD_TEST options", ""); CHKERRQ(ierr);
    {
      bool flag;
      ierr = PISMOptionsInt("-Mx", "Number of grid points in the X direction",
                            grid.Mx, flag); CHKERRQ(ierr);
      ierr = PISMOptionsInt("-My", "Number of grid points in the X direction",
                            grid.My, flag); CHKERRQ(ierr);
      ierr = PISMOptionsString("-o", "Set the output file name",
                               output_file, flag); CHKERRQ(ierr);
    }
    ierr = PetscOptionsEnd(); CHKERRQ(ierr);

    grid.compute_nprocs();
    grid.compute_ownership_ranges();
    ierr = grid.compute_vertical_levels(); CHKERRQ(ierr); 
    ierr = grid.compute_horizontal_spacing(); CHKERRQ(ierr);
    ierr = grid.createDA(); CHKERRQ(ierr);

    ierr = setVerbosityLevel(5); CHKERRQ(ierr);
    ierr = grid.printInfo(1); CHKERRQ(ierr);
    //ierr = grid.printVertLevels(1); CHKERRQ(ierr); 

    CustomGlenIce ice(grid.com, "", config);

    IceBasalResistancePlasticLaw basal(
           config.get("plastic_regularization") / secpera, 
           config.get_flag("do_pseudo_plastic_till"),
           config.get("pseudo_plastic_q"),
           config.get("pseudo_plastic_uthreshold") / secpera);
    ierr = basal.printInfo(1,grid.com); CHKERRQ(ierr);

    EnthalpyConverter EC(config);

    IceModelVec2S vh, vH, vbed, vtauc;
    IceModelVec2Mask vMask;
    IceModelVec3 enthalpy;
    IceModelVec2V vel_bc;
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

    ierr = vel_bc.create(grid, "_bc", false); CHKERRQ(ierr); // u_bc and v_bc
    ierr = vel_bc.set_attrs("intent",
                            "X-component of the SSA velocity boundary conditions",
                            "m s-1", "", 0); CHKERRQ(ierr);
    ierr = vel_bc.set_attrs("intent",
                            "Y-component of the SSA velocity boundary conditions",
                            "m s-1", "", 1); CHKERRQ(ierr);
    ierr = vel_bc.set_glaciological_units("m year-1"); CHKERRQ(ierr);
    ierr = vel_bc.set_attr("valid_min", -1e6/secpera, 0); CHKERRQ(ierr); 
    ierr = vel_bc.set_attr("valid_max",  1e6/secpera, 0); CHKERRQ(ierr); 
    ierr = vel_bc.set_attr("valid_min", -1e6/secpera, 1); CHKERRQ(ierr); 
    ierr = vel_bc.set_attr("valid_max",  1e6/secpera, 1); CHKERRQ(ierr); 
    ierr = vel_bc.set_attr("_FillValue", 2e6/secpera, 0); CHKERRQ(ierr); 
    ierr = vel_bc.set_attr("_FillValue", 2e6/secpera, 1); CHKERRQ(ierr); 
    vel_bc.write_in_glaciological_units = true;
    ierr = vel_bc.set(2e6/secpera); CHKERRQ(ierr); 

    // Create the SSA solver object:
    SSAFD ssa(grid, basal, ice, EC, config);

    // Allocate the SSA solver:
    ierr = ssa.init(vars); CHKERRQ(ierr);

    // fill the fields:
    ierr = setInitStateJ(config, grid, ssa, ice,
                         &vbed, &vMask, &vh, &vH,
                         &vel_bc, &vtauc, &enthalpy); CHKERRQ(ierr);

    ierr = ssa.set_boundary_conditions(vMask, vel_bc); CHKERRQ(ierr); 

    // Solve (fast==true means "no update"):
    bool fast = false;
    ierr = ssa.update(fast); CHKERRQ(ierr); 

    // Report errors relative to the exact solution:
    IceModelVec2V *vel_ssa;
    ierr = ssa.get_advective_2d_velocity(vel_ssa); CHKERRQ(ierr); 

    ierr = reportErrors(grid, config, *vel_ssa); CHKERRQ(ierr);

    // Write results to an output file:
    PISMIO pio(&grid);

    ierr = pio.open_for_writing(output_file, false, true); CHKERRQ(ierr);
    ierr = pio.append_time(0.0);
    ierr = pio.close(); CHKERRQ(ierr); 

    ierr = vh.write(output_file.c_str()); CHKERRQ(ierr);
    ierr = vH.write(output_file.c_str()); CHKERRQ(ierr);
    ierr = vMask.write(output_file.c_str()); CHKERRQ(ierr);
    ierr = vtauc.write(output_file.c_str()); CHKERRQ(ierr);
    ierr = vbed.write(output_file.c_str()); CHKERRQ(ierr);
    ierr = enthalpy.write(output_file.c_str()); CHKERRQ(ierr);
    ierr = vel_bc.write(output_file.c_str()); CHKERRQ(ierr);
    ierr = vel_ssa->write(output_file.c_str()); CHKERRQ(ierr);
  }
  ierr = PetscFinalize(); CHKERRQ(ierr);
  return 0;
}
