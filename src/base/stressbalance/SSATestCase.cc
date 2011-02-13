// Copyright (C) 2009--2011 Ed Bueler, Constantine Khroulev, and David Maxwell
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

#include "SSATestCase.hh"
#include "PISMIO.hh"

#include "SSAFD.hh"
#include "SSAFEM.hh"

SSA *SSAFEMFactory(IceGrid &g, IceBasalResistancePlasticLaw &b, 
                IceFlowLaw &i, EnthalpyConverter &ec, 
                const NCConfigVariable &c)
{
  return new SSAFEM(g,b,i,ec,c);
}

SSA *SSAFDFactory(IceGrid &g, IceBasalResistancePlasticLaw &b, 
                IceFlowLaw &i, EnthalpyConverter &ec, 
                const NCConfigVariable &c)
{
  return new SSAFD(g,b,i,ec,c);
}

//! Initialize the storage for the various coefficients used as input to the SSA
//! (ice elevation, thickness, etc.)  
PetscErrorCode SSATestCase::buildSSACoefficients()
{
  PetscErrorCode ierr;

  const PetscInt WIDE_STENCIL = 2;
  
  // ice surface elevation
  ierr = surface.create(grid, "usurf", true, WIDE_STENCIL); CHKERRQ(ierr);
  ierr = surface.set_attrs("diagnostic", "ice upper surface elevation", "m", 
                                      "surface_altitude"); CHKERRQ(ierr);
  ierr = vars.add(surface); CHKERRQ(ierr);
  
  // land ice thickness
  ierr = thickness.create(grid, "thk", true, WIDE_STENCIL); CHKERRQ(ierr);
  ierr = thickness.set_attrs("model_state", "land ice thickness", "m", 
                                  "land_ice_thickness"); CHKERRQ(ierr);
  ierr = thickness.set_attr("valid_min", 0.0); CHKERRQ(ierr);
  ierr = vars.add(thickness); CHKERRQ(ierr);

  // bedrock surface elevation
  ierr = bed.create(grid, "topg", true, WIDE_STENCIL); CHKERRQ(ierr);
  ierr = bed.set_attrs("model_state", "bedrock surface elevation", "m", 
                                          "bedrock_altitude"); CHKERRQ(ierr);
  ierr = vars.add(bed); CHKERRQ(ierr);

  // yield stress for basal till (plastic or pseudo-plastic model)
  ierr = tauc.create(grid, "tauc", true, WIDE_STENCIL); CHKERRQ(ierr);
  ierr = tauc.set_attrs("diagnostic",  
  "yield stress for basal till (plastic or pseudo-plastic model)", "Pa", ""); 
      CHKERRQ(ierr);
  ierr = vars.add(tauc); CHKERRQ(ierr);

  // enthalpy
  ierr = enthalpy.create(grid, "enthalpy", true, WIDE_STENCIL); CHKERRQ(ierr);
  ierr = enthalpy.set_attrs("model_state",
              "ice enthalpy (includes sensible heat, latent heat, pressure)",
              "J kg-1", ""); CHKERRQ(ierr);
  ierr = vars.add(enthalpy); CHKERRQ(ierr);


  // dirichlet boundary condition (FIXME: perhaps unused!)
  ierr = vel_bc.create(grid, "_bc", true,WIDE_STENCIL); CHKERRQ(ierr); // u_bc and v_bc
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
  
  // grounded_dragging_floating integer mask
  ierr = mask.create(grid, "mask", true, WIDE_STENCIL); CHKERRQ(ierr);
  ierr = mask.set_attrs("model_state", 
          "grounded_dragging_floating integer mask", "", ""); CHKERRQ(ierr);
  vector<double> mask_values(6);
  mask_values[0] = MASK_ICE_FREE_BEDROCK;
  mask_values[1] = MASK_SHEET;
  mask_values[2] = MASK_DRAGGING_SHEET;
  mask_values[3] = MASK_FLOATING;
  mask_values[4] = MASK_ICE_FREE_OCEAN;
  mask_values[5] = MASK_OCEAN_AT_TIME_0;
  ierr = mask.set_attr("flag_values", mask_values); CHKERRQ(ierr);
  ierr = mask.set_attr("flag_meanings",
    "ice_free_bedrock sheet dragging_sheet floating ice_free_ocean ocean_at_time_zero");
             CHKERRQ(ierr);
  mask.output_data_type = NC_BYTE;
  ierr = vars.add(mask); CHKERRQ(ierr);  
  
  return 0;
}


//! Initialize the test case at the start of a run
PetscErrorCode SSATestCase::init(PetscInt Mx, PetscInt My, SSAFactory ssafactory)
{
  PetscErrorCode ierr;
  
  // Subclass builds grid.
  ierr = initializeGrid(Mx,My);
  
  // Subclass builds ice flow law, basal resistance, etc.
  ierr = initializeSSAModel(); CHKERRQ(ierr);

  // We setup storage for the coefficients.
  ierr = buildSSACoefficients(); CHKERRQ(ierr);

  // Allocate the actual SSA solver.
  ssa = ssafactory(grid, *basal, *ice, *enthalpyconverter, config);
  ierr = ssa->init(vars); CHKERRQ(ierr); // vars was setup preivouisly with buildSSACoefficients

  // Allow the subclass to setup the coefficients.
  ierr = initializeSSACoefficients(); CHKERRQ(ierr);

  return 0;
}

//! Solve the SSA
PetscErrorCode SSATestCase::run()
{
  PetscErrorCode ierr;
  // Solve (fast==true means "no update"):
  ierr = verbPrintf(2,grid.com,"* Solving the SSA stress balance ...\n"); CHKERRQ(ierr);

  bool fast = false;
  ierr = ssa->update(fast); CHKERRQ(ierr); 

  return 0;
}

//! Report on the generated solution
PetscErrorCode SSATestCase::report()
{
  PetscErrorCode  ierr;
    
  string ssa_stdout;
  ierr = ssa->stdout_report(ssa_stdout); CHKERRQ(ierr);
  ierr = verbPrintf(3,grid.com,ssa_stdout.c_str()); CHKERRQ(ierr);
  
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


  IceModelVec2V *vel_ssa;
  ierr = ssa->get_advective_2d_velocity(vel_ssa); CHKERRQ(ierr);
  ierr = vel_ssa->begin_access(); CHKERRQ(ierr);

  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; i++) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; j++) {
      PetscScalar junk1, junk2, uexact, vexact;
      PetscScalar myr,myx,myy;
      grid.mapcoords(i,j,myx,myy,myr);
      exactSolution(i,j,myx,myy,&uexact,&vexact);

      // compute maximum errors
      const PetscScalar uerr = PetscAbsReal((*vel_ssa)(i,j).u - uexact);
      const PetscScalar verr = PetscAbsReal((*vel_ssa)(i,j).v - vexact);
      avuerr = avuerr + uerr;      
      avverr = avverr + verr;      
      maxuerr = PetscMax(maxuerr,uerr);
      maxverr = PetscMax(maxverr,verr);
      const PetscScalar vecerr = sqrt(uerr * uerr + verr * verr);
      maxvecerr = PetscMax(maxvecerr,vecerr);
      avvecerr = avvecerr + vecerr;
    }
  }
  ierr = vel_ssa->end_access(); CHKERRQ(ierr);
   
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

PetscErrorCode SSATestCase::exactSolution(PetscInt i, PetscInt j, 
  PetscReal x, PetscReal y, PetscReal *u, PetscReal *v )
{
  *u=0; *v=0;
}

//! Save the computation and data to a file.
PetscErrorCode SSATestCase::write(const string &filename)
{
  PetscErrorCode ierr;

  // Write results to an output file:
  PISMIO pio(&grid);
  ierr = pio.open_for_writing(filename, false, true); CHKERRQ(ierr);
  ierr = pio.append_time(0.0); CHKERRQ(ierr);
  ierr = pio.close(); CHKERRQ(ierr); 

  ierr = surface.write(filename.c_str()); CHKERRQ(ierr);
  ierr = thickness.write(filename.c_str()); CHKERRQ(ierr);
  ierr = mask.write(filename.c_str()); CHKERRQ(ierr);
  ierr = tauc.write(filename.c_str()); CHKERRQ(ierr);
  ierr = bed.write(filename.c_str()); CHKERRQ(ierr);
  ierr = enthalpy.write(filename.c_str()); CHKERRQ(ierr);
  ierr = vel_bc.write(filename.c_str()); CHKERRQ(ierr);

  IceModelVec2V *vel_ssa;
  ierr = ssa->get_advective_2d_velocity(vel_ssa); CHKERRQ(ierr);
  ierr = vel_ssa->write(filename.c_str()); CHKERRQ(ierr);

  return 0;
}


/*! Initialize a uniform, shallow (3 z-levels), doubly periodic grid with 
half-widths (Lx,Ly) and Mx by My nodes for time-independent computations.*/
PetscErrorCode init_shallow_periodic_grid(IceGrid &grid, PetscReal Lx, 
                                      PetscReal Ly, PetscInt Mx, PetscInt My)
{
  PetscErrorCode ierr;
  
  grid.Lx = Lx;
  grid.Ly = Ly;
  grid.periodicity = XY_PERIODIC;
  grid.start_year = grid.year = 0.0;
  grid.Mx = Mx; grid.My=My; grid.Mz = 3;
  
  grid.compute_nprocs();
  grid.compute_ownership_ranges();
  ierr = grid.compute_vertical_levels(); CHKERRQ(ierr); 
  ierr = grid.compute_horizontal_spacing(); CHKERRQ(ierr);
  ierr = grid.createDA(); CHKERRQ(ierr);

  return 0;
}

