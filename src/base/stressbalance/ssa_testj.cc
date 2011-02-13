// Copyright (C) 2010--2011 Ed Bueler, Constantine Khroulov, and David Maxwell
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
  "  Testing program for the finite element implementation of the SSA.\n"
  "  Does a time-independent calculation.  Does not run IceModel or a derived\n"
  "  class thereof. Uses verification test J. Also may be used in a PISM\n"
  "  software (regression) test.\n\n";

#include "pism_const.hh"
#include "iceModelVec.hh"
#include "flowlaws.hh" // IceFlowLaw
#include "materials.hh" // IceBasalResistancePlasticLaw
#include "PISMIO.hh"
#include "NCVariable.hh"
#include "SSAFEM.hh"
#include "SSAFD.hh"
#include "exactTestsIJ.h"
#include "SSATestCase.hh"

class SSATestCaseJ: public SSATestCase
{
public:
  SSATestCaseJ( MPI_Comm com, PetscMPIInt rank, 
                 PetscMPIInt size, NCConfigVariable &config ): 
                 SSATestCase(com,rank,size,config)
  { };
  
protected:
  virtual PetscErrorCode initializeGrid(PetscInt Mx,PetscInt My);

  virtual PetscErrorCode initializeSSAModel();

  virtual PetscErrorCode initializeSSACoefficients();

  virtual PetscErrorCode exactSolution(PetscInt i, PetscInt j, 
    PetscReal x, PetscReal y, PetscReal *u, PetscReal *v );

};

PetscErrorCode SSATestCaseJ::initializeGrid(PetscInt Mx,PetscInt My)
{
  PetscReal halfWidth = 300.0e3;  // 300.0 km half-width
  PetscReal Lx = halfWidth, Ly = halfWidth;
  init_shallow_periodic_grid(grid,Lx,Ly,Mx,My);
}

PetscErrorCode SSATestCaseJ::initializeSSAModel()
{
  basal = new IceBasalResistancePlasticLaw(
         config.get("plastic_regularization") / secpera,
         config.get_flag("do_pseudo_plastic_till"),
         config.get("pseudo_plastic_q"),
         config.get("pseudo_plastic_uthreshold") / secpera);

  ice = new CustomGlenIce(grid.com, "", config);
  enthalpyconverter = new EnthalpyConverter(config);
  return 0;
}

PetscErrorCode SSATestCaseJ::initializeSSACoefficients()
{
  PetscErrorCode ierr;
  ierr = tauc.set(0.0); CHKERRQ(ierr);    // irrelevant for test J
  ierr = bed.set(0.0); CHKERRQ(ierr); // assures shelf is floating
  ierr = mask.set(MASK_FLOATING); CHKERRQ(ierr);
  ierr = enthalpy.set(528668.35); 
        CHKERRQ(ierr); // arbitrary; corresponds to 263.15 Kelvin at depth=0.

  /* use Ritz et al (2001) value of 30 MPa yr for typical vertically-averaged viscosity */
  double ocean_rho = config.get("sea_water_density");
  const PetscScalar nu0 = 30.0 * 1.0e6 * secpera; /* = 9.45e14 Pa s */
  const PetscScalar H0 = 500.0;       /* 500 m typical thickness */

  // Test J has a viscosity that is independent of velocity.  So we force a 
  // constant viscosity by settting the strength_extension
  // thickness larger than the given ice thickness. (max = 770m).
  ssa->strength_extension->set_notional_strength(nu0 * H0);
  ssa->strength_extension->set_min_thickness(800);

  ierr = thickness.begin_access(); CHKERRQ(ierr);
  ierr = surface.begin_access(); CHKERRQ(ierr);
  ierr = mask.begin_access(); CHKERRQ(ierr);
  ierr = vel_bc.begin_access(); CHKERRQ(ierr);

  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      PetscScalar junk1, myu, myv, H;
      const PetscInt ifrom0 = i - (grid.Mx)/2,
                     jfrom0 = j - (grid.My)/2;
      const PetscScalar myx = grid.dx * ifrom0,
                        myy = grid.dy * jfrom0;

      // set H,h on regular grid
      ierr = exactJ(myx, myy, &H, &junk1, &myu, &myv); CHKERRQ(ierr);

      thickness(i,j) = H;
      surface(i,j) = (1.0 - ice->rho / ocean_rho) * H; // FIXME task #7297

      // special case at center point: here we set vel_bc at (i,j) by marking
      // this grid point as SHEET and setting vel_bc approriately
      if ( (i == (grid.Mx)/2) && (j == (grid.My)/2) ) {
        mask(i,j) = MASK_SHEET; // FIXME: replace with MASK_BC
        vel_bc(i,j).u = myu;
        vel_bc(i,j).v = myv;
      }
    }
  }  

  ierr = surface.end_access(); CHKERRQ(ierr);    
  ierr = thickness.end_access(); CHKERRQ(ierr);
  ierr = mask.end_access(); CHKERRQ(ierr);
  ierr = vel_bc.end_access(); CHKERRQ(ierr);

  // communicate what we have set
  ierr = surface.beginGhostComm(); CHKERRQ(ierr);
  ierr = surface.endGhostComm(); CHKERRQ(ierr);

  ierr = thickness.beginGhostComm(); CHKERRQ(ierr);
  ierr = thickness.endGhostComm(); CHKERRQ(ierr);

  ierr = mask.beginGhostComm(); CHKERRQ(ierr);
  ierr = mask.endGhostComm(); CHKERRQ(ierr);

  ierr = vel_bc.beginGhostComm(); CHKERRQ(ierr);
  ierr = vel_bc.endGhostComm(); CHKERRQ(ierr);

  ierr = ssa->set_boundary_conditions(mask, vel_bc); 
      CHKERRQ(ierr); 

  return 0;    
}

PetscErrorCode SSATestCaseJ::exactSolution(PetscInt i, PetscInt j, 
 PetscReal x, PetscReal y, PetscReal *u, PetscReal *v)
{
  PetscReal junk1, junk2;
  const PetscInt ifrom0 = i - (grid.Mx)/2, jfrom0 = j - (grid.My)/2;
  PetscReal myx = grid.dx * ifrom0;
  PetscReal myy = grid.dy * jfrom0;
  exactJ(myx, myy, &junk1, &junk2, u, v);
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

    ierr = setVerbosityLevel(5); CHKERRQ(ierr);

    PetscTruth usage_set, help_set;
    ierr = PetscOptionsHasName(PETSC_NULL, "-usage", &usage_set); CHKERRQ(ierr);
    ierr = PetscOptionsHasName(PETSC_NULL, "-help", &help_set); CHKERRQ(ierr);
    if ((usage_set==PETSC_TRUE) || (help_set==PETSC_TRUE)) {
      PetscPrintf(com,
                  "\n"
                  "usage of SSA_TESTJ:\n"
                  "  run ssafe_test -Mx <number> -My <number> -ssa <fd|fem>\n"
                  "\n");
    }

    // Parameters that can be overridden by command line options
    PetscInt Mx=61;
    PetscInt My=61;
    string output_file = "ssa_testj.nc";
    string driver = "fem";

    ierr = PetscOptionsBegin(com, "", "SSAFD_TEST options", ""); CHKERRQ(ierr);
    {
      bool flag;
      ierr = PISMOptionsInt("-Mx", "Number of grid points in the X direction", 
                                                      Mx, flag); CHKERRQ(ierr);
      ierr = PISMOptionsInt("-My", "Number of grid points in the Y direction", 
                                                      My, flag); CHKERRQ(ierr);
      ierr = PISMOptionsString("-ssa_method", "Algorithm for computing the SSA solution",
                                                  driver, flag); CHKERRQ(ierr);                                                      
      ierr = PISMOptionsString("-o", "Set the output file name", 
                                              output_file, flag); CHKERRQ(ierr);
    }
    ierr = PetscOptionsEnd(); CHKERRQ(ierr);

    // Determine the kind of solver to use.
    SSAFactory ssafactory;
    if(driver.compare("fem") == 0) ssafactory = SSAFEMFactory;
    else if(driver.compare("fd") == 0) ssafactory = SSAFDFactory;
    else SETERRQ(1,"SSA algorithm argument should be one of -ssa fe or -ssa fem");

    SSATestCaseJ testcase(com,rank,size,config);
    ierr = testcase.init(Mx,My,ssafactory); CHKERRQ(ierr);
    ierr = testcase.run(); CHKERRQ(ierr);
    ierr = testcase.report(); CHKERRQ(ierr);
    ierr = testcase.write(output_file);
  }

  ierr = PetscFinalize(); CHKERRQ(ierr);
  return 0;
}
