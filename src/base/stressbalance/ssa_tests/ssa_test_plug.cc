// Copyright (C) 2010--2011 Ed Bueler, Constantine Khroulev, and David Maxwell
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

/* This file implements a test case for the ssa:  plug flow.  The geometry consists of a constant
   surface slope in the positive x-direction, and the ice is pinned on the y-boundaries.  There is
   no basal shear stress, and hence the the only nonzero terms in the SSA are the "p-laplacian" 
   and the driving stress. */

static char help[] =
  "\nSSA_TEST_PLUG\n"
  "  Testing program for the finite element implementation of the SSA.\n"
  "  Does a time-independent calculation.  Does not run IceModel or a derived\n"
  "  class thereof.\n\n";

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
#include <math.h>
class SSATestCasePlug: public SSATestCase
{
public:
  SSATestCasePlug( MPI_Comm com, PetscMPIInt rank, 
                 PetscMPIInt size, NCConfigVariable &config ): 
                 SSATestCase(com,rank,size,config)
  { 
    H0 = 2000.; //m
    L=50.e3; // 50km half-width
    dhdx = 0.001; // pure number, slope of surface & bed
    tauc0 = 0.; // No basal shear stress
    B0 = 3.7e8; // Pa s^{1/3}; hardness 
               // given on p. 239 of Schoof; why so big?      
  }
  
protected:
  virtual PetscErrorCode initializeGrid(PetscInt Mx,PetscInt My);

  virtual PetscErrorCode initializeSSAModel();

  virtual PetscErrorCode initializeSSACoefficients();

  virtual PetscErrorCode exactSolution(PetscInt i, PetscInt j, 
    PetscReal x, PetscReal y, PetscReal *u, PetscReal *v );


  PetscScalar H0; // Thickness
  PetscScalar L;  // Half-width
  PetscScalar dhdx; // surface slope
  PetscScalar tauc0; // zero basal shear stress
  PetscScalar B0;  // hardness

  bool dimensionless;
  
};


PetscErrorCode SSATestCasePlug::initializeGrid(PetscInt Mx,PetscInt My)
{
  PetscReal Lx=L, Ly = L; 
  init_shallow_grid(grid,Lx,Ly,Mx,My,NONE);
}


PetscErrorCode SSATestCasePlug::initializeSSAModel()
{
  // The following is irrelevant because tauc=0
  PetscScalar linear_q = 1.;
  basal = new IceBasalResistancePlasticLaw(
         config.get("plastic_regularization") / secpera,
         true, // do not force a pure-plastic law
         linear_q,
         config.get("pseudo_plastic_uthreshold") / secpera);

  // Use constant hardness
  CustomGlenIce *glenIce = new CustomGlenIce(grid.com, "", config);
  glenIce->setHardness(B0);
  ice=glenIce;

  // Enthalpy converter is irrelevant for this test.
  enthalpyconverter = new EnthalpyConverter(config);
  return 0;
}

PetscErrorCode SSATestCasePlug::initializeSSACoefficients()
{
  PetscErrorCode ierr;
  PetscScalar    **ph, **pbed;
  
  // The finite difference code uses the following flag to treat the non-periodic grid correctly.
  config.set_flag("compute_surf_grad_inward_ssa", true);
  config.set("epsilon_ssa", 0.);

  // Ensure we never use the strength extension.
  ssa->strength_extension->set_min_thickness(H0/2);

  // Set constant coefficients.
  ierr = thickness.set(H0); CHKERRQ(ierr);
  ierr = tauc.set(tauc0); CHKERRQ(ierr);
  

  // Set boundary conditions (Dirichlet all the way around).
  ierr = mask.set(MASK_DRAGGING_SHEET); CHKERRQ(ierr);
  ierr = vel_bc.begin_access(); CHKERRQ(ierr);
  ierr = mask.begin_access(); CHKERRQ(ierr);
  ierr = bed.begin_access(); CHKERRQ(ierr);
  ierr = surface.begin_access(); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      PetscScalar junk, myu, myv;
      const PetscScalar myx = grid.x[i], myy=grid.y[j];

      bed(i,j) = -myx*(dhdx);
      surface(i,j) = bed(i,j) + H0;
      
      bool edge = ( (j == 0) || (j == grid.My - 1) ) || ( (i==0) || (i==grid.Mx-1) );
      if (edge) {
        mask(i,j) = MASK_BC;
        exactSolution(i,j,myx,myy,&myu,&myv);
        vel_bc(i,j).u = myu;
        vel_bc(i,j).v = myv;
      }
    }
  } 
  ierr = vel_bc.end_access(); CHKERRQ(ierr);
  ierr = mask.end_access(); CHKERRQ(ierr);
  ierr = bed.end_access(); CHKERRQ(ierr);
  ierr = surface.end_access(); CHKERRQ(ierr);
  
  
  ierr = vel_bc.beginGhostComm(); CHKERRQ(ierr);
  ierr = vel_bc.endGhostComm(); CHKERRQ(ierr);
  ierr = mask.beginGhostComm(); CHKERRQ(ierr);
  ierr = mask.endGhostComm(); CHKERRQ(ierr);
  ierr = bed.beginGhostComm(); CHKERRQ(ierr);
  ierr = bed.endGhostComm(); CHKERRQ(ierr);  
  ierr = surface.beginGhostComm(); CHKERRQ(ierr);
  ierr = surface.endGhostComm(); CHKERRQ(ierr);  

  ierr = ssa->set_boundary_conditions(mask, vel_bc); CHKERRQ(ierr); 

  return 0;
}

PetscErrorCode SSATestCasePlug::exactSolution(PetscInt i, PetscInt j, 
 PetscReal x, PetscReal y, PetscReal *u, PetscReal *v)
{
  PetscScalar earth_grav = config.get("standard_gravity");
  PetscScalar f = ice->rho * earth_grav * H0* dhdx;
  PetscScalar ynd = y/L;
  
  *u = 0.5*pow(f,3)*pow(L,4)/pow(B0*H0,3)*(1-pow(ynd,4));
  *v = 0;
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
    PetscInt Mx=11;
    PetscInt My=61;
    string output_file = "ssa_test_plug.nc";
    string driver = "fem";
    PetscScalar H0dim = 1.;

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

    SSATestCasePlug testcase(com,rank,size,config);
    ierr = testcase.init(Mx,My,ssafactory); CHKERRQ(ierr);
    ierr = testcase.run(); CHKERRQ(ierr);
    ierr = testcase.report(); CHKERRQ(ierr);
    ierr = testcase.write(output_file); CHKERRQ(ierr);
  }

  ierr = PetscFinalize(); CHKERRQ(ierr);
  return 0;
}
