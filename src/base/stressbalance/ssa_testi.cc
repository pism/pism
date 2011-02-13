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

class SSATestCaseI: public SSATestCase
{
public:
  SSATestCaseI( MPI_Comm com, PetscMPIInt rank, 
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

const PetscScalar m_schoof = 10; // (pure number)
const PetscScalar L_schoof = 40e3; // meters
const PetscScalar aspect_schoof = 0.05; // (pure)
const PetscScalar H0_schoof = aspect_schoof * L_schoof; 
                                       // = 2000 m THICKNESS
const PetscScalar B_schoof = 3.7e8; // Pa s^{1/3}; hardness 
                                     // given on p. 239 of Schoof; why so big?
const PetscScalar p_schoof = 4.0/3.0; // = 1 + 1/n


PetscErrorCode SSATestCaseI::initializeGrid(PetscInt Mx,PetscInt My)
{
  PetscReal Ly = 3*L_schoof;  // 300.0 km half-width (L=40.0km in Schoof's choice of variables)
  PetscReal Lx = PetscMax(60.0e3, ((Mx - 1) / 2) * (2.0 * Ly / (My - 1)) );


  PetscErrorCode ierr;
  
  grid.Lx = Lx;
  grid.Ly = Ly;
  grid.periodicity = NONE;
  grid.start_year = grid.year = 0.0;
  grid.Mx = Mx; grid.My=My; grid.Mz = 3;
  
  grid.compute_nprocs();
  grid.compute_ownership_ranges();
  ierr = grid.compute_vertical_levels(); CHKERRQ(ierr); 
  ierr = grid.compute_horizontal_spacing(); CHKERRQ(ierr);
  ierr = grid.createDA(); CHKERRQ(ierr);

  return 0;


  // // Fixme: make only y-periodic
  // init_shallow_periodic_grid(grid,Lx,Ly,Mx,My);
}


PetscErrorCode SSATestCaseI::initializeSSAModel()
{
  PetscReal nearPlasticQ = 0.05;
  basal = new IceBasalResistancePlasticLaw(
         config.get("plastic_regularization") / secpera,
         config.get_flag("do_pseudo_plastic_till"),
         nearPlasticQ,
         config.get("pseudo_plastic_uthreshold") / secpera);

  ice = new CustomGlenIce(grid.com, "", config);
  enthalpyconverter = new EnthalpyConverter(config);
  return 0;
}

PetscErrorCode SSATestCaseI::initializeSSACoefficients()
{
  PetscErrorCode ierr;
  PetscScalar    **ph, **pbed;
  
  ierr = mask.set(MASK_DRAGGING_SHEET); CHKERRQ(ierr);
  ierr = thickness.set(H0_schoof); CHKERRQ(ierr);

  // set h, bed everywhere
  // on edges y = +- 3 L_schoof, set velocity and make mask=SHEET
  // ierr = vMask.begin_access(); CHKERRQ(ierr);
  // Test J has a viscosity that is independent of velocity.  So we force a 
  // constant viscosity by settting the strength_extension
  // thickness larger than the given ice thickness. (max = 770m).
  const PetscScalar nu0 = 30.0 * 1.0e6 * secpera; /* = 9.45e14 Pa s */
  const PetscScalar H0 = 500.0;       /* 500 m typical thickness */
  ssa->strength_extension->set_notional_strength(nu0 * H0);
  ssa->strength_extension->set_min_thickness(4000);

  // The finite difference code uses the following flag to treat the non-periodic grid correctly.
  config.set_flag("compute_surf_grad_inward_ssa", true);
  config.set("epsilon_ssa", 0.0);  // don't use this lower bound

  // The finite difference version 

  PetscScalar **ptauc;

  ierr = tauc.get_array(ptauc); CHKERRQ(ierr);

  PetscScalar standard_gravity = config.get("standard_gravity");

  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      const PetscInt jfrom0 = j - (grid.My - 1)/2;
      const PetscScalar y = grid.dy * jfrom0;
      const PetscScalar theta = atan(0.001);   /* a slope of 1/1000, a la Siple streams */
      const PetscScalar f = ice->rho * standard_gravity * H0_schoof * tan(theta);
      ptauc[i][j] = f * pow(PetscAbs(y / L_schoof), m_schoof);
    }
  }
  ierr = tauc.end_access(); CHKERRQ(ierr);
  ierr = tauc.beginGhostComm(); CHKERRQ(ierr);
  ierr = tauc.endGhostComm(); CHKERRQ(ierr);


  ierr = surface.get_array(ph); CHKERRQ(ierr);    
  ierr = bed.get_array(pbed); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      PetscScalar junk, myu, myv, bedval;
      const PetscInt ifrom0 = i - (grid.Mx - 1)/2,
                     jfrom0 = j - (grid.My - 1)/2;
      const PetscScalar myx = grid.dx * ifrom0,
        myy = grid.dy * jfrom0;
      // eval exact solution; will only use exact vels if at edge
      exactI(m_schoof, myx, myy, &(pbed[i][j]), &junk, &myu, &myv); 
      ph[i][j] = pbed[i][j] + H0_schoof;
    }
  }  
  ierr = vel_bc.end_access(); CHKERRQ(ierr);
  ierr = surface.end_access(); CHKERRQ(ierr);    
  ierr = bed.end_access(); CHKERRQ(ierr);

  // communicate what we have set
  ierr = surface.beginGhostComm(); CHKERRQ(ierr);
  ierr = surface.endGhostComm(); CHKERRQ(ierr);
  ierr = bed.beginGhostComm(); CHKERRQ(ierr);
  ierr = bed.endGhostComm(); CHKERRQ(ierr);

  // ierr = ssa->set_boundary_conditions(mask, &NONE); CHKERRQ(ierr); 

  return 0;
}


PetscErrorCode SSATestCaseI::exactSolution(PetscInt i, PetscInt j, 
 PetscReal x, PetscReal y, PetscReal *u, PetscReal *v)
{
  PetscReal junk1, junk2;
  const PetscInt ifrom0 = i - (grid.Mx)/2, jfrom0 = j - (grid.My)/2;
  PetscReal myx = grid.dx * ifrom0;
  PetscReal myy = grid.dy * jfrom0;

  exactI(m_schoof, myx,myy, &junk1, &junk2,u,v); 
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
    string output_file = "ssa_testi.nc";
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

    SSATestCaseI testcase(com,rank,size,config);
    ierr = testcase.init(Mx,My,ssafactory); CHKERRQ(ierr);
    ierr = testcase.run(); CHKERRQ(ierr);
    ierr = testcase.report(); CHKERRQ(ierr);
    ierr = testcase.write(output_file);
  }

  ierr = PetscFinalize(); CHKERRQ(ierr);
  return 0;
}
