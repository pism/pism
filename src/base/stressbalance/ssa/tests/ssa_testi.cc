// Copyright (C) 2010--2014 Ed Bueler, Constantine Khroulev, and David Maxwell
//
// This file is part of PISM.
//
// PISM is free software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the Free Software
// Foundation; either version 3 of the License, or (at your option) any later
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
  "\nSSA_TESTI\n"
  "  Testing program for the finite element implementation of the SSA.\n"
  "  Does a time-independent calculation.  Does not run IceModel or a derived\n"
  "  class thereof. Uses verification test I. Also may be used in a PISM\n"
  "  software (regression) test.\n\n";

#include "pism_const.hh"
#include "pism_options.hh"
#include "iceModelVec.hh"
#include "flowlaws.hh" // IceFlowLaw
#include "basal_resistance.hh" // IceBasalResistancePlasticLaw
#include "PIO.hh"
#include "NCVariable.hh"
#include "SSAFEM.hh"
#include "SSAFD.hh"
#include "exactTestsIJ.h"
#include "SSATestCase.hh"

class SSATestCaseI: public SSATestCase
{
public:
  SSATestCaseI(MPI_Comm com, PISMConfig &c): 
                 SSATestCase(com,c)
  { };

protected:
  virtual PetscErrorCode initializeGrid(int Mx,int My);

  virtual PetscErrorCode initializeSSAModel();

  virtual PetscErrorCode initializeSSACoefficients();

  virtual PetscErrorCode exactSolution(int i, int j, 
    double x, double y, double *u, double *v );

};

const double m_schoof = 10; // (pure number)
const double L_schoof = 40e3; // meters
const double aspect_schoof = 0.05; // (pure)
const double H0_schoof = aspect_schoof * L_schoof; 
                                       // = 2000 m THICKNESS
const double B_schoof = 3.7e8; // Pa s^{1/3}; hardness 
                                     // given on p. 239 of Schoof; why so big?
const double p_schoof = 4.0/3.0; // = 1 + 1/n


PetscErrorCode SSATestCaseI::initializeGrid(int Mx,int My)
{
  double Ly = 3*L_schoof;  // 300.0 km half-width (L=40.0km in Schoof's choice of variables)
  double Lx = PetscMax(60.0e3, ((Mx - 1) / 2) * (2.0 * Ly / (My - 1)) );
  init_shallow_grid(grid,Lx,Ly,Mx,My,NONE);
  return 0;
}


PetscErrorCode SSATestCaseI::initializeSSAModel()
{
  enthalpyconverter = new EnthalpyConverter(config);

  config.set_flag("do_pseudo_plastic_till", false);

  config.set_string("ssa_flow_law", "isothermal_glen");
  config.set_double("ice_softness", pow(B_schoof, -config.get("Glen_exponent")));

  return 0;
}

PetscErrorCode SSATestCaseI::initializeSSACoefficients()
{
  PetscErrorCode ierr;

  ierr = bc_mask.set(0); CHKERRQ(ierr);
  ierr = thickness.set(H0_schoof); CHKERRQ(ierr);

  // ssa->strength_extension->set_min_thickness(2*H0_schoof);

  // The finite difference code uses the following flag to treat the non-periodic grid correctly.
  config.set_flag("compute_surf_grad_inward_ssa", true);
  config.set_double("epsilon_ssa", 0.0);  // don't use this lower bound

  ierr = tauc.begin_access(); CHKERRQ(ierr);

  double standard_gravity = config.get("standard_gravity"),
    ice_rho = config.get("ice_density");

  for (int i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (int j=grid.ys; j<grid.ys+grid.ym; ++j) {
      const double y = grid.y[j];
      const double theta = atan(0.001);   /* a slope of 1/1000, a la Siple streams */
      const double f = ice_rho * standard_gravity * H0_schoof * tan(theta);
      tauc(i,j) = f * pow(PetscAbs(y / L_schoof), m_schoof);
    }
  }
  ierr = tauc.end_access(); CHKERRQ(ierr);
  ierr = tauc.update_ghosts(); CHKERRQ(ierr);



  ierr = vel_bc.begin_access(); CHKERRQ(ierr);
  ierr = bc_mask.begin_access(); CHKERRQ(ierr);
  ierr = surface.begin_access(); CHKERRQ(ierr);
  ierr = bed.begin_access(); CHKERRQ(ierr);
  for (int i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (int j=grid.ys; j<grid.ys+grid.ym; ++j) {
      double junk, myu, myv;
      const double myx = grid.x[i], myy=grid.y[j];
      // eval exact solution; will only use exact vels if at edge
      exactI(m_schoof, myx, myy, &(bed(i,j)), &junk, &myu, &myv); 
      surface(i,j) = bed(i,j) + H0_schoof;

      bool edge = ( (j == 0) || (j == grid.My - 1) ) || ( (i==0) || (i==grid.Mx-1) );
      if (edge) {
        bc_mask(i,j) = 1;
        vel_bc(i,j).u = myu;
        vel_bc(i,j).v = myv;
      }
    }
  }
  ierr = vel_bc.end_access(); CHKERRQ(ierr);
  ierr = bc_mask.end_access(); CHKERRQ(ierr);
  ierr = surface.end_access(); CHKERRQ(ierr);
  ierr = bed.end_access(); CHKERRQ(ierr);

  // communicate what we have set
  ierr = surface.update_ghosts(); CHKERRQ(ierr);

  ierr = bed.update_ghosts(); CHKERRQ(ierr);

  ierr = bc_mask.update_ghosts(); CHKERRQ(ierr);

  ierr = vel_bc.update_ghosts(); CHKERRQ(ierr);


  ierr = ssa->set_boundary_conditions(bc_mask, vel_bc); CHKERRQ(ierr); 

  return 0;
}


PetscErrorCode SSATestCaseI::exactSolution(int /*i*/, int /*j*/, 
                                           double x, double y,
                                           double *u, double *v)
{
  double junk1, junk2;
  exactI(m_schoof, x,y, &junk1, &junk2,u,v); 
  return 0;
}


int main(int argc, char *argv[]) {
  PetscErrorCode  ierr;

  MPI_Comm    com;

  ierr = PetscInitialize(&argc, &argv, PETSC_NULL, help); CHKERRQ(ierr);

  com = PETSC_COMM_WORLD;
  
  /* This explicit scoping forces destructors to be called before PetscFinalize() */
  {  
    PISMUnitSystem unit_system(NULL);
    PISMConfig config(com, "pism_config", unit_system),
      overrides(com, "pism_overrides", unit_system);
    ierr = init_config(com, config, overrides); CHKERRQ(ierr);

    ierr = setVerbosityLevel(5); CHKERRQ(ierr);

    PetscBool usage_set, help_set;
    ierr = PetscOptionsHasName(PETSC_NULL, "-usage", &usage_set); CHKERRQ(ierr);
    ierr = PetscOptionsHasName(PETSC_NULL, "-help", &help_set); CHKERRQ(ierr);
    if ((usage_set==PETSC_TRUE) || (help_set==PETSC_TRUE)) {
      PetscPrintf(com,
                  "\n"
                  "usage of SSA_TESTi:\n"
                  "  run ssa_testi -Mx <number> -My <number> -ssa_method <fd|fem>\n"
                  "\n");
    }

    // Parameters that can be overridden by command line options
    int Mx=11;
    int My=61;
    std::string output_file = "ssa_test_i.nc";

    std::set<std::string> ssa_choices;
    ssa_choices.insert("fem");
    ssa_choices.insert("fd");
    std::string driver = "fem";

    ierr = PetscOptionsBegin(com, "", "SSA_TESTI options", ""); CHKERRQ(ierr);
    {
      bool flag;
      int my_verbosity_level;
      ierr = PISMOptionsInt("-Mx", "Number of grid points in the X direction", 
                                                      Mx, flag); CHKERRQ(ierr);
      ierr = PISMOptionsInt("-My", "Number of grid points in the Y direction", 
                                                      My, flag); CHKERRQ(ierr);
      ierr = PISMOptionsList(com, "-ssa_method", "Algorithm for computing the SSA solution",
                             ssa_choices, driver, driver, flag); CHKERRQ(ierr);
             
      ierr = PISMOptionsString("-o", "Set the output file name", 
                                              output_file, flag); CHKERRQ(ierr);
      ierr = PISMOptionsInt("-verbose", "Verbosity level",
                            my_verbosity_level, flag); CHKERRQ(ierr);
      if (flag) setVerbosityLevel(my_verbosity_level);
    }
    ierr = PetscOptionsEnd(); CHKERRQ(ierr);

    // Determine the kind of solver to use.
    SSAFactory ssafactory = NULL;
    if(driver == "fem") ssafactory = SSAFEMFactory;
    else if(driver == "fd") ssafactory = SSAFDFactory;
    else { /* can't happen */ }

    SSATestCaseI testcase(com, config);
    ierr = testcase.init(Mx,My,ssafactory); CHKERRQ(ierr);
    ierr = testcase.run(); CHKERRQ(ierr);
    ierr = testcase.report("I"); CHKERRQ(ierr);
    ierr = testcase.write(output_file); CHKERRQ(ierr);
  }

  ierr = PetscFinalize(); CHKERRQ(ierr);
  return 0;
}
