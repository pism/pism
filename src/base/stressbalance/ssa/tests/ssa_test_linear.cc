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

/* This file implements a test case for the ssa: linear flow. The rheology is
   linear (i.e. n=1 in the Glen flow law) and the basal shear stress is also
   linear viscous flow. The geometry consists of a constant surface slope in
   the positive x-direction, and dirichlet conditions leading to an exponential
   solution are imposed along the entire boundary.
*/


static char help[] =
  "\nSSA_TEST_EXP\n"
  "  Testing program for the finite element implementation of the SSA.\n"
  "  Does a time-independent calculation.  Does not run IceModel or a derived\n"
  "  class thereof.Also may be used in a PISM\n"
  "  software (regression) test.\n\n";

#include "pism_const.hh"
#include "iceModelVec.hh"
#include "flowlaws.hh" // IceFlowLaw
#include "basal_resistance.hh" // IceBasalResistancePlasticLaw
#include "PIO.hh"
#include "NCVariable.hh"
#include "SSAFEM.hh"
#include "SSAFD.hh"
#include "exactTestsIJ.h"
#include "SSATestCase.hh"
#include <math.h>
#include "pism_options.hh"

#include "PetscInitializer.hh"
#include "error_handling.hh"

using namespace pism;

class SSATestCaseExp: public SSATestCase
{
public:
  SSATestCaseExp(MPI_Comm com, Config &c)
    : SSATestCase(com, c)
  {
    UnitSystem s = c.get_unit_system();

    L     = s.convert(50, "km", "m"); // 50km half-width
    H0    = 500;                      // meters
    dhdx  = 0.005;                    // pure number
    nu0   = s.convert(30.0, "MPa year", "Pa s");
    tauc0 = 1.e4;               // 1kPa
  }

protected:
  virtual PetscErrorCode initializeGrid(int Mx,int My);

  virtual PetscErrorCode initializeSSAModel();

  virtual PetscErrorCode initializeSSACoefficients();

  virtual PetscErrorCode exactSolution(int i, int j,
    double x, double y, double *u, double *v);

  double L, H0, dhdx, nu0, tauc0;
};


PetscErrorCode SSATestCaseExp::initializeGrid(int Mx,int My)
{
  double Lx=L, Ly = L;
  grid = IceGrid::Shallow(m_com, config, Lx, Ly,
                          0.0, 0.0, // center: (x0,y0)
                          Mx, My, NOT_PERIODIC);
  return 0;
}


PetscErrorCode SSATestCaseExp::initializeSSAModel()
{
  // Use a pseudo-plastic law with linear till
  config.set_flag("do_pseudo_plastic_till", true);
  config.set_double("pseudo_plastic_q", 1.0);

  // The following is irrelevant because we will force linear rheology later.
  enthalpyconverter = new EnthalpyConverter(config);

  return 0;
}

PetscErrorCode SSATestCaseExp::initializeSSACoefficients()
{

  // Force linear rheology
  ssa->strength_extension->set_notional_strength(nu0 * H0);
  ssa->strength_extension->set_min_thickness(4000*10);

  // The finite difference code uses the following flag to treat the non-periodic grid correctly.
  config.set_flag("compute_surf_grad_inward_ssa", true);

  // Set constants for most coefficients.
  thickness.set(H0);
  surface.set(H0);
  bed.set(0.);
  // double threshold_velocity = config.get("pseudo_plastic_uthreshold", "m/year", "m/second");
  // double tauc0 = 4*nu0*H0*threshold_velocity*log(2)*log(2)/(4*L*L);
  // printf("tauc0=%g\n",tauc0);
  tauc.set(tauc0);


  // Set boundary conditions (Dirichlet all the way around).
  bc_mask.set(0.0);

  IceModelVec::AccessList list;
  list.add(vel_bc);
  list.add(bc_mask);
  for (Points p(*grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    double myu, myv;
    const double myx = grid->x(i), myy=grid->y(j);

    bool edge = ((j == 0) || (j == (int)grid->My() - 1) ||
                 (i == 0) || (i == (int)grid->Mx() - 1));
    if (edge) {
      bc_mask(i,j) = 1;
      exactSolution(i,j,myx,myy,&myu,&myv);
      vel_bc(i,j).u = myu;
      vel_bc(i,j).v = myv;
    }
  }

  vel_bc.update_ghosts();
  bc_mask.update_ghosts();

  ssa->set_boundary_conditions(bc_mask, vel_bc);

  return 0;
}


PetscErrorCode SSATestCaseExp::exactSolution(int /*i*/, int /*j*/,
                                             double x, double /*y*/,
                                             double *u, double *v)
{
  double tauc_threshold_velocity = config.get("pseudo_plastic_uthreshold",
                                                   "m/year", "m/second");
  double v0 = grid->convert(100.0, "m/year", "m/second");
  // double alpha=log(2.)/(2*L);
  double alpha = sqrt((tauc0/tauc_threshold_velocity) / (4*nu0*H0));
  *u = v0*exp(-alpha*(x-L));
  *v = 0;
  return 0;
}


int main(int argc, char *argv[]) {
  PetscErrorCode  ierr;

  MPI_Comm com = MPI_COMM_WORLD;  // won't be used except for rank,size

  PetscInitializer petsc(argc, argv, help);

  com = PETSC_COMM_WORLD;

  /* This explicit scoping forces destructors to be called before PetscFinalize() */
  try {
    UnitSystem unit_system;
    Config config(com, "pism_config", unit_system),
      overrides(com, "pism_overrides", unit_system);
    init_config(com, config, overrides);

    setVerbosityLevel(5);

    PetscBool usage_set, help_set;
    ierr = PetscOptionsHasName(NULL, "-usage", &usage_set);
    PISM_PETSC_CHK(ierr, "PetscOptionsHasName");
    ierr = PetscOptionsHasName(NULL, "-help", &help_set);
    PISM_PETSC_CHK(ierr, "PetscOptionsHasName");
    if ((usage_set==PETSC_TRUE) || (help_set==PETSC_TRUE)) {
      PetscPrintf(com,
                  "\n"
                  "usage:\n"
                  "  run %s -Mx <number> -My <number> -ssa_method <fd|fem>\n"
                  "\n",argv[0]);
    }

    // Parameters that can be overridden by command line options
    int Mx=61;
    int My=61;
    std::string output_file = "ssa_test_linear.nc";

    std::set<std::string> ssa_choices;
    ssa_choices.insert("fem");
    ssa_choices.insert("fd");
    std::string driver = "fem";

    ierr = PetscOptionsBegin(com, "", "SSA_TEST_LINEAR options", "");
    PISM_PETSC_CHK(ierr, "PetscOptionsBegin");
    {
      bool flag;
      int my_verbosity_level;
      OptionsInt("-Mx", "Number of grid points in the X direction",
                 Mx, flag);
      OptionsInt("-My", "Number of grid points in the Y direction",
                 My, flag);
      OptionsList("-ssa_method", "Algorithm for computing the SSA solution",
                  ssa_choices, driver, driver, flag);

      OptionsString("-o", "Set the output file name",
                    output_file, flag);

      OptionsInt("-verbose", "Verbosity level",
                 my_verbosity_level, flag);
      if (flag) {
        setVerbosityLevel(my_verbosity_level);
      }
    }
    ierr = PetscOptionsEnd();
    PISM_PETSC_CHK(ierr, "PetscOptionsEnd");

    // Determine the kind of solver to use.
    SSAFactory ssafactory = NULL;
    if (driver == "fem") {
      ssafactory = SSAFEMFactory;
    } else if (driver == "fd") {
      ssafactory = SSAFDFactory;
    } else {
      /* can't happen */
    }

    SSATestCaseExp testcase(com,config);
    testcase.init(Mx,My,ssafactory);
    testcase.run();
    testcase.report("linear");
    testcase.write(output_file);
  }
  catch (...) {
    handle_fatal_errors(com);
  }

  return 0;
}
