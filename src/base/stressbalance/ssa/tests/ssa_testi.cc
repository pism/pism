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

#include "PetscInitializer.hh"
#include "error_handling.hh"

using namespace pism;

class SSATestCaseI: public SSATestCase
{
public:
  SSATestCaseI(MPI_Comm com, Config &c):
                 SSATestCase(com,c)
  { };

protected:
  virtual PetscErrorCode initializeGrid(int Mx,int My);

  virtual PetscErrorCode initializeSSAModel();

  virtual PetscErrorCode initializeSSACoefficients();

  virtual PetscErrorCode exactSolution(int i, int j,
    double x, double y, double *u, double *v);

};

const double m_schoof = 10; // (pure number)
const double L_schoof = 40e3; // meters
const double aspect_schoof = 0.05; // (pure)
const double H0_schoof = aspect_schoof * L_schoof;
                                       // = 2000 m THICKNESS
const double B_schoof = 3.7e8; // Pa s^{1/3}; hardness
                                     // given on p. 239 of Schoof; why so big?

PetscErrorCode SSATestCaseI::initializeGrid(int Mx,int My)
{
  double Ly = 3*L_schoof;  // 300.0 km half-width (L=40.0km in Schoof's choice of variables)
  double Lx = std::max(60.0e3, ((Mx - 1) / 2) * (2.0 * Ly / (My - 1)));
  grid = IceGrid::Shallow(m_com, config, Lx, Ly,
                          0.0, 0.0, // center: (x0,y0)
                          Mx, My, NOT_PERIODIC);
  return 0;
}


PetscErrorCode SSATestCaseI::initializeSSAModel()
{
  enthalpyconverter = new EnthalpyConverter(config);

  config.set_flag("do_pseudo_plastic_till", false);

  config.set_string("ssa_flow_law", "isothermal_glen");
  config.set_double("ice_softness", pow(B_schoof, -config.get("ssa_Glen_exponent")));

  return 0;
}

PetscErrorCode SSATestCaseI::initializeSSACoefficients()
{

  bc_mask.set(0);
  thickness.set(H0_schoof);

  // ssa->strength_extension->set_min_thickness(2*H0_schoof);

  // The finite difference code uses the following flag to treat the non-periodic grid correctly.
  config.set_flag("compute_surf_grad_inward_ssa", true);
  config.set_double("epsilon_ssa", 0.0);  // don't use this lower bound

  IceModelVec::AccessList list;
  list.add(tauc);

  double standard_gravity = config.get("standard_gravity"),
    ice_rho = config.get("ice_density");

  for (Points p(*grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    const double y = grid->y(j);
    const double theta = atan(0.001);   /* a slope of 1/1000, a la Siple streams */
    const double f = ice_rho * standard_gravity * H0_schoof * tan(theta);
    tauc(i,j) = f * pow(fabs(y / L_schoof), m_schoof);
  }
  tauc.update_ghosts();

  list.add(vel_bc);
  list.add(bc_mask);
  list.add(surface);
  list.add(bed);
  for (Points p(*grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    double junk, myu, myv;
    const double myx = grid->x(i), myy=grid->y(j);
    // eval exact solution; will only use exact vels if at edge
    exactI(m_schoof, myx, myy, &(bed(i,j)), &junk, &myu, &myv);
    surface(i,j) = bed(i,j) + H0_schoof;

    bool edge = ((j == 0) || (j == (int)grid->My() - 1) ||
                 (i == 0) || (i == (int)grid->Mx() - 1));
    if (edge) {
      bc_mask(i,j) = 1;
      vel_bc(i,j).u = myu;
      vel_bc(i,j).v = myv;
    }
  }

  // communicate what we have set
  surface.update_ghosts();
  bed.update_ghosts();
  bc_mask.update_ghosts();
  vel_bc.update_ghosts();

  ssa->set_boundary_conditions(bc_mask, vel_bc);

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

  MPI_Comm com = MPI_COMM_WORLD;

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

    ierr = PetscOptionsBegin(com, "", "SSA_TESTI options", "");
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

    SSATestCaseI testcase(com, config);
    testcase.init(Mx,My,ssafactory);
    testcase.run();
    testcase.report("I");
    testcase.write(output_file);
  }
  catch (...) {
    handle_fatal_errors(com);
  }

  return 0;
}
