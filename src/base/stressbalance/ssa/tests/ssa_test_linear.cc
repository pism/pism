// Copyright (C) 2010--2015 Ed Bueler, Constantine Khroulev, and David Maxwell
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
  m_grid = IceGrid::Shallow(m_com, m_config, Lx, Ly,
                          0.0, 0.0, // center: (x0,y0)
                          Mx, My, NOT_PERIODIC);
  return 0;
}


PetscErrorCode SSATestCaseExp::initializeSSAModel()
{
  // Use a pseudo-plastic law with linear till
  m_config.set_flag("do_pseudo_plastic_till", true);
  m_config.set_double("pseudo_plastic_q", 1.0);

  // The following is irrelevant because we will force linear rheology later.
  m_enthalpyconverter = new EnthalpyConverter(m_config);

  return 0;
}

PetscErrorCode SSATestCaseExp::initializeSSACoefficients()
{

  // Force linear rheology
  m_ssa->strength_extension->set_notional_strength(nu0 * H0);
  m_ssa->strength_extension->set_min_thickness(4000*10);

  // The finite difference code uses the following flag to treat the non-periodic grid correctly.
  m_config.set_flag("compute_surf_grad_inward_ssa", true);

  // Set constants for most coefficients.
  m_thickness.set(H0);
  m_surface.set(H0);
  m_bed.set(0.);
  // double threshold_velocity = config.get("pseudo_plastic_uthreshold", "m/year", "m/second");
  // double tauc0 = 4*nu0*H0*threshold_velocity*log(2)*log(2)/(4*L*L);
  // printf("tauc0=%g\n",tauc0);
  m_tauc.set(tauc0);


  // Set boundary conditions (Dirichlet all the way around).
  m_bc_mask.set(0.0);

  IceModelVec::AccessList list;
  list.add(m_vel_bc);
  list.add(m_bc_mask);
  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    double myu, myv;
    const double myx = m_grid->x(i), myy=m_grid->y(j);

    bool edge = ((j == 0) || (j == (int)m_grid->My() - 1) ||
                 (i == 0) || (i == (int)m_grid->Mx() - 1));
    if (edge) {
      m_bc_mask(i,j) = 1;
      exactSolution(i,j,myx,myy,&myu,&myv);
      m_vel_bc(i,j).u = myu;
      m_vel_bc(i,j).v = myv;
    }
  }

  m_vel_bc.update_ghosts();
  m_bc_mask.update_ghosts();

  m_ssa->set_boundary_conditions(m_bc_mask, m_vel_bc);

  return 0;
}


PetscErrorCode SSATestCaseExp::exactSolution(int /*i*/, int /*j*/,
                                             double x, double /*y*/,
                                             double *u, double *v)
{
  double tauc_threshold_velocity = m_config.get("pseudo_plastic_uthreshold",
                                                   "m/year", "m/second");
  double v0 = m_grid->convert(100.0, "m/year", "m/second");
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
    if ((usage_set==true) || (help_set==true)) {
      PetscPrintf(com,
                  "\n"
                  "usage:\n"
                  "  run %s -Mx <number> -My <number> -ssa_method <fd|fem>\n"
                  "\n",argv[0]);
    }

    // Parameters that can be overridden by command line options
    options::Integer Mx("-Mx", "Number of grid points in the X direction", 61);
    options::Integer My("-My", "Number of grid points in the Y direction", 61);

    options::Keyword method("-ssa_method", "Algorithm for computing the SSA solution",
                            "fem,fd", "fem");

    options::String output_file("-o", "Set the output file name", "ssa_test_linear.nc");

    options::Integer my_verbosity_level("-verbose", "Verbosity level", 2);
    if (my_verbosity_level.is_set()) {
      setVerbosityLevel(my_verbosity_level);
    }

    // Determine the kind of solver to use.
    SSAFactory ssafactory = NULL;
    if (method == "fem") {
      ssafactory = SSAFEMFactory;
    } else if (method == "fd") {
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
