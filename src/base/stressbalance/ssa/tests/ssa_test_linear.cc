// Copyright (C) 2010--2016 Ed Bueler, Constantine Khroulev, and David Maxwell
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

#include <cmath>

#include "base/basalstrength/basal_resistance.hh" // IceBasalResistancePlasticLaw
#include "base/stressbalance/ssa/SSAFD.hh"
#include "base/stressbalance/ssa/SSAFEM.hh"
#include "base/stressbalance/ssa/SSATestCase.hh"
#include "base/util/Context.hh"
#include "base/util/VariableMetadata.hh"
#include "base/util/error_handling.hh"
#include "base/util/iceModelVec.hh"
#include "base/util/io/PIO.hh"
#include "base/util/petscwrappers/PetscInitializer.hh"
#include "base/util/pism_const.hh"
#include "base/util/pism_options.hh"
#include "verif/tests/exactTestsIJ.h"

namespace pism {
namespace stressbalance {

class SSATestCaseExp: public SSATestCase
{
public:
  SSATestCaseExp(Context::Ptr ctx)
    : SSATestCase(ctx) {
    L     = units::convert(ctx->unit_system(), 50, "km", "m"); // 50km half-width
    H0    = 500;                      // meters
    dhdx  = 0.005;                    // pure number
    nu0   = units::convert(ctx->unit_system(), 30.0, "MPa year", "Pa s");
    tauc0 = 1.e4;               // 1kPa
  }

protected:
  virtual void initializeGrid(int Mx,int My);

  virtual void initializeSSAModel();

  virtual void initializeSSACoefficients();

  virtual void exactSolution(int i, int j,
    double x, double y, double *u, double *v);

  double L, H0, dhdx, nu0, tauc0;
};


void SSATestCaseExp::initializeGrid(int Mx,int My) {
  double Lx=L, Ly = L;
  m_grid = IceGrid::Shallow(m_ctx, Lx, Ly,
                            0.0, 0.0, // center: (x0,y0)
                            Mx, My, NOT_PERIODIC);
}


void SSATestCaseExp::initializeSSAModel() {
  // Use a pseudo-plastic law with linear till
  m_config->set_boolean("do_pseudo_plastic_till", true);
  m_config->set_double("pseudo_plastic_q", 1.0);

  // The following is irrelevant because we will force linear rheology later.
  m_enthalpyconverter = EnthalpyConverter::Ptr(new EnthalpyConverter(*m_config));
}

void SSATestCaseExp::initializeSSACoefficients() {

  // Force linear rheology
  m_ssa->strength_extension->set_notional_strength(nu0 * H0);
  m_ssa->strength_extension->set_min_thickness(4000*10);

  // The finite difference code uses the following flag to treat the non-periodic grid correctly.
  m_config->set_boolean("compute_surf_grad_inward_ssa", true);

  // Set constants for most coefficients.
  m_thickness.set(H0);
  m_surface.set(H0);
  m_bed.set(0.);
  m_tauc.set(tauc0);


  // Set boundary conditions (Dirichlet all the way around).
  m_bc_mask.set(0.0);

  IceModelVec::AccessList list;
  list.add(m_bc_values);
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
      m_bc_values(i,j).u = myu;
      m_bc_values(i,j).v = myv;
    }
  }

  m_bc_values.update_ghosts();
  m_bc_mask.update_ghosts();

  m_ssa->set_boundary_conditions(m_bc_mask, m_bc_values);
}


void SSATestCaseExp::exactSolution(int /*i*/, int /*j*/, double x, double /*y*/,
                                   double *u, double *v) {
  double tauc_threshold_velocity = m_config->get_double("pseudo_plastic_uthreshold",
                                                        "m second-1");
  double v0 = units::convert(m_sys, 100.0, "m year-1", "m second-1");
  // double alpha=log(2.)/(2*L);
  double alpha = sqrt((tauc0/tauc_threshold_velocity) / (4*nu0*H0));
  *u = v0*exp(-alpha*(x-L));
  *v = 0;
}

} // end of namespace stressbalance
} // end of namespace pism

int main(int argc, char *argv[]) {

  using namespace pism;
  using namespace pism::stressbalance;

  MPI_Comm com = MPI_COMM_WORLD;  // won't be used except for rank,size
  petsc::Initializer petsc(argc, argv, help);
  PetscErrorCode ierr;

  com = PETSC_COMM_WORLD;

  /* This explicit scoping forces destructors to be called before PetscFinalize() */
  try {
    verbosityLevelFromOptions();
    Context::Ptr ctx = context_from_options(com, "ssa_test_linear");
    Config::Ptr config = ctx->config();

    setVerbosityLevel(5);

    bool
      usage_set = options::Bool("-usage", "print usage info"),
      help_set  = options::Bool("-help", "print help info");
    if (usage_set or help_set) {
      ierr = PetscPrintf(com,
                         "\n"
                         "usage:\n"
                         "  run %s -Mx <number> -My <number> -ssa_method <fd|fem>\n"
                         "\n",argv[0]);
      PISM_CHK(ierr, "PetscPrintf");
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

    SSATestCaseExp testcase(ctx);
    testcase.init(Mx,My,ssafactory);
    testcase.run();
    testcase.report("linear");
    testcase.write(output_file);
  }
  catch (...) {
    handle_fatal_errors(com);
    return 1;
  }

  return 0;
}
