// Copyright (C) 2010--2018 Ed Bueler, Constantine Khroulev, and David Maxwell
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

#include "pism/basalstrength/basal_resistance.hh" // IceBasalResistancePlasticLaw
#include "pism/stressbalance/ssa/SSAFD.hh"
#include "pism/stressbalance/ssa/SSAFEM.hh"
#include "pism/stressbalance/ssa/SSATestCase.hh"
#include "pism/util/Context.hh"
#include "pism/util/VariableMetadata.hh"
#include "pism/util/error_handling.hh"
#include "pism/util/iceModelVec.hh"
#include "pism/util/io/PIO.hh"
#include "pism/util/petscwrappers/PetscInitializer.hh"
#include "pism/util/pism_utilities.hh"
#include "pism/util/pism_options.hh"
#include "pism/verification/tests/exactTestsIJ.h"

namespace pism {
namespace stressbalance {

class SSATestCaseExp: public SSATestCase
{
public:
  SSATestCaseExp(Context::Ptr ctx, int Mx, int My, SSAFactory ssafactory)
    : SSATestCase(ctx, Mx, My, 50e3, 50e3, CELL_CORNER, NOT_PERIODIC) {
    L     = units::convert(ctx->unit_system(), 50, "km", "m"); // 50km half-width
    H0    = 500;                      // meters
    dhdx  = 0.005;                    // pure number
    nu0   = units::convert(ctx->unit_system(), 30.0, "MPa year", "Pa s");
    tauc0 = 1.e4;               // 1kPa

    // Use a pseudo-plastic law with linear till
    m_config->set_boolean("basal_resistance.pseudo_plastic.enabled", true);
    m_config->set_double("basal_resistance.pseudo_plastic.q", 1.0);

    // The following is irrelevant because we will force linear rheology later.
    m_enthalpyconverter = EnthalpyConverter::Ptr(new EnthalpyConverter(*m_config));

    m_ssa = ssafactory(m_grid);
  }

protected:
  virtual void initializeSSACoefficients();

  virtual void exactSolution(int i, int j,
    double x, double y, double *u, double *v);

  double L, H0, dhdx, nu0, tauc0;
};

void SSATestCaseExp::initializeSSACoefficients() {

  // Force linear rheology
  m_ssa->strength_extension->set_notional_strength(nu0 * H0);
  m_ssa->strength_extension->set_min_thickness(4000*10);

  // The finite difference code uses the following flag to treat the non-periodic grid correctly.
  m_config->set_boolean("stress_balance.ssa.compute_surface_gradient_inward", true);

  // Set constants for most coefficients.
  m_geometry.ice_thickness.set(H0);
  m_geometry.ice_surface_elevation.set(H0);
  m_geometry.bed_elevation.set(0.0);
  m_tauc.set(tauc0);


  // Set boundary conditions (Dirichlet all the way around).
  m_bc_mask.set(0.0);

  IceModelVec::AccessList list{&m_bc_values, &m_bc_mask};

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
}


void SSATestCaseExp::exactSolution(int /*i*/, int /*j*/, double x, double /*y*/,
                                   double *u, double *v) {
  double tauc_threshold_velocity = m_config->get_double("basal_resistance.pseudo_plastic.u_threshold",
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

  com = PETSC_COMM_WORLD;

  /* This explicit scoping forces destructors to be called before PetscFinalize() */
  try {
    Context::Ptr ctx = context_from_options(com, "ssa_test_linear");
    Config::Ptr config = ctx->config();

    std::string usage = "\n"
      "usage:\n"
      "  run ssa_test_linear -Mx <number> -My <number> -ssa_method <fd|fem>\n"
      "\n";

    bool stop = show_usage_check_req_opts(*ctx->log(), "ssa_test_linear", {}, usage);

    if (stop) {
      return 0;
    }

    // Parameters that can be overridden by command line options
    unsigned int Mx = config->get_double("grid.Mx");
    unsigned int My = config->get_double("grid.My");

    auto method      = config->get_string("stress_balance.ssa.method");
    auto output_file = config->get_string("output.file_name");

    // Determine the kind of solver to use.
    SSAFactory ssafactory = NULL;
    if (method == "fem") {
      ssafactory = SSAFEMFactory;
    } else if (method == "fd") {
      ssafactory = SSAFDFactory;
    } else {
      /* can't happen */
    }

    SSATestCaseExp testcase(ctx, Mx, My, ssafactory);
    testcase.init();
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
