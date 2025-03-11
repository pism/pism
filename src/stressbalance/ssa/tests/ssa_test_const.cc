// Copyright (C) 2010--2018, 2021, 2022, 2024, 2025 Ed Bueler, Constantine Khroulev, and David Maxwell
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

/* This file implements a test case for the ssa: constant flow. The rheology is
   nonlinear (i.e. n=3 in the Glen flow law) and the basal shear stress is a
   nonlinear function of velocity (peseudo-plastic flow with parameter q
   specified at runtime).

   The geometry consists of a constant surface slope in the positive
   x-direction, and a constant velocity is specified as a Dirichlet condition
   on the boundary that should lead to a constant solution in the interior.
   Because the solution is constant, the nonzero terms in the SSA are only the
   basal shear stress and the driving stress.
 */

#include <memory>
static char help[] =
  "\nSSA_TEST_CONST\n"
  "  Testing program for the finite element implementation of the SSA.\n"
  "  Does a time-independent calculation.  Does not run IceModel or a derived\n"
  "  class thereof.Also may be used in a PISM\n"
  "  software (regression) test.\n\n";

#include <cmath>

#include "pism/basalstrength/basal_resistance.hh" // IceBasalResistancePlasticLaw
#include "pism/stressbalance/ssa/SSAFD.hh"
#include "pism/stressbalance/ssa/SSAFEM.hh"
#include "pism/stressbalance/ssa/SSATestCase.hh"
#include "pism/util/Mask.hh"
#include "pism/util/Context.hh"
#include "pism/util/VariableMetadata.hh"
#include "pism/util/error_handling.hh"
#include "pism/util/io/File.hh"
#include "pism/util/petscwrappers/PetscInitializer.hh"
#include "pism/util/pism_utilities.hh"
#include "pism/util/pism_options.hh"
#include "pism/verification/tests/exactTestsIJ.h"

namespace pism {
namespace stressbalance {

std::shared_ptr<Grid> ssa_test_const_grid(std::shared_ptr<Context> ctx, int Mx, int My) {
  return SSATestCase::grid(ctx, Mx, My, 50e3, 50e3, grid::CELL_CORNER, grid::NOT_PERIODIC);
}

class SSATestCaseConst: public SSATestCase
{
public:
  SSATestCaseConst(std::shared_ptr<Context> ctx, std::shared_ptr<SSA> ssa) : SSATestCase(ssa) {
    m_L     = units::convert(m_sys, 50.0, "km", "m"); // 50km half-width
    m_H0    = 500;                                    // m
    m_dhdx  = 0.005;                                  // pure number
    m_nu0   = units::convert(m_sys, 30.0, "MPa year", "Pa s");
    m_tauc0 = 1.e4; // Pa

    auto config = ctx->config();
    m_basal_q = config->get_number("basal_resistance.pseudo_plastic.q");
  };

protected:
  void initializeSSACoefficients();

  void exactSolution(int i, int j, double x, double y, double *u, double *v);

  double m_basal_q;
  double m_L;
  double m_H0;
  double m_dhdx;
  double m_nu0;
  double m_tauc0;
};

void SSATestCaseConst::initializeSSACoefficients() {

  // Force linear rheology
  m_ssa->strength_extension->set_notional_strength(m_nu0 * m_H0);
  m_ssa->strength_extension->set_min_thickness(0.5*m_H0);

  // Set constant thickness, tauc
  m_bc_mask.set(0);
  m_geometry.ice_thickness.set(m_H0);
  m_tauc.set(m_tauc0);

  array::AccessScope list{&m_bc_values, &m_bc_mask,
      &m_geometry.bed_elevation, &m_geometry.ice_surface_elevation};

  for (auto p = m_grid->points(); p; p.next()) {
    const int i = p.i(), j = p.j();

    double u, v;
    const double x = m_grid->x(i), y=m_grid->y(j);

    m_geometry.bed_elevation(i, j) = -x*(m_dhdx);
    m_geometry.ice_surface_elevation(i, j) = m_geometry.bed_elevation(i, j) + m_H0;

    bool edge = ((j == 0) || (j == (int)m_grid->My() - 1) ||
                 (i == 0) || (i == (int)m_grid->Mx() - 1));
    if (edge) {
      m_bc_mask(i,j) = 1;
      exactSolution(i, j, x, y, &u, &v);
      m_bc_values(i, j).u = u;
      m_bc_values(i, j).v = v;
    }
  }

  m_bc_values.update_ghosts();
  m_bc_mask.update_ghosts();
  m_geometry.bed_elevation.update_ghosts();
  m_geometry.ice_surface_elevation.update_ghosts();
}


void SSATestCaseConst::exactSolution(int /*i*/, int /*j*/,
                                     double /*x*/, double /*y*/,
                                     double *u, double *v) {
  double earth_grav = m_config->get_number("constants.standard_gravity"),
    tauc_threshold_velocity = m_config->get_number("basal_resistance.pseudo_plastic.u_threshold",
                                                   "m second-1"),
    ice_rho = m_config->get_number("constants.ice.density");

  *u = pow(ice_rho * earth_grav * m_H0 * m_dhdx / m_tauc0, 1./m_basal_q)*tauc_threshold_velocity;
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
    std::shared_ptr<Context> ctx = context_from_options(com, "ssa_test_const");
    auto config = ctx->config();

    std::string usage = "\n"
      "usage of SSA_TEST_CONST:\n"
      "  run ssa_test_const -Mx <number> -My <number> -ssa_method <fd|fem>\n"
      "\n";

    bool stop = show_usage_check_req_opts(*ctx->log(), "ssa_test_const", {}, usage);

    if (stop) {
      return 0;
    }

    // Parameters that can be overridden by command line options
    unsigned int Mx = config->get_number("grid.Mx");
    unsigned int My = config->get_number("grid.My");

    double basal_q = 1.0;

    auto method = config->get_string("stress_balance.ssa.method");
    auto output_file = config->get_string("output.file");

    bool write_output = config->get_string("output.size") != "none";

    config->set_flag("basal_resistance.pseudo_plastic.enabled", true);
    config->set_number("basal_resistance.pseudo_plastic.q", basal_q);

    // Use a pseudo-plastic law with a constant q determined at run time
    config->set_flag("basal_resistance.pseudo_plastic.enabled", true);

    // The finite difference code uses the following flag to treat the non-periodic grid correctly.
    config->set_flag("stress_balance.ssa.compute_surface_gradient_inward", true);

    auto grid = ssa_test_const_grid(ctx, Mx, My);
    SSATestCaseConst testcase(ctx, SSATestCase::solver(grid, method));
    testcase.init();
    testcase.run();
    testcase.report("const");
    if (write_output) {
      testcase.write(output_file);
    }
  }
  catch (...) {
    handle_fatal_errors(com);
    return 1;
  }

  return 0;
}
