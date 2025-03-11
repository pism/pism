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

/* This file implements a test case for the ssa: plug flow. The geometry
   consists of a constant surface slope in the positive x-direction, and the
   ice is pinned on the y-boundaries. There is no basal shear stress, and hence
   the the only nonzero terms in the SSA are the "p-laplacian" and the driving
   stress.
*/

static char help[] =
  "\nSSA_TEST_PLUG\n"
  "  Testing program for the finite element implementation of the SSA.\n"
  "  Does a time-independent calculation.  Does not run IceModel or a derived\n"
  "  class thereof.\n\n";

#include <cmath>

#include "pism/stressbalance/ssa/SSAFD.hh"
#include "pism/stressbalance/ssa/SSAFEM.hh"
#include "pism/stressbalance/ssa/SSATestCase.hh"
#include "pism/util/Context.hh"
#include "pism/util/error_handling.hh"
#include "pism/util/petscwrappers/PetscInitializer.hh"
#include "pism/util/pism_options.hh"

namespace pism {
namespace stressbalance {

std::shared_ptr<Grid> ssa_test_plug_grid(std::shared_ptr<Context> ctx, int Mx, int My) {
  return SSATestCase::grid(ctx, Mx, My, 50e3, 50e3, grid::CELL_CORNER, grid::NOT_PERIODIC);
}

class SSATestCasePlug : public SSATestCase {
public:
  SSATestCasePlug(std::shared_ptr<SSA> ssa)
    : SSATestCase(ssa) {
  }

  static const double H0;    // Thickness, meters
  static const double L;     // Half-width, meters
  static const double dhdx;  // surface slope, pure number
  static const double tauc0; // zero basal shear stress, Pa
  static const double B0;    // Pa s^{1/3}; hardness given on p. 239 of Schoof; why so big?

protected:
  virtual void initializeSSACoefficients();

  virtual void exactSolution(int i, int j, double x, double y, double *u, double *v);
};

const double SSATestCasePlug::H0    = 2000.0;
const double SSATestCasePlug::L     = 50e3;
const double SSATestCasePlug::dhdx  = 0.001;
const double SSATestCasePlug::tauc0 = 0.0;
const double SSATestCasePlug::B0    = 3.7e8;

void SSATestCasePlug::initializeSSACoefficients() {

  // Ensure we never use the strength extension.
  m_ssa->strength_extension->set_min_thickness(H0 / 2);

  // Set constant coefficients.
  m_geometry.ice_thickness.set(H0);
  m_tauc.set(tauc0);

  // Set boundary conditions (Dirichlet all the way around).
  m_bc_mask.set(0.0);

  array::AccessScope list{ &m_bc_values, &m_bc_mask, &m_geometry.bed_elevation,
                           &m_geometry.ice_surface_elevation };

  for (auto p = m_grid->points(); p; p.next()) {
    const int i = p.i(), j = p.j();

    double myu, myv;
    const double myx = m_grid->x(i), myy = m_grid->y(j);

    m_geometry.bed_elevation(i, j)         = -myx * (dhdx);
    m_geometry.ice_surface_elevation(i, j) = m_geometry.bed_elevation(i, j) + H0;

    bool edge =
        ((j == 0) || (j == (int)m_grid->My() - 1) || (i == 0) || (i == (int)m_grid->Mx() - 1));
    if (edge) {
      m_bc_mask(i, j) = 1;
      exactSolution(i, j, myx, myy, &myu, &myv);
      m_bc_values(i, j).u = myu;
      m_bc_values(i, j).v = myv;
    }
  }

  m_bc_values.update_ghosts();
  m_bc_mask.update_ghosts();
  m_geometry.bed_elevation.update_ghosts();
  m_geometry.ice_surface_elevation.update_ghosts();
}

void SSATestCasePlug::exactSolution(int /*i*/, int /*j*/, double /*x*/, double y, double *u,
                                    double *v) {
  double earth_grav = m_config->get_number("constants.standard_gravity"),
         ice_rho    = m_config->get_number("constants.ice.density");
  double f          = ice_rho * earth_grav * H0 * dhdx;
  double ynd        = y / L;

  *u = 0.5 * pow(f, 3) * pow(L, 4) / pow(B0 * H0, 3) * (1 - pow(ynd, 4));
  *v = 0;
}

} // end of namespace stressbalance
} // end of namespace pism

int main(int argc, char *argv[]) {

  using namespace pism;
  using namespace pism::stressbalance;

  MPI_Comm com = MPI_COMM_WORLD; // won't be used except for rank,size
  petsc::Initializer petsc(argc, argv, help);

  com = PETSC_COMM_WORLD;

  /* This explicit scoping forces destructors to be called before PetscFinalize() */
  try {
    std::shared_ptr<Context> ctx = context_from_options(com, "ssa_test_plug");
    auto config           = ctx->config();

    std::string usage = "\n"
                        "usage of SSA_TEST_PLUG:\n"
                        "  run ssa_test_plug -Mx <number> -My <number> -ssa_method <fd|fem>\n"
                        "\n";

    bool stop = show_usage_check_req_opts(*ctx->log(), "ssa_test_plug", {}, usage);

    if (stop) {
      return 0;
    }

    // Parameters that can be overridden by command line options

    unsigned int Mx = config->get_number("grid.Mx");
    unsigned int My = config->get_number("grid.My");

    auto method      = config->get_string("stress_balance.ssa.method");
    auto output_file = config->get_string("output.file");

    bool write_output = config->get_string("output.size") != "none";

    // we have to set parameters *before* `ssa` is allocated
    //
    // Use constant hardness
    config->set_string("stress_balance.ssa.flow_law", "isothermal_glen");
    config->set_number(
        "flow_law.isothermal_Glen.ice_softness",
        pow(SSATestCasePlug::B0, -config->get_number("stress_balance.ssa.Glen_exponent")));

    // The finite difference code uses the following flag to treat the non-periodic grid correctly.
    config->set_flag("stress_balance.ssa.compute_surface_gradient_inward", true);
    config->set_number("stress_balance.ssa.epsilon", 0.0);

    auto grid = ssa_test_plug_grid(ctx, Mx, My);
    SSATestCasePlug testcase(SSATestCase::solver(grid, method));
    testcase.init();
    testcase.run();
    testcase.report("plug");
    if (write_output) {
      testcase.write(output_file);
    }
  } catch (...) {
    handle_fatal_errors(com);
    return 1;
  }

  return 0;
}
