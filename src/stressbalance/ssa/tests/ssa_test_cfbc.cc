// Copyright (C) 2010--2018, 2021, 2022, 2023, 2024 Ed Bueler, Constantine Khroulev, and David Maxwell
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
  "\nSSA_TESTCFBC\n"
  "  Testing program for PISM's implementations of the SSA.\n"
  "  Does a time-independent calculation.  Does not run IceModel or a derived\n"
  "  class thereof. Uses the van der Veen flow-line shelf geometry. Also may be\n"
  "  used in a PISM software (regression) test.\n\n";

#include "pism/stressbalance/ssa/SSATestCase.hh"
#include "pism/util/Context.hh"
#include "pism/util/error_handling.hh"
#include "pism/util/petscwrappers/PetscInitializer.hh"
#include "pism/util/pism_options.hh"

#include "pism/stressbalance/ssa/SSAFD.hh"
#include "pism/stressbalance/ssa/SSAFEM.hh"

namespace pism {
namespace stressbalance {

// thickness profile in the van der Veen solution
static double H_exact(double V0, double H0, double C, double x) {
  const double Q0 = V0*H0;
  return pow(4 * C / Q0 * x + 1/pow(H0, 4), -0.25);
}

// velocity profile; corresponds to constant flux
static double u_exact(double V0, double H0, double C, double x) {
  const double Q0 = V0*H0;
  return Q0 / H_exact(V0, H0, C, x);
}

std::shared_ptr<Grid> ssa_test_cfbc_grid(std::shared_ptr<Context> ctx, int Mx, int My) {
  return SSATestCase::grid(ctx, Mx, My, 250e3, 250e3, grid::CELL_CENTER, grid::Y_PERIODIC);
}

class SSATestCaseCFBC: public SSATestCase {
public:
  SSATestCaseCFBC(std::shared_ptr<SSA> ssa)
    : SSATestCase(ssa) {
    m_V0 = units::convert(m_sys, 300.0, "m year-1", "m second-1");
    m_H0 = 600.0;                 // meters
    m_C  = 2.45e-18;

    EnthalpyConverter EC(*m_config);
    // 0.01 water fraction
    m_ice_enthalpy.set(EC.enthalpy(273.15, 0.01, 0.0));
  }

  void write_nuH(const std::string &filename);

protected:
  void initializeSSACoefficients();

  void exactSolution(int i, int j, double x, double y, double *u, double *v);

  double m_V0, //!< grounding line vertically-averaged velocity
    m_H0,      //!< grounding line thickness (meters)
    m_C;       //!< "typical constant ice parameter"
};

void SSATestCaseCFBC::write_nuH(const std::string &filename) {
  auto diagnostics = m_ssa->diagnostics();

  if (diagnostics.find("nuH") != diagnostics.end()) {
    diagnostics["nuH"]->compute()->write(filename);
  }
}

void SSATestCaseCFBC::initializeSSACoefficients() {

  m_tauc.set(0.0);    // irrelevant
  m_geometry.bed_elevation.set(-1000.0); // assures shelf is floating

  array::AccessScope list{&m_geometry.ice_thickness,
      &m_geometry.ice_surface_elevation, &m_bc_mask, &m_bc_values, &m_geometry.cell_type};

  const double x_min = m_grid->x(0);

  for (auto p = m_grid->points(); p; p.next()) {
    const int i = p.i(), j = p.j();

    const double x = m_grid->x(i);

    if (i != (int)m_grid->Mx() - 1) {
      m_geometry.ice_thickness(i, j) = H_exact(m_V0, m_H0, m_C, x - x_min);
    } else {
      m_geometry.ice_thickness(i, j) = 0.0;
    }

    if (i == 0) {
      m_bc_mask(i, j)   = 1;
      m_bc_values(i, j) = {m_V0, 0.0};
    } else {
      m_bc_mask(i, j)   = 0;
      m_bc_values(i, j) = {0.0, 0.0};
    }
  }

  // communicate what we have set
  m_geometry.ensure_consistency(0.0);

  m_bc_mask.update_ghosts();
  m_bc_values.update_ghosts();
}

void SSATestCaseCFBC::exactSolution(int i, int /*j*/,
                                    double x, double /*y*/,
                                    double *u, double *v) {
  const double x_min = m_grid->x(0);

  if (i != (int)m_grid->Mx() - 1) {
    *u = u_exact(m_V0, m_H0, m_C, x - x_min);
  } else {
    *u = 0;
  }

  *v = 0;
}

} // end of namespace stressbalance
} // end of namespace pism

int main(int argc, char *argv[]) {

  using namespace pism;
  using namespace pism::stressbalance;

  MPI_Comm com = MPI_COMM_WORLD;
  petsc::Initializer petsc(argc, argv, help);

  com = PETSC_COMM_WORLD;

  /* This explicit scoping forces destructors to be called before PetscFinalize() */
  try {
    std::shared_ptr<Context> ctx = context_from_options(com, "ssa_test_cfbc");
    Config::Ptr config = ctx->config();

    std::string usage = "\n"
      "usage of SSA_TEST_CFBC:\n"
      "  run ssa_test_cfbc -Mx <number> -My <number>\n"
      "\n";

    bool stop = show_usage_check_req_opts(*ctx->log(), "ssa_test_cfbc", {}, usage);

    if (stop) {
      return 0;
    }

    // Parameters that can be overridden by command line options
    unsigned int Mx = config->get_number("grid.Mx");
    unsigned int My = config->get_number("grid.My");

    auto method = config->get_string("stress_balance.ssa.method");
    auto output_file = config->get_string("output.file");

    bool write_output = config->get_string("output.size") != "none";

    // we have to set parameters *before* `ssa` is allocated
    config->set_number("flow_law.isothermal_Glen.ice_softness",
                       pow(1.9e8, -config->get_number("stress_balance.ssa.Glen_exponent")));
    config->set_flag("stress_balance.ssa.compute_surface_gradient_inward", false);
    config->set_flag("stress_balance.calving_front_stress_bc", true);
    config->set_flag("stress_balance.ssa.fd.flow_line_mode", true);
    config->set_flag("stress_balance.ssa.fd.extrapolate_at_margins", false);
    config->set_string("stress_balance.ssa.flow_law", "isothermal_glen");

    auto grid = ssa_test_cfbc_grid(ctx, Mx, My);
    SSATestCaseCFBC testcase(SSATestCase::solver(grid, method));
    testcase.init();
    testcase.run();
    testcase.report("V");
    if (write_output) {
      testcase.write(output_file);
      testcase.write_nuH(output_file);
    }
  }
  catch (...) {
    handle_fatal_errors(com);
    return 1;
  }

  return 0;
}
