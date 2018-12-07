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
#include "pism/util/iceModelVec.hh"
#include "pism/util/io/PIO.hh"
#include "pism/util/petscwrappers/PetscInitializer.hh"
#include "pism/util/pism_utilities.hh"
#include "pism/util/pism_options.hh"
#include "pism/verification/tests/exactTestsIJ.h"

namespace pism {
namespace stressbalance {

class SSATestCaseConst: public SSATestCase
{
public:
  SSATestCaseConst(Context::Ptr ctx, int Mx, int My, double q,
                   SSAFactory ssafactory):
    SSATestCase(ctx, Mx, My, 50e3, 50e3, CELL_CORNER, NOT_PERIODIC),
    basal_q(q)
  {
    L     = units::convert(ctx->unit_system(), 50.0, "km", "m"); // 50km half-width
    H0    = 500;                        // m
    dhdx  = 0.005;                      // pure number
    nu0   = units::convert(ctx->unit_system(), 30.0, "MPa year", "Pa s");
    tauc0 = 1.e4;               // Pa

    m_config->set_boolean("basal_resistance.pseudo_plastic.enabled", true);
    m_config->set_double("basal_resistance.pseudo_plastic.q", basal_q);

    // Use a pseudo-plastic law with a constant q determined at run time
    m_config->set_boolean("basal_resistance.pseudo_plastic.enabled", true);

    // The following is irrelevant because we will force linear rheology later.
    m_enthalpyconverter = EnthalpyConverter::Ptr(new EnthalpyConverter(*m_config));

    m_ssa = ssafactory(m_grid);
  };

protected:
  virtual void initializeSSACoefficients();

  virtual void exactSolution(int i, int j,
    double x, double y, double *u, double *v);

  double basal_q,
    L, H0, dhdx, nu0, tauc0;
};

void SSATestCaseConst::initializeSSACoefficients() {

  // Force linear rheology
  m_ssa->strength_extension->set_notional_strength(nu0 * H0);
  m_ssa->strength_extension->set_min_thickness(0.5*H0);

  // The finite difference code uses the following flag to treat the non-periodic grid correctly.
  m_config->set_boolean("stress_balance.ssa.compute_surface_gradient_inward", true);

  // Set constant thickness, tauc
  m_bc_mask.set(0);
  m_geometry.ice_thickness.set(H0);
  m_tauc.set(tauc0);

  IceModelVec::AccessList list{&m_bc_values, &m_bc_mask,
      &m_geometry.bed_elevation, &m_geometry.ice_surface_elevation};

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    double u, v;
    const double x = m_grid->x(i), y=m_grid->y(j);

    m_geometry.bed_elevation(i, j) = -x*(dhdx);
    m_geometry.ice_surface_elevation(i, j) = m_geometry.bed_elevation(i, j) + H0;

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
  double earth_grav = m_config->get_double("constants.standard_gravity"),
    tauc_threshold_velocity = m_config->get_double("basal_resistance.pseudo_plastic.u_threshold",
                                                   "m second-1"),
    ice_rho = m_config->get_double("constants.ice.density");

  *u = pow(ice_rho * earth_grav * H0 * dhdx / tauc0, 1./basal_q)*tauc_threshold_velocity;
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
    Context::Ptr ctx = context_from_options(com, "ssa_test_const");
    Config::Ptr config = ctx->config();

    std::string usage = "\n"
      "usage of SSA_TEST_CONST:\n"
      "  run ssa_test_const -Mx <number> -My <number> -ssa_method <fd|fem>\n"
      "\n";

    bool stop = show_usage_check_req_opts(*ctx->log(), "ssa_test_const", {}, usage);

    if (stop) {
      return 0;
    }

    // Parameters that can be overridden by command line options
    unsigned int Mx = config->get_double("grid.Mx");
    unsigned int My = config->get_double("grid.My");

    config->set_double("basal_resistance.pseudo_plastic.q", 1.0);
    double basal_q = config->get_double("basal_resistance.pseudo_plastic.q");

    auto method = config->get_string("stress_balance.ssa.method");
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

    SSATestCaseConst testcase(ctx, Mx, My, basal_q, ssafactory);
    testcase.init();
    testcase.run();
    testcase.report("const");
    testcase.write(output_file);
  }
  catch (...) {
    handle_fatal_errors(com);
    return 1;
  }

  return 0;
}
