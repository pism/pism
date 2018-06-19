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

class SSATestCasePlug: public SSATestCase {
public:
  SSATestCasePlug(Context::Ptr ctx, int Mx, int My,
                  double n, SSAFactory ssafactory)
    : SSATestCase(ctx, Mx, My, 50e3, 50e3, CELL_CORNER, NOT_PERIODIC) {
    H0    = 2000.;              // m
    L     = 50.e3;              // 50km half-width
    dhdx  = 0.001;              // pure number, slope of surface & bed
    tauc0 = 0.;                 // No basal shear stress

    // Pa s^{1/3}; hardness given on p. 239 of Schoof; why so big?
    B0           = 3.7e8;

    this->glen_n = n;

    // Basal sliding law parameters are irrelevant because tauc=0

    // Enthalpy converter is irrelevant (but still required) for this test.
    m_enthalpyconverter = EnthalpyConverter::Ptr(new EnthalpyConverter(*m_config));

    // Use constant hardness
    m_config->set_string("stress_balance.ssa.flow_law", "isothermal_glen");
    m_config->set_double("flow_law.isothermal_Glen.ice_softness", pow(B0, -glen_n));

    m_ssa = ssafactory(m_grid);
  }

protected:
  virtual void initializeSSACoefficients();

  virtual void exactSolution(int i, int j,
    double x, double y, double *u, double *v);


  double H0; // Thickness
  double L;  // Half-width
  double dhdx; // surface slope
  double tauc0; // zero basal shear stress
  double B0;  // hardness
  double glen_n;

  bool dimensionless;

};

void SSATestCasePlug::initializeSSACoefficients() {

  // The finite difference code uses the following flag to treat the non-periodic grid correctly.
  m_config->set_boolean("stress_balance.ssa.compute_surface_gradient_inward", true);
  m_config->set_double("stress_balance.ssa.epsilon", 0.0);

  // Ensure we never use the strength extension.
  m_ssa->strength_extension->set_min_thickness(H0/2);

  // Set constant coefficients.
  m_geometry.ice_thickness.set(H0);
  m_tauc.set(tauc0);

  // Set boundary conditions (Dirichlet all the way around).
  m_bc_mask.set(0.0);

  IceModelVec::AccessList list{&m_bc_values, &m_bc_mask, &m_geometry.bed_elevation, &m_geometry.ice_surface_elevation};

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    double myu, myv;
    const double myx = m_grid->x(i), myy=m_grid->y(j);

    m_geometry.bed_elevation(i,j) = -myx*(dhdx);
    m_geometry.ice_surface_elevation(i,j) = m_geometry.bed_elevation(i,j) + H0;

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
  m_geometry.bed_elevation.update_ghosts();
  m_geometry.ice_surface_elevation.update_ghosts();
}

void SSATestCasePlug::exactSolution(int /*i*/, int /*j*/,
                                    double /*x*/, double y,
                                    double *u, double *v) {
  double earth_grav = m_config->get_double("constants.standard_gravity"),
    ice_rho = m_config->get_double("constants.ice.density");
  double f = ice_rho * earth_grav * H0* dhdx;
  double ynd = y/L;

  *u = 0.5*pow(f,3)*pow(L,4)/pow(B0*H0,3)*(1-pow(ynd,4));
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
    Context::Ptr ctx = context_from_options(com, "ssa_test_plug");
    Config::Ptr config = ctx->config();

    std::string usage = "\n"
      "usage of SSA_TEST_PLUG:\n"
      "  run ssa_test_plug -Mx <number> -My <number> -ssa_method <fd|fem>\n"
      "\n";

    bool stop = show_usage_check_req_opts(*ctx->log(), "ssa_test_plug", {}, usage);

    if (stop) {
      return 0;
    }

    // Parameters that can be overridden by command line options

    unsigned int Mx = config->get_double("grid.Mx");
    unsigned int My = config->get_double("grid.My");

    auto method      = config->get_string("stress_balance.ssa.method");
    auto output_file = config->get_string("output.file_name");
    auto glen_n      = config->get_double("stress_balance.ssa.Glen_exponent");

    // Determine the kind of solver to use.
    SSAFactory ssafactory = NULL;
    if (method == "fem") {
      ssafactory = SSAFEMFactory;
    } else if (method == "fd") {
      ssafactory = SSAFDFactory;
    } else {
      /* can't happen */
    }

    SSATestCasePlug testcase(ctx, Mx, My, glen_n, ssafactory);
    testcase.init();
    testcase.run();
    testcase.report("plug");
    testcase.write(output_file);
  }
  catch (...) {
    handle_fatal_errors(com);
    return 1;
  }

  return 0;
}
