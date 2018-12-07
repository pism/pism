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

static char help[] =
  "\nSSA_TESTJ\n"
  "  Testing program for the finite element implementation of the SSA.\n"
  "  Does a time-independent calculation.  Does not run IceModel or a derived\n"
  "  class thereof. Uses verification test J. Also may be used in a PISM\n"
  "  software (regression) test.\n\n";

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

class SSATestCaseJ: public SSATestCase
{
public:
  SSATestCaseJ(Context::Ptr ctx, int Mx, int My, SSAFactory ssafactory)
    : SSATestCase(ctx, Mx, My, 300e3, 300e3, CELL_CENTER, XY_PERIODIC) {
  m_config->set_boolean("basal_resistance.pseudo_plastic.enabled", false);

  m_enthalpyconverter = EnthalpyConverter::Ptr(new EnthalpyConverter(*m_config));
  m_config->set_string("stress_balance.ssa.flow_law", "isothermal_glen");

  m_ssa = ssafactory(m_grid);
  }

protected:
  virtual void initializeSSACoefficients();

  virtual void exactSolution(int i, int j,
                             double x, double y, double *u, double *v);
};

void SSATestCaseJ::initializeSSACoefficients() {
  m_tauc.set(0.0);    // irrelevant for test J
  m_geometry.bed_elevation.set(-1000.0); // assures shelf is floating (maximum ice thickness is 770 m)
  m_geometry.cell_type.set(MASK_FLOATING);

  double enth0  = m_enthalpyconverter->enthalpy(273.15, 0.01, 0.0); // 0.01 water fraction
  m_ice_enthalpy.set(enth0);

  /* use Ritz et al (2001) value of 30 MPa year for typical vertically-averaged viscosity */
  double ocean_rho = m_config->get_double("constants.sea_water.density"),
    ice_rho = m_config->get_double("constants.ice.density");
  const double nu0 = units::convert(m_sys, 30.0, "MPa year", "Pa s"); /* = 9.45e14 Pa s */
  const double H0 = 500.;       /* 500 m typical thickness */

  // Test J has a viscosity that is independent of velocity.  So we force a
  // constant viscosity by settting the strength_extension
  // thickness larger than the given ice thickness. (max = 770m).
  m_ssa->strength_extension->set_notional_strength(nu0 * H0);
  m_ssa->strength_extension->set_min_thickness(800);

  IceModelVec::AccessList list{&m_geometry.ice_thickness, &m_geometry.ice_surface_elevation, &m_bc_mask, &m_bc_values};

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    const double myx = m_grid->x(i), myy = m_grid->y(j);

    // set H,h on regular grid
    struct TestJParameters J_parameters = exactJ(myx, myy);

    m_geometry.ice_thickness(i,j) = J_parameters.H;
    m_geometry.ice_surface_elevation(i,j) = (1.0 - ice_rho / ocean_rho) * J_parameters.H; // FIXME issue #15

    // special case at center point: here we set bc_values at (i,j) by
    // setting bc_mask and bc_values appropriately
    if ((i == ((int)m_grid->Mx()) / 2) and
        (j == ((int)m_grid->My()) / 2)) {
      m_bc_mask(i,j) = 1;
      m_bc_values(i,j).u = J_parameters.u;
      m_bc_values(i,j).v = J_parameters.v;
    }
  }

  // communicate what we have set
  m_geometry.ice_surface_elevation.update_ghosts();
  m_geometry.ice_thickness.update_ghosts();
  m_bc_mask.update_ghosts();
  m_bc_values.update_ghosts();
}

void SSATestCaseJ::exactSolution(int /*i*/, int /*j*/,
                                 double x, double y,
                                 double *u, double *v) {
  struct TestJParameters J_parameters = exactJ(x, y);
  *u = J_parameters.u;
  *v = J_parameters.v;
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
    Context::Ptr ctx = context_from_options(com, "ssa_testj");
    Config::Ptr config = ctx->config();

    std::string usage = "\n"
      "usage of SSA_TESTJ:\n"
      "  run ssafe_test -Mx <number> -My <number> -ssa_method <fd|fem>\n"
      "\n";

    bool stop = show_usage_check_req_opts(*ctx->log(), "ssa_testj", {}, usage);

    if (stop) {
      return 0;
    }

    // Parameters that can be overridden by command line options

    unsigned int Mx = config->get_double("grid.Mx");
    unsigned int My = config->get_double("grid.My");

    auto method = config->get_string("stress_balance.ssa.method");
    auto output = config->get_string("output.file_name");

    // Determine the kind of solver to use.
    SSAFactory ssafactory = NULL;
    if (method == "fem") {
      ssafactory = SSAFEMFactory;
    } else if (method == "fd") {
      ssafactory = SSAFDFactory;
    } else {
      /* can't happen */
    }

    SSATestCaseJ testcase(ctx, Mx, My, ssafactory);
    testcase.init();
    testcase.run();
    testcase.report("J");
    testcase.write(output);
  }
  catch (...) {
    handle_fatal_errors(com);
    return 1;
  }

  return 0;
}
