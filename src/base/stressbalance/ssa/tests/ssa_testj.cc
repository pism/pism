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

static char help[] =
  "\nSSA_TESTJ\n"
  "  Testing program for the finite element implementation of the SSA.\n"
  "  Does a time-independent calculation.  Does not run IceModel or a derived\n"
  "  class thereof. Uses verification test J. Also may be used in a PISM\n"
  "  software (regression) test.\n\n";

#include "base/basalstrength/basal_resistance.hh" // IceBasalResistancePlasticLaw
#include "base/stressbalance/ssa/SSAFD.hh"
#include "base/stressbalance/ssa/SSAFEM.hh"
#include "base/stressbalance/ssa/SSATestCase.hh"
#include "base/util/Mask.hh"
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

class SSATestCaseJ: public SSATestCase
{
public:
  SSATestCaseJ(Context::Ptr ctx)
    : SSATestCase(ctx) {
    // empty
  }

protected:
  virtual void initializeGrid(int Mx,int My);

  virtual void initializeSSAModel();

  virtual void initializeSSACoefficients();

  virtual void exactSolution(int i, int j,
    double x, double y, double *u, double *v);
};

void SSATestCaseJ::initializeGrid(int Mx,int My) {

  double halfWidth = 300.0e3;  // 300.0 km half-width
  double Lx = halfWidth, Ly = halfWidth;
  m_grid = IceGrid::Shallow(m_ctx, Lx, Ly,
                            0.0, 0.0, // center: (x0,y0)
                            Mx, My, XY_PERIODIC);
}

void SSATestCaseJ::initializeSSAModel() {
  m_config->set_boolean("do_pseudo_plastic_till", false);

  m_enthalpyconverter = EnthalpyConverter::Ptr(new EnthalpyConverter(*m_config));
  m_config->set_string("ssa_flow_law", "isothermal_glen");
}

void SSATestCaseJ::initializeSSACoefficients() {
  m_tauc.set(0.0);    // irrelevant for test J
  m_bed.set(-1000.0); // assures shelf is floating (maximum ice thickness is 770 m)
  m_ice_mask.set(MASK_FLOATING);

  double enth0  = m_enthalpyconverter->enthalpy(273.15, 0.01, 0.0); // 0.01 water fraction
  m_ice_enthalpy.set(enth0);

  /* use Ritz et al (2001) value of 30 MPa year for typical vertically-averaged viscosity */
  double ocean_rho = m_config->get_double("sea_water_density"),
    ice_rho = m_config->get_double("ice_density");
  const double nu0 = units::convert(m_sys, 30.0, "MPa year", "Pa s"); /* = 9.45e14 Pa s */
  const double H0 = 500.;       /* 500 m typical thickness */

  // Test J has a viscosity that is independent of velocity.  So we force a
  // constant viscosity by settting the strength_extension
  // thickness larger than the given ice thickness. (max = 770m).
  m_ssa->strength_extension->set_notional_strength(nu0 * H0);
  m_ssa->strength_extension->set_min_thickness(800);

  IceModelVec::AccessList list;
  list.add(m_thickness);
  list.add(m_surface);
  list.add(m_bc_mask);
  list.add(m_bc_values);

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    const double myx = m_grid->x(i), myy = m_grid->y(j);

    // set H,h on regular grid
    struct TestJParameters J_parameters = exactJ(myx, myy);

    m_thickness(i,j) = J_parameters.H;
    m_surface(i,j) = (1.0 - ice_rho / ocean_rho) * J_parameters.H; // FIXME issue #15

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
  m_surface.update_ghosts();
  m_thickness.update_ghosts();
  m_bc_mask.update_ghosts();
  m_bc_values.update_ghosts();

  m_ssa->set_boundary_conditions(m_bc_mask, m_bc_values);
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
  PetscErrorCode ierr;

  com = PETSC_COMM_WORLD;

  /* This explicit scoping forces destructors to be called before PetscFinalize() */
  try {
    verbosityLevelFromOptions();
    Context::Ptr ctx = context_from_options(com, "ssa_testj");
    Config::Ptr config = ctx->config();

    setVerbosityLevel(5);

    bool
      usage_set = options::Bool("-usage", "print usage info"),
      help_set  = options::Bool("-help", "print help info");
    if (usage_set or help_set) {
      ierr = PetscPrintf(com,
                         "\n"
                         "usage of SSA_TESTJ:\n"
                         "  run ssafe_test -Mx <number> -My <number> -ssa_method <fd|fem>\n"
                         "\n");
      PISM_CHK(ierr, "PetscPrintf");
    }

    // Parameters that can be overridden by command line options

    options::Integer Mx("-Mx", "Number of grid points in the X direction", 61);
    options::Integer My("-My", "Number of grid points in the Y direction", 61);

    options::Keyword method("-ssa_method", "Algorithm for computing the SSA solution",
                            "fem,fd", "fem");

    options::String output("-o", "Set the output file name",
                           "ssa_test_j.nc", options::DONT_ALLOW_EMPTY);

    options::Integer verbose("-verbose", "Verbosity level", 2);
    if (verbose.is_set()) {
      setVerbosityLevel(verbose);
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

    SSATestCaseJ testcase(ctx);
    testcase.init(Mx,My,ssafactory);
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
