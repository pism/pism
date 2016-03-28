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
  "\nSSA_TESTCFBC\n"
  "  Testing program for PISM's implementations of the SSA.\n"
  "  Does a time-independent calculation.  Does not run IceModel or a derived\n"
  "  class thereof. Uses the van der Veen flow-line shelf geometry. Also may be\n"
  "  used in a PISM software (regression) test.\n\n";

#include "base/basalstrength/basal_resistance.hh" // IceBasalResistancePlasticLaw
#include "base/stressbalance/ssa/SSAFD.hh"
#include "base/stressbalance/ssa/SSAFD_diagnostics.hh"
#include "base/stressbalance/ssa/SSATestCase.hh"
#include "base/stressbalance/ssa/SSAFEM.hh"
#include "base/util/Mask.hh"
#include "base/util/Context.hh"
#include "base/util/VariableMetadata.hh"
#include "base/util/error_handling.hh"
#include "base/util/iceModelVec.hh"
#include "base/util/io/PIO.hh"
#include "base/util/petscwrappers/PetscInitializer.hh"
#include "base/util/pism_const.hh"
#include "base/util/pism_options.hh"

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

class SSATestCaseCFBC: public SSATestCase {
public:
  SSATestCaseCFBC(Context::Ptr ctx)
    : SSATestCase(ctx) {
    V0 = units::convert(ctx->unit_system(), 300.0, "m year-1", "m second-1");
    H0 = 600.0;                 // meters
    C  = 2.45e-18;
  }

  virtual void write_nuH(const std::string &filename);

protected:
  virtual void initializeGrid(int Mx, int My);

  virtual void initializeSSAModel();

  virtual void initializeSSACoefficients();

  virtual void exactSolution(int i, int j,
    double x, double y, double *u, double *v);

  double V0, //!< grounding line vertically-averaged velocity
    H0,      //!< grounding line thickness (meters)
    C;       //!< "typical constant ice parameter"
};

void SSATestCaseCFBC::write_nuH(const std::string &filename) {

  SSAFD *ssafd = dynamic_cast<SSAFD*>(m_ssa);
  if (ssafd != NULL) {
    SSAFD_nuH(ssafd).compute()->write(filename);
  }
}

void SSATestCaseCFBC::initializeGrid(int Mx, int My) {

  double halfWidth = 250.0e3;  // 500.0 km length
  double Lx = halfWidth, Ly = halfWidth;
  m_grid = IceGrid::Shallow(m_ctx, Lx, Ly,
                            0.0, 0.0, // center: (x0,y0)
                            Mx, My, Y_PERIODIC);
}

void SSATestCaseCFBC::initializeSSAModel() {

  m_config->set_double("ice_softness", pow(1.9e8, -m_config->get_double("ssa_Glen_exponent")));
  m_config->set_boolean("compute_surf_grad_inward_ssa", false);
  m_config->set_boolean("calving_front_stress_boundary_condition", true);
  m_config->set_string("ssa_flow_law", "isothermal_glen");
  m_config->set_string("output_variable_order", "zyx");

  m_enthalpyconverter = EnthalpyConverter::Ptr(new EnthalpyConverter(*m_config));
}

void SSATestCaseCFBC::initializeSSACoefficients() {

  m_tauc.set(0.0);    // irrelevant
  m_bed.set(-1000.0); // assures shelf is floating


  double enth0  = m_enthalpyconverter->enthalpy(273.15, 0.01, 0.0); // 0.01 water fraction
  m_ice_enthalpy.set(enth0);

  IceModelVec::AccessList list;
  list.add(m_thickness);
  list.add(m_surface);
  list.add(m_bc_mask);
  list.add(m_bc_values);
  list.add(m_ice_mask);

  double ocean_rho = m_config->get_double("sea_water_density"),
    ice_rho = m_config->get_double("ice_density");

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    const double x = m_grid->x(i);

    if (i != (int)m_grid->Mx() - 1) {
      m_thickness(i, j) = H_exact(V0, H0, C, x + m_grid->Lx());
      m_ice_mask(i, j)  = MASK_FLOATING;
    } else {
      m_thickness(i, j) = 0;
      m_ice_mask(i, j)  = MASK_ICE_FREE_OCEAN;
    }

    m_surface(i,j) = (1.0 - ice_rho / ocean_rho) * m_thickness(i, j);

    if (i == 0) {
      m_bc_mask(i, j)  = 1;
      m_bc_values(i, j).u = V0;
      m_bc_values(i, j).v = 0;
    } else {
      m_bc_mask(i, j)  = 0;
      m_bc_values(i, j).u = 0;
      m_bc_values(i, j).v = 0;
    }
  }


  // communicate what we have set
  m_surface.update_ghosts();

  m_thickness.update_ghosts();

  m_bc_mask.update_ghosts();

  m_ice_mask.update_ghosts();

  m_bc_values.update_ghosts();

  m_ssa->set_boundary_conditions(m_bc_mask, m_bc_values);
}

void SSATestCaseCFBC::exactSolution(int i, int /*j*/,
                                    double x, double /*y*/,
                                    double *u, double *v) {
  if (i != (int)m_grid->Mx() - 1) {
    *u = u_exact(V0, H0, C, x + m_grid->Lx());
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
  PetscErrorCode ierr;

  com = PETSC_COMM_WORLD;

  /* This explicit scoping forces destructors to be called before PetscFinalize() */
  try {
    verbosityLevelFromOptions();
    Context::Ptr ctx = context_from_options(com, "ssa_test_cfbc");

    bool
      usage_set = options::Bool("-usage", "print usage info"),
      help_set  = options::Bool("-help", "print help info");
    if (usage_set or help_set) {
      ierr = PetscPrintf(com,
                         "\n"
                         "usage of SSA_TEST_CFBC:\n"
                         "  run ssa_test_cfbc -Mx <number> -My <number>\n"
                         "\n");
      PISM_CHK(ierr, "PetscPrintf");
    }

    // Parameters that can be overridden by command line options
    options::Integer Mx("-Mx", "Number of grid points in the X direction", 61);
    options::Integer My("-My", "Number of grid points in the Y direction", 61);

    options::Keyword method("-ssa_method", "Algorithm for computing the SSA solution",
                            "fem,fd", "fd");

    options::String output_file("-o", "Set the output file name", "ssa_test_cfbc.nc");

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

    SSATestCaseCFBC testcase(ctx);
    testcase.init(Mx,My,ssafactory);
    testcase.run();
    testcase.report("V");
    testcase.write(output_file);
    testcase.write_nuH(output_file);
  }
  catch (...) {
    handle_fatal_errors(com);
    return 1;
  }

  return 0;
}
