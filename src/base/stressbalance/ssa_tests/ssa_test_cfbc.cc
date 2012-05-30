// Copyright (C) 2010--2012 Ed Bueler, Constantine Khroulev, and David Maxwell
//
// This file is part of PISM.
//
// PISM is free software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the Free Software
// Foundation; either version 2 of the License, or (at your option) any later
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

#include "pism_const.hh"
#include "iceModelVec.hh"
#include "flowlaws.hh" // IceFlowLaw
#include "basal_resistance.hh" // IceBasalResistancePlasticLaw
#include "PIO.hh"
#include "NCVariable.hh"
#include "SSAFD.hh"
#include "SSATestCase.hh"
#include "pism_options.hh"

const double H0 = 600;          // meters
const double V0 = 300/secpera;  // 300 meters/year
const double C  = 2.45e-18;     // "typical constant ice parameter"
const double T  = 400;          // time used to compute the calving front location

// thickness profile in the van der Veen solution
static double H_exact(double x)
{
  const double Q0 = V0*H0;
  return pow(4 * C / Q0 * x + 1/pow(H0, 4), -0.25);
}

// theoretical location of the calving front
// static double x_cfbc(double t)
// {
//   const double Q0 = V0*H0;
//   return Q0 / (4*C) * (pow(3*C*t + 1/pow(H0,3),4.0/3.0) - 1/pow(H0,4));
// }

// velocity profile; corresponds to constant flux
static double u_exact(double x)
{
  const double Q0 = V0*H0;
  return Q0 / H_exact(x);
}

class SSATestCaseCFBC: public SSATestCase
{
public:
  SSATestCaseCFBC( MPI_Comm com, PetscMPIInt rank,
                 PetscMPIInt size, NCConfigVariable &c ):
                 SSATestCase(com, rank, size, c)
  { };

protected:
  virtual PetscErrorCode initializeGrid(PetscInt Mx, PetscInt My);

  virtual PetscErrorCode initializeSSAModel();

  virtual PetscErrorCode initializeSSACoefficients();

  virtual PetscErrorCode exactSolution(PetscInt i, PetscInt j,
    PetscReal x, PetscReal y, PetscReal *u, PetscReal *v );

};

PetscErrorCode SSATestCaseCFBC::initializeGrid(PetscInt Mx, PetscInt My)
{

  PetscReal halfWidth = 250.0e3;  // 500.0 km length
  PetscReal Lx = halfWidth, Ly = halfWidth;
  init_shallow_grid(grid, Lx, Ly, Mx, My, Y_PERIODIC);
  return 0;
}

PetscErrorCode SSATestCaseCFBC::initializeSSAModel()
{

  config.set_flag("compute_surf_grad_inward_ssa", true); 
  config.set_flag("calving_front_stress_boundary_condition", true); 
  config.set_string("ssa_flow_law", "isothermal_glen");

  basal = new IceBasalResistancePlasticLaw(
         config.get("plastic_regularization", "1/year", "1/second"),
         config.get_flag("do_pseudo_plastic_till"),
         config.get("pseudo_plastic_q"),
         config.get("pseudo_plastic_uthreshold", "m/year", "m/second"));

  enthalpyconverter = new EnthalpyConverter(config);

  config.set_flag("ssa_flow_law", "isothermal_glen");
  config.set("ice_softness", pow(1.9e8, -config.get("Glen_exponent")));

  return 0;
}

PetscErrorCode SSATestCaseCFBC::initializeSSACoefficients()
{
  PetscErrorCode ierr;

  ierr = tauc.set(0.0); CHKERRQ(ierr);    // irrelevant
  ierr = bed.set(-1000.0); CHKERRQ(ierr); // assures shelf is floating
  ierr = enthalpy.set(528668.35); // arbitrary; corresponds to 263.15 Kelvin at depth=0.
  CHKERRQ(ierr);

  ierr = thickness.begin_access(); CHKERRQ(ierr);
  ierr = surface.begin_access(); CHKERRQ(ierr);
  ierr = bc_mask.begin_access(); CHKERRQ(ierr);
  ierr = vel_bc.begin_access(); CHKERRQ(ierr);
  ierr = ice_mask.begin_access(); CHKERRQ(ierr);

  double ocean_rho = config.get("sea_water_density"),
    ice_rho = config.get("ice_density");

  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      const PetscScalar x = grid.x[i];

      if (x <= 0) {
        thickness(i, j) = H_exact(x + grid.Lx);
        ice_mask(i, j)  = MASK_FLOATING;
      } else {
        thickness(i, j) = 0;
        ice_mask(i, j)  = MASK_ICE_FREE_OCEAN;
      }

      surface(i,j) = (1.0 - ice_rho / ocean_rho) * thickness(i, j); // FIXME issue #15

      if (i == 0) {
        bc_mask(i, j)  = 1;
        vel_bc(i, j).u = V0;
        vel_bc(i, j).v = 0;
      } else {
        bc_mask(i, j)  = 0;
        vel_bc(i, j).u = 0;
        vel_bc(i, j).v = 0;
      }
    }
  }

  ierr = ice_mask.end_access(); CHKERRQ(ierr);
  ierr = surface.end_access(); CHKERRQ(ierr);
  ierr = thickness.end_access(); CHKERRQ(ierr);
  ierr = bc_mask.end_access(); CHKERRQ(ierr);
  ierr = vel_bc.end_access(); CHKERRQ(ierr);

  // communicate what we have set
  ierr = surface.beginGhostComm(); CHKERRQ(ierr);
  ierr = surface.endGhostComm(); CHKERRQ(ierr);

  ierr = thickness.beginGhostComm(); CHKERRQ(ierr);
  ierr = thickness.endGhostComm(); CHKERRQ(ierr);

  ierr = bc_mask.beginGhostComm(); CHKERRQ(ierr);
  ierr = bc_mask.endGhostComm(); CHKERRQ(ierr);

  ierr = vel_bc.beginGhostComm(); CHKERRQ(ierr);
  ierr = vel_bc.endGhostComm(); CHKERRQ(ierr);

  ierr = ssa->set_boundary_conditions(bc_mask, vel_bc); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode SSATestCaseCFBC::exactSolution(PetscInt /*i*/, PetscInt /*j*/,
                                              PetscReal x, PetscReal /*y*/,
                                              PetscReal *u, PetscReal *v)
{
  if (x <= 0) {
    *u = u_exact(x + grid.Lx);
  } else {
    *u = 0;
  }

  *v = 0;
  return 0;
}


int main(int argc, char *argv[]) {
  PetscErrorCode  ierr;

  MPI_Comm    com;  // won't be used except for rank,size
  PetscMPIInt rank, size;

  ierr = PetscInitialize(&argc, &argv, PETSC_NULL, help); CHKERRQ(ierr);

  com = PETSC_COMM_WORLD;
  ierr = MPI_Comm_rank(com, &rank); CHKERRQ(ierr);
  ierr = MPI_Comm_size(com, &size); CHKERRQ(ierr);

  /* This explicit scoping forces destructors to be called before PetscFinalize() */
  {
    NCConfigVariable config, overrides;
    ierr = init_config(com, rank, config, overrides); CHKERRQ(ierr);

    ierr = setVerbosityLevel(5); CHKERRQ(ierr);

    PetscBool usage_set, help_set;
    ierr = PetscOptionsHasName(PETSC_NULL, "-usage", &usage_set); CHKERRQ(ierr);
    ierr = PetscOptionsHasName(PETSC_NULL, "-help", &help_set); CHKERRQ(ierr);
    if ((usage_set==PETSC_TRUE) || (help_set==PETSC_TRUE)) {
      PetscPrintf(com,
                  "\n"
                  "usage of SSA_TEST_CFBC:\n"
                  "  run ssa_test_cfbc -Mx <number> -My <number>\n"
                  "\n");
    }

    // Parameters that can be overridden by command line options
    PetscInt Mx=61;
    PetscInt My=61;
    string output_file = "ssa_test_cfbc.nc";

    ierr = PetscOptionsBegin(com, "", "SSA_TESTCFBC options", ""); CHKERRQ(ierr);
    {
      bool flag;
      int my_verbosity_level;
      ierr = PISMOptionsInt("-Mx", "Number of grid points in the X direction",
                                                      Mx, flag); CHKERRQ(ierr);
      ierr = PISMOptionsInt("-My", "Number of grid points in the Y direction",
                                                      My, flag); CHKERRQ(ierr);
      ierr = PISMOptionsString("-o", "Set the output file name",
                                              output_file, flag); CHKERRQ(ierr);
      ierr = PISMOptionsInt("-verbose", "Verbosity level",
                            my_verbosity_level, flag); CHKERRQ(ierr);
      if (flag) setVerbosityLevel(my_verbosity_level);
    }
    ierr = PetscOptionsEnd(); CHKERRQ(ierr);

    SSAFactory ssafactory = SSAFDFactory;

    SSATestCaseCFBC testcase(com,rank,size,config);
    ierr = testcase.init(Mx,My,ssafactory); CHKERRQ(ierr);
    ierr = testcase.run(); CHKERRQ(ierr);
    ierr = testcase.report("CFBC"); CHKERRQ(ierr);
    ierr = testcase.write(output_file);
  }

  ierr = PetscFinalize(); CHKERRQ(ierr);
  return 0;
}
