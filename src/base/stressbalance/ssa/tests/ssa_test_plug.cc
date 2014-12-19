// Copyright (C) 2010--2014 Ed Bueler, Constantine Khroulev, and David Maxwell
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

#include "pism_const.hh"
#include "pism_options.hh"
#include "iceModelVec.hh"
#include "flowlaws.hh" // IceFlowLaw
#include "basal_resistance.hh" // IceBasalResistancePlasticLaw
#include "PIO.hh"
#include "NCVariable.hh"
#include "SSAFEM.hh"
#include "SSAFD.hh"
#include "exactTestsIJ.h"
#include "SSATestCase.hh"
#include <math.h>

#include "PetscInitializer.hh"
#include "error_handling.hh"

using namespace pism;

class SSATestCasePlug: public SSATestCase
{
public:
  SSATestCasePlug(MPI_Comm com, Config &c, double n)
    : SSATestCase(com, c)
  {
    H0 = 2000.; //m
    L=50.e3; // 50km half-width
    dhdx = 0.001; // pure number, slope of surface & bed
    tauc0 = 0.; // No basal shear stress
    B0 = 3.7e8; // Pa s^{1/3}; hardness
               // given on p. 239 of Schoof; why so big?
    this->glen_n = n;
  }

protected:
  virtual PetscErrorCode initializeGrid(int Mx,int My);

  virtual PetscErrorCode initializeSSAModel();

  virtual PetscErrorCode initializeSSACoefficients();

  virtual PetscErrorCode exactSolution(int i, int j,
    double x, double y, double *u, double *v);


  double H0; // Thickness
  double L;  // Half-width
  double dhdx; // surface slope
  double tauc0; // zero basal shear stress
  double B0;  // hardness
  double glen_n;

  bool dimensionless;

};


PetscErrorCode SSATestCasePlug::initializeGrid(int Mx,int My)
{
  double Lx=L, Ly = L;
  grid = IceGrid::Shallow(m_com, config, Lx, Ly,
                          0.0, 0.0, // center: (x0,y0)
                          Mx, My, NOT_PERIODIC);
  return 0;
}


PetscErrorCode SSATestCasePlug::initializeSSAModel()
{
  // Basal sliding law parameters are irrelevant because tauc=0

  // Enthalpy converter is irrelevant (but still required) for this test.
  enthalpyconverter = new EnthalpyConverter(config);

  // Use constant hardness
  config.set_string("ssa_flow_law", "isothermal_glen");
  config.set_double("ice_softness", pow(B0, -glen_n));
  return 0;
}

PetscErrorCode SSATestCasePlug::initializeSSACoefficients()
{

  // The finite difference code uses the following flag to treat the non-periodic grid correctly.
  config.set_flag("compute_surf_grad_inward_ssa", true);
  config.set_double("epsilon_ssa", 0.0);

  // Ensure we never use the strength extension.
  ssa->strength_extension->set_min_thickness(H0/2);

  // Set constant coefficients.
  thickness.set(H0);
  tauc.set(tauc0);


  // Set boundary conditions (Dirichlet all the way around).
  bc_mask.set(0.0);

  IceModelVec::AccessList list;
  list.add(vel_bc);
  list.add(bc_mask);
  list.add(bed);
  list.add(surface);

  for (Points p(*grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    double myu, myv;
    const double myx = grid->x(i), myy=grid->y(j);

    bed(i,j) = -myx*(dhdx);
    surface(i,j) = bed(i,j) + H0;

    bool edge = ((j == 0) || (j == (int)grid->My() - 1) ||
                 (i == 0) || (i == (int)grid->Mx() - 1));
    if (edge) {
      bc_mask(i,j) = 1;
      exactSolution(i,j,myx,myy,&myu,&myv);
      vel_bc(i,j).u = myu;
      vel_bc(i,j).v = myv;
    }
  }

  vel_bc.update_ghosts();
  bc_mask.update_ghosts();
  bed.update_ghosts();
  surface.update_ghosts();

  ssa->set_boundary_conditions(bc_mask, vel_bc);

  return 0;
}

PetscErrorCode SSATestCasePlug::exactSolution(int /*i*/, int /*j*/,
                                              double /*x*/, double y,
                                              double *u, double *v)
{
  double earth_grav = config.get("standard_gravity"),
    ice_rho = config.get("ice_density");
  double f = ice_rho * earth_grav * H0* dhdx;
  double ynd = y/L;

  *u = 0.5*pow(f,3)*pow(L,4)/pow(B0*H0,3)*(1-pow(ynd,4));
  *v = 0;
  return 0;
}

int main(int argc, char *argv[]) {

  MPI_Comm com = MPI_COMM_WORLD;  // won't be used except for rank,size

  PetscInitializer petsc(argc, argv, help);

  com = PETSC_COMM_WORLD;

  /* This explicit scoping forces destructors to be called before PetscFinalize() */
  try {
    UnitSystem unit_system;
    Config config(com, "pism_config", unit_system),
      overrides(com, "pism_overrides", unit_system);
    init_config(com, config, overrides);

    setVerbosityLevel(5);

    bool
      usage_set = options::Bool("-usage", "show the usage message"),
      help_set  = options::Bool("-help", "show the help message");
    if (usage_set or help_set) {
      PetscPrintf(com,
                  "\n"
                  "usage of SSA_TEST_PLUG:\n"
                  "  run ssa_test_plug -Mx <number> -My <number> -ssa_method <fd|fem>\n"
                  "\n");
    }

    // Parameters that can be overridden by command line options

    options::Integer Mx("-Mx", "Number of grid points in the X direction", 11);
    options::Integer My("-My", "Number of grid points in the Y direction", 61);

    options::Keyword method("-ssa_method", "Algorithm for computing the SSA solution",
                            "fem,fd", "fem");

    options::String output_file("-o", "Set the output file name", "ssa_test_plug.nc");
    options::Real glen_n("-ssa_glen_n", "Glen exponent for the SSA", 3.0);

    options::Integer my_verbosity_level("-verbose", "Verbosity level", 2);
    if (my_verbosity_level.is_set()) {
      setVerbosityLevel(my_verbosity_level);
    }

    // Determine the kind of solver to use.
    SSAFactory ssafactory = NULL;
    if (method.value() == "fem") {
      ssafactory = SSAFEMFactory;
    } else if (method.value() == "fd") {
      ssafactory = SSAFDFactory;
    } else {
      /* can't happen */
    }

    SSATestCasePlug testcase(com,config,glen_n);
    testcase.init(Mx,My,ssafactory);
    testcase.run();
    testcase.report("plug");
    testcase.write(output_file);
  }
  catch (...) {
    handle_fatal_errors(com);
  }

  return 0;
}
