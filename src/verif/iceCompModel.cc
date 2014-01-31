// Copyright (C) 2004-2014 Jed Brown, Ed Bueler and Constantine Khroulev
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

#include <cmath>
#include <cstring>
#include <petscdmda.h>
#include <assert.h>

#include <vector>     // STL vector container; sortable; used in test L
#include <algorithm>  // required by sort(...) in test L

#include "tests/exactTestsABCDE.h"
#include "tests/exactTestsFG.h"
#include "tests/exactTestH.h"
#include "tests/exactTestL.h"

#include "iceCompModel.hh"
#include "SIA_Sliding.hh"
#include "SIAFD.hh"
#include "flowlaw_factory.hh"
#include "PISMStressBalance.hh"
#include "enthalpyConverter.hh"
#include "PIO.hh"
#include "pism_options.hh"
#include "POConstant.hh"
#include "PSVerification.hh"

const double IceCompModel::secpera = 3.15569259747e7;

IceCompModel::IceCompModel(IceGrid &g, PISMConfig &conf, PISMConfig &conf_overrides, int mytest)
  : IceModel(g, conf, conf_overrides) {

  // note lots of defaults are set by the IceModel constructor

  // defaults for IceCompModel:
  testname = mytest;
  exactOnly = PETSC_FALSE;
  bedrock_is_ice_forK = PETSC_FALSE;

  // Override some defaults from parent class
  config.set_double("sia_enhancement_factor", 1.0);
  // none use bed smoothing & bed roughness parameterization
  config.set_double("bed_smoother_range", 0.0);

  // set values of flags in run()
  config.set_flag("do_mass_conserve", true);
  config.set_flag("include_bmr_in_continuity", false);

  if (testname == 'V') {
    config.set_string("ssa_flow_law", "isothermal_glen");
    config.set_double("ice_softness", pow(1.9e8, -config.get("Glen_exponent")));
  } else {
    // Set the default for IceCompModel:
    config.set_string("sia_flow_law", "arr");
  }
}

PetscErrorCode IceCompModel::createVecs() {
  PetscErrorCode ierr;

  ierr = IceModel::createVecs(); CHKERRQ(ierr);

  ierr = vHexactL.create(grid, "HexactL", WITH_GHOSTS, 2); CHKERRQ(ierr);

  ierr = strain_heating3_comp.create(grid,"strain_heating_comp", WITHOUT_GHOSTS); CHKERRQ(ierr);
  ierr = strain_heating3_comp.set_attrs("internal","rate of compensatory strain heating in ice",
			      "W m-3", ""); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode IceCompModel::set_grid_defaults() {
  PetscErrorCode ierr;

  // This sets the defaults for each test; command-line options can override this.

  // equal spacing is the default for all the tests except K
  grid.ice_vertical_spacing = EQUAL;

  switch (testname) {
  case 'A':
  case 'E':
    // use 1600km by 1600km by 4000m rectangular domain
    grid.Lx = grid.Ly = 800e3;
    grid.Lz = 4000;
    break;
  case 'B':
  case 'H':
    // use 2400km by 2400km by 4000m rectangular domain
    grid.Lx = grid.Ly = 1200e3;
    grid.Lz = 4000;
    break;
  case 'C':
  case 'D':
    // use 2000km by 2000km by 4000m rectangular domain
    grid.Lx = grid.Ly = 1000e3;
    grid.Lz = 4000;
    break;
  case 'F':
  case 'G':
  case 'L':
    // use 1800km by 1800km by 4000m rectangular domain
    grid.Lx = grid.Ly = 900e3;
    grid.Lz = 4000;
    break;
  case 'K':
  case 'O':
    // use 2000km by 2000km by 4000m rectangular domain, but make truely periodic
    config.set_double("grid_Mbz", 2);
    config.set_double("grid_Lbz", 1000);
    grid.Lx = grid.Ly = 1000e3;
    grid.Lz = 4000;
    grid.periodicity = XY_PERIODIC;
    grid.ice_vertical_spacing = QUADRATIC;
    break;
  case 'V':
    grid.My = 3;                // it's a flow-line setup
    grid.Lx = 500e3;            // 500 km long
    grid.periodicity = Y_PERIODIC;
    break;
  default:
    ierr = PetscPrintf(grid.com, "IceCompModel ERROR : desired test not implemented\n");
    CHKERRQ(ierr);
    PISMEnd();
  }

  ierr =  grid.time->init(); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode IceCompModel::setFromOptions() {
  PetscErrorCode ierr;

  ierr = verbPrintf(2, grid.com, "starting Test %c ...\n", testname);  CHKERRQ(ierr);

  /* This switch turns off actual numerical evolution and simply reports the
     exact solution. */
  bool flag;
  ierr = PISMOptionsIsSet("-eo", flag); CHKERRQ(ierr);
  if (flag) {
    exactOnly = PETSC_TRUE;
    ierr = verbPrintf(1,grid.com, "!!EXACT SOLUTION ONLY, NO NUMERICAL SOLUTION!!\n");
    CHKERRQ(ierr);
  }

  // These ifs are here (and not in the constructor or later) because
  // testname actually comes from a command-line *and* because command-line
  // options should be able to override parameter values set here.

  if (testname == 'H') {
    config.set_string("bed_deformation_model", "iso");
  } else
    config.set_string("bed_deformation_model", "none");

  if ((testname == 'F') || (testname == 'G') || (testname == 'K') || (testname == 'O')) {
    config.set_flag("do_energy", true);
    // essentially turn off run-time reporting of extremely low computed
    // temperatures; *they will be reported as errors* anyway
    config.set_double("global_min_allowed_temp", 0.0);
    config.set_double("max_low_temp_count", 1000000);
  } else
    config.set_flag("do_energy", false);

  config.set_flag("is_dry_simulation", true);

  // special considerations for K and O wrt thermal bedrock and pressure-melting
  if ((testname == 'K') || (testname == 'O')) {
    config.set_flag("temperature_allow_above_melting", false);
  } else {
    // note temps are generally allowed to go above pressure melting in verify
    config.set_flag("temperature_allow_above_melting", true);
  }

  if (testname == 'V') {
    // no sub-shelf melting
    config.set_flag("include_bmr_in_continuity", false);

    // this test is isothermal
    config.set_flag("do_energy", false);

    // do not use the SIA stress balance
    config.set_flag("do_sia", false);

    // do use the SSA solver
    config.set_string("stress_balance_model", "ssa");

    // this certainly is not a "dry simulation"
    config.set_flag("is_dry_simulation", false);

    config.set_flag("ssa_dirichlet_bc", true);
  }

  config.set_flag("do_cold_ice_methods", true);

  ierr = IceModel::setFromOptions();CHKERRQ(ierr);

  return 0;
}

PetscErrorCode IceCompModel::allocate_enthalpy_converter() {

  if (EC != NULL)
    return 0;

  // allocate the "special" enthalpy converter;
  EC = new ICMEnthalpyConverter(config);

  return 0;
}

PetscErrorCode IceCompModel::allocate_bedrock_thermal_unit() {
  PetscErrorCode ierr;

  if (btu != NULL)
    return 0;

  // this switch changes Test K to make material properties for bedrock the same as for ice
  bool biiSet;
  ierr = PISMOptionsIsSet("-bedrock_is_ice", biiSet); CHKERRQ(ierr);
  if (biiSet == PETSC_TRUE) {
    if (testname == 'K') {
      ierr = verbPrintf(1,grid.com,
                        "setting material properties of bedrock to those of ice in Test K\n"); CHKERRQ(ierr);
      config.set_double("bedrock_thermal_density", config.get("ice_density"));
      config.set_double("bedrock_thermal_conductivity", config.get("ice_thermal_conductivity"));
      config.set_double("bedrock_thermal_specific_heat_capacity", config.get("ice_specific_heat_capacity"));
      bedrock_is_ice_forK = PETSC_TRUE;
    } else {
      ierr = verbPrintf(1,grid.com,
                        "IceCompModel WARNING: option -bedrock_is_ice ignored; only applies to Test K\n");
      CHKERRQ(ierr);
    }
  }

  if (testname != 'K') {
    // now make bedrock have same material properties as ice
    // (note Mbz=1 also, by default, but want ice/rock interface to see
    // pure ice from the point of view of applying geothermal boundary
    // condition, especially in tests F and G)
    config.set_double("bedrock_thermal_density", config.get("ice_density"));
    config.set_double("bedrock_thermal_conductivity", config.get("ice_thermal_conductivity"));
    config.set_double("bedrock_thermal_specific_heat_capacity", config.get("ice_specific_heat_capacity"));
  }

  btu = new BTU_Verification(grid, config, testname, bedrock_is_ice_forK);

  return 0;
}

PetscErrorCode IceCompModel::allocate_stressbalance() {
  PetscErrorCode ierr;

  if (stress_balance != NULL)
    return 0;

  if (testname == 'E') {
    config.set_flag("sia_sliding_verification_mode", true);
    ShallowStressBalance *ssb = new SIA_Sliding(grid, *EC, config);
    SIAFD *sia = new SIAFD(grid, *EC, config);

    stress_balance = new PISMStressBalance(grid, ssb, sia, config);
    ierr = stress_balance->init(variables); CHKERRQ(ierr);
  } else {
    ierr = IceModel::allocate_stressbalance(); CHKERRQ(ierr);
  }

  if (testname != 'V') {
    // check on whether the options (already checked) chose the right
    // IceFlowLaw for verification (we need to have the right flow law for
    // errors to make sense)

    IceFlowLaw *ice = stress_balance->get_ssb_modifier()->get_flow_law();

    if (IceFlowLawIsPatersonBuddCold(ice, config, EC) == PETSC_FALSE) {
      ierr = verbPrintf(1, grid.com,
                        "WARNING: SIA flow law should be '-sia_flow_law arr' for the selected pismv test.\n");
      CHKERRQ(ierr);
    }
  }

  return 0;
}

PetscErrorCode IceCompModel::allocate_bed_deformation() {
  PetscErrorCode ierr;

  ierr = IceModel::allocate_bed_deformation(); CHKERRQ(ierr);

  f = config.get("ice_density") / config.get("lithosphere_density");  // for simple isostasy

  std::string bed_def_model = config.get_string("bed_deformation_model");

  if ( (testname == 'H') && bed_def_model != "iso" ) {
    ierr = verbPrintf(1,grid.com,
                      "IceCompModel WARNING: Test H should be run with option\n"
                      "  '-bed_def iso'  for the reported errors to be correct.\n"); CHKERRQ(ierr);
  }

  return 0;
}

PetscErrorCode IceCompModel::allocate_couplers() {
  // Climate will always come from verification test formulas.
  surface = new PSVerification(grid, config, EC, testname);
  ocean   = new POConstant(grid, config);

  return 0;
}

PetscErrorCode IceCompModel::set_vars_from_options() {
  PetscErrorCode ierr;

  // -boot_file command-line option is not allowed here.
  ierr = stop_if_set(grid.com, "-boot_file"); CHKERRQ(ierr);

  ierr = strain_heating3_comp.set(0.0); CHKERRQ(ierr);

  ierr = verbPrintf(3,grid.com,
                    "initializing Test %c from formulas ...\n",testname);  CHKERRQ(ierr);

  // all have no uplift
  ierr = bed_uplift_rate.set(0.0); CHKERRQ(ierr);

  // this is the correct initialization for Test O (and every other
  // test; they all generate zero basal melt rate)
  ierr = basal_melt_rate.set(0.0); CHKERRQ(ierr);

  // Test-specific initialization:
  switch (testname) {
  case 'A':
  case 'B':
  case 'C':
  case 'D':
  case 'E':
  case 'H':
    ierr = initTestABCDEH(); CHKERRQ(ierr);
    break;
  case 'F':
  case 'G':
    ierr = initTestFG(); CHKERRQ(ierr);  // see iCMthermo.cc
    break;
  case 'K':
  case 'O':
    ierr = initTestsKO(); CHKERRQ(ierr);  // see iCMthermo.cc
    break;
  case 'L':
    ierr = initTestL(); CHKERRQ(ierr);
    break;
  case 'V':
    ierr = test_V_init(); CHKERRQ(ierr);
    break;
  default:
    SETERRQ(grid.com, 1,"Desired test not implemented by IceCompModel.\n");
  }

  ierr = compute_enthalpy_cold(T3, Enth3); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode IceCompModel::initTestABCDEH() {
  PetscErrorCode  ierr;
  PetscScalar     A0, T0, H, accum, dummy1, dummy2, dummy3;

  ThermoGlenArrIce tgaIce(grid.com, "sia_", config, EC);

  const double time = grid.time->current();

  // compute T so that A0 = A(T) = Acold exp(-Qcold/(R T))  (i.e. for ThermoGlenArrIce);
  // set all temps to this constant
  A0 = 1.0e-16/secpera;    // = 3.17e-24  1/(Pa^3 s);  (EISMINT value) flow law parameter
  T0 = tgaIce.tempFromSoftness(A0);

  ierr = T3.set(T0);               CHKERRQ(ierr);
  ierr = geothermal_flux.set(Ggeo);           CHKERRQ(ierr);
  ierr = vMask.set(MASK_GROUNDED); CHKERRQ(ierr);

  ierr = ice_thickness.begin_access(); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      PetscScalar xx = grid.x[i], yy = grid.y[j],
        r = grid.radius(i, j);
      switch (testname) {
        case 'A':
          exactA(r, &H, &accum);
          ice_thickness(i, j)   = H;
          break;
        case 'B':
          exactB(time, r, &H, &accum);
          ice_thickness(i, j)   = H;
          break;
        case 'C':
          exactC(time, r, &H, &accum);
          ice_thickness(i, j)   = H;
          break;
        case 'D':
          exactD(time, r, &H, &accum);
          ice_thickness(i, j)   = H;
          break;
        case 'E':
          exactE(xx, yy, &H, &accum, &dummy1, &dummy2, &dummy3);
          ice_thickness(i, j)   = H;
          break;
        case 'H':
          exactH(f, time, r, &H, &accum);
          ice_thickness(i, j)   = H;
          break;
        default:  SETERRQ(grid.com, 1, "test must be A, B, C, D, E, or H");
      }
    }
  }
  ierr = ice_thickness.end_access(); CHKERRQ(ierr);

  ierr = ice_thickness.update_ghosts(); CHKERRQ(ierr);

  if (testname == 'H') {
    ierr = ice_thickness.copy_to(bed_topography); CHKERRQ(ierr);
    ierr = bed_topography.scale(-f); CHKERRQ(ierr);
  } else {  // flat bed case otherwise
    ierr = bed_topography.set(0.0); CHKERRQ(ierr);
  }

  return 0;
}

//! Class used initTestL() in generating sorted list for ODE solver.
class rgrid {
public:
  double r;
  int    i,j;
};

//! Comparison used initTestL() in generating sorted list for ODE solver.
struct rgridReverseSort {
  bool operator()(rgrid a, rgrid b) { return (a.r > b.r); }
};

PetscErrorCode IceCompModel::initTestL() {
  PetscErrorCode  ierr;
  PetscScalar     A0, T0;

  assert(testname == 'L');

  ThermoGlenArrIce tgaIce(grid.com, "sia_", config, EC);

  // compute T so that A0 = A(T) = Acold exp(-Qcold/(R T))  (i.e. for ThermoGlenArrIce);
  // set all temps to this constant
  A0 = 1.0e-16/secpera;    // = 3.17e-24  1/(Pa^3 s);  (EISMINT value) flow law parameter
  T0 = tgaIce.tempFromSoftness(A0);

  ierr = T3.set(T0);                CHKERRQ(ierr); 
  ierr = geothermal_flux.set(Ggeo); CHKERRQ(ierr); 

  // setup to evaluate test L; requires solving an ODE numerically
  //   using sorted list of radii, sorted in decreasing radius order
  const int MM = grid.xm * grid.ym;

  std::vector<rgrid> rrv(MM);  // destructor at end of scope
  int k = 0;
  for (PetscInt   i = grid.xs; i < grid.xs+grid.xm; ++i) {
    for (PetscInt j = grid.ys; j < grid.ys+grid.ym; ++j) {
      rrv[k].i = i;
      rrv[k].j = j;
      rrv[k].r = grid.radius(i,j);

      k += 1;
    }
  }

  std::sort(rrv.begin(), rrv.end(), rgridReverseSort()); // so rrv[k].r > rrv[k+1].r

  // get soln to test L at these radii; solves ODE only once (on each processor)
  std::vector<double> rr(MM), HH(MM), bb(MM), aa(MM);

  for (k = 0; k < MM; k++)
    rr[k] = rrv[k].r;

  ierr = exactL_list(&rr[0], MM, &HH[0], &bb[0], &aa[0]);
  switch (ierr) {
     case TESTL_NOT_DONE:
       verbPrintf(1,grid.com,
          "\n\nTest L ERROR: exactL_list() returns 'NOT_DONE' ...\n\n\n",ierr);
       break;
     case TESTL_NOT_DECREASING:
       verbPrintf(1,grid.com,
          "\n\nTest L ERROR: exactL_list() returns 'NOT_DECREASING' ...\n\n\n",ierr);
       break;
     case TESTL_INVALID_METHOD:
       verbPrintf(1,grid.com,
          "\n\nTest L ERROR: exactL_list() returns 'INVALID_METHOD' ...\n\n\n",ierr);
       break;
     case TESTL_NO_LIST:
       verbPrintf(1,grid.com,
          "\n\nTest L ERROR: exactL_list() returns 'NO_LIST' ...\n\n\n",ierr);
       break;
     default:
       break;
  }
  CHKERRQ(ierr);

  ierr = ice_thickness.begin_access();  CHKERRQ(ierr); 
  ierr = bed_topography.begin_access(); CHKERRQ(ierr); 
  for (k = 0; k < MM; k++) {
    ice_thickness(rrv[k].i, rrv[k].j)  = HH[k];
    bed_topography(rrv[k].i, rrv[k].j) = bb[k];
  }
  ierr = ice_thickness.end_access();  CHKERRQ(ierr); 
  ierr = bed_topography.end_access(); CHKERRQ(ierr); 

  ierr = ice_thickness.update_ghosts();  CHKERRQ(ierr); 
  ierr = bed_topography.update_ghosts(); CHKERRQ(ierr); 

  // store copy of ice_thickness for "-eo" runs and for evaluating geometry errors
  ierr = ice_thickness.copy_to(vHexactL); CHKERRQ(ierr);

  return 0;
}

//! \brief Tests A and E have a thickness B.C. (ice_thickness == 0 outside a circle of radius 750km).
PetscErrorCode IceCompModel::reset_thickness_tests_AE() {
  PetscErrorCode ierr;
  const PetscScalar LforAE = 750e3; // m

  ierr = ice_thickness.begin_access(); CHKERRQ(ierr);
  for (PetscInt   i = grid.xs; i < grid.xs+grid.xm; ++i) {
    for (PetscInt j = grid.ys; j < grid.ys+grid.ym; ++j) {
      if (grid.radius(i, j) > LforAE)
        ice_thickness(i, j) = 0;
    }
  }
  ierr = ice_thickness.end_access(); CHKERRQ(ierr);

  ierr = ice_thickness.update_ghosts(); CHKERRQ(ierr);
  return 0;
}



PetscErrorCode IceCompModel::fillSolnTestABCDH() {
  PetscErrorCode  ierr;
  PetscScalar     H, accum;

  const double time = grid.time->current();

  ierr = ice_thickness.begin_access(); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      PetscScalar r = grid.radius(i, j);
      switch (testname) {
        case 'A':
          exactA(r, &H, &accum);
          ice_thickness(i, j)   = H;
          break;
        case 'B':
          exactB(time, r, &H, &accum);
          ice_thickness(i, j)   = H;
          break;
        case 'C':
          exactC(time, r, &H, &accum);
          ice_thickness(i, j)   = H;
          break;
        case 'D':
          exactD(time, r, &H, &accum);
          ice_thickness(i, j)   = H;
          break;
        case 'H':
          exactH(f, time, r, &H, &accum);
          ice_thickness(i, j)   = H;
          break;
        default:  SETERRQ(grid.com, 1, "test must be A, B, C, D, or H");
      }
    }
  }
  ierr = ice_thickness.end_access(); CHKERRQ(ierr);

  ierr = ice_thickness.update_ghosts(); CHKERRQ(ierr);

  if (testname == 'H') {
    ierr = ice_thickness.copy_to(bed_topography); CHKERRQ(ierr);
    ierr = bed_topography.scale(-f); CHKERRQ(ierr);
    ierr = bed_topography.update_ghosts(); CHKERRQ(ierr);
  }
  return 0;
}


PetscErrorCode IceCompModel::fillSolnTestE() {
  PetscErrorCode  ierr;
  PetscScalar     H, accum, dummy;
  PISMVector2     bvel;
  IceModelVec2V *vel_adv;
  ierr = stress_balance->get_2D_advective_velocity(vel_adv); CHKERRQ(ierr);

  ierr = ice_thickness.begin_access(); CHKERRQ(ierr); 
  ierr = vel_adv->begin_access();      CHKERRQ(ierr); 
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      PetscScalar xx = grid.x[i], yy = grid.y[j];
      exactE(xx, yy, &H, &accum, &dummy, &bvel.u, &bvel.v);
      ice_thickness(i,j) = H;
      (*vel_adv)(i,j)    = bvel;
    }
  }
  ierr = ice_thickness.end_access(); CHKERRQ(ierr); 
  ierr = vel_adv->end_access();      CHKERRQ(ierr); 

  ierr = ice_thickness.update_ghosts(); CHKERRQ(ierr); 
  return 0;
}


PetscErrorCode IceCompModel::fillSolnTestL() {
  PetscErrorCode  ierr;

  ierr = vHexactL.update_ghosts(); CHKERRQ(ierr);
  ierr = ice_thickness.copy_from(vHexactL);

  // note bed was filled at initialization and hasn't changed
  return 0;
}


PetscErrorCode IceCompModel::computeGeometryErrors(PetscScalar &gvolexact, PetscScalar &gareaexact,
                                                   PetscScalar &gdomeHexact, PetscScalar &volerr,
                                                   PetscScalar &areaerr, PetscScalar &gmaxHerr,
                                                   PetscScalar &gavHerr, PetscScalar &gmaxetaerr,
                                                   PetscScalar &centerHerr) {
  // compute errors in thickness, eta=thickness^{(2n+2)/n}, volume, area

  const double time = grid.time->current();

  PetscErrorCode ierr;
  double
    Hexact     = 0.0,
    vol        = 0.0,
    area       = 0.0,
    domeH      = 0.0,
    volexact   = 0.0,
    areaexact  = 0.0,
    domeHexact = 0.0;
  double
    Herr   = 0.0,
    avHerr = 0.0,
    etaerr = 0.0;

  PetscScalar     dummy, z, dummy1, dummy2, dummy3, dummy4, dummy5;

  ierr = ice_thickness.begin_access(); CHKERRQ(ierr);
  if (testname == 'L') {
    ierr = vHexactL.begin_access(); CHKERRQ(ierr);
  }

  double
    seawater_density = config.get("sea_water_density"),
    ice_density      = config.get("ice_density"),
    Glen_n           = config.get("Glen_exponent"),
    standard_gravity = config.get("standard_gravity"),
    // enthalpy and pressure do not matter here
    B0, C,
    H0               = 600.0,
    v0               = grid.convert(300.0, "m/year", "m/second"),
    Q0               = H0 * v0;

  if (testname == 'V') {
    B0 = stress_balance->get_stressbalance()->get_flow_law()->hardness_parameter(0, 0);
    C  = pow(ice_density * standard_gravity * (1.0 - ice_density/seawater_density) / (4 * B0), 3);
  }

  // area of grid square in square km:
  const PetscScalar   a = grid.dx * grid.dy * 1e-3 * 1e-3;
  const PetscScalar   m = (2.0 * Glen_n + 2.0) / Glen_n;
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      if (ice_thickness(i,j) > 0) {
        area += a;
        vol += a * ice_thickness(i,j) * 1e-3;
      }
      PetscScalar xx = grid.x[i], yy = grid.y[j],
        r = grid.radius(i,j);
      switch (testname) {
      case 'A':
        exactA(r,&Hexact,&dummy);
        break;
      case 'B':
        exactB(time,r,&Hexact,&dummy);
        break;
      case 'C':
        exactC(time,r,&Hexact,&dummy);
        break;
      case 'D':
        exactD(time,r,&Hexact,&dummy);
        break;
      case 'E':
        exactE(xx,yy,&Hexact,&dummy,&dummy1,&dummy2,&dummy3);
        break;
      case 'F':
        if (r > LforFG - 1.0) {  // outside of sheet
          Hexact=0.0;
        } else {
          r=PetscMax(r,1.0);
          z=0.0;
          bothexact(0.0,r,&z,1,0.0,
                    &Hexact,&dummy,&dummy5,&dummy1,&dummy2,&dummy3,&dummy4);
        }
        break;
      case 'G':
        if (r > LforFG -1.0) {  // outside of sheet
          Hexact=0.0;
        } else {
          r=PetscMax(r,1.0);
          z=0.0;
          bothexact(time,r,&z,1,ApforG,
                    &Hexact,&dummy,&dummy5,&dummy1,&dummy2,&dummy3,&dummy4);
        }
        break;
      case 'H':
        exactH(f,time,r,&Hexact,&dummy);
        break;
      case 'K':
      case 'O':
        Hexact = 3000.0;
        break;
      case 'L':
        Hexact = vHexactL(i,j);
        break;
      case 'V':
        {
          Hexact = pow(4 * C / Q0 * xx + 1/pow(H0, 4), -0.25);
        }
        break;
      default:
        SETERRQ(grid.com, 1, "test must be A, B, C, D, E, F, G, H, K, L, or O");
      }

      if (Hexact > 0) {
        areaexact += a;
        volexact += a * Hexact * 1e-3;
      }
      if (i == (grid.Mx - 1)/2 && j == (grid.My - 1)/2) {
        domeH = ice_thickness(i,j);
        domeHexact = Hexact;
      }
      // compute maximum errors
      Herr = PetscMax(Herr,PetscAbsReal(ice_thickness(i,j) - Hexact));
      etaerr = PetscMax(etaerr,PetscAbsReal(pow(ice_thickness(i,j),m) - pow(Hexact,m)));
      // add to sums for average errors
      avHerr += PetscAbsReal(ice_thickness(i,j) - Hexact);
    }
  }

  ierr = ice_thickness.end_access(); CHKERRQ(ierr);
  if (testname == 'L') {
    ierr = vHexactL.end_access(); CHKERRQ(ierr);
  }

  // globalize (find errors over all processors)
  PetscScalar gvol, garea, gdomeH;
  ierr = PISMGlobalSum(&volexact, &gvolexact, grid.com); CHKERRQ(ierr);
  ierr = PISMGlobalMax(&domeHexact, &gdomeHexact, grid.com); CHKERRQ(ierr);
  ierr = PISMGlobalSum(&areaexact, &gareaexact, grid.com); CHKERRQ(ierr);

  ierr = PISMGlobalSum(&vol, &gvol, grid.com); CHKERRQ(ierr);
  ierr = PISMGlobalSum(&area, &garea, grid.com); CHKERRQ(ierr);
  volerr = PetscAbsReal(gvol - gvolexact);
  areaerr = PetscAbsReal(garea - gareaexact);

  ierr = PISMGlobalMax(&Herr, &gmaxHerr, grid.com); CHKERRQ(ierr);
  ierr = PISMGlobalSum(&avHerr, &gavHerr, grid.com); CHKERRQ(ierr);
  gavHerr = gavHerr/(grid.Mx*grid.My);
  ierr = PISMGlobalMax(&etaerr, &gmaxetaerr, grid.com); CHKERRQ(ierr);

  ierr = PISMGlobalMax(&domeH, &gdomeH, grid.com); CHKERRQ(ierr);
  centerHerr = PetscAbsReal(gdomeH - gdomeHexact);

  return 0;
}


PetscErrorCode IceCompModel::computeBasalVelocityErrors(PetscScalar &exactmaxspeed, PetscScalar &gmaxvecerr,
                                                        PetscScalar &gavvecerr, PetscScalar &gmaxuberr,
                                                        PetscScalar &gmaxvberr) {

  PetscErrorCode ierr;
  PetscScalar    maxvecerr, avvecerr, maxuberr, maxvberr;
  PetscScalar    ubexact,vbexact, dummy1,dummy2,dummy3;
  PISMVector2    bvel;

  if (testname != 'E')
    SETERRQ(grid.com, 1,"basal velocity errors only computable for test E\n");

  IceModelVec2V *vel_adv;
  ierr = stress_balance->get_2D_advective_velocity(vel_adv); CHKERRQ(ierr);

  ierr = vel_adv->begin_access(); CHKERRQ(ierr);
  ierr = ice_thickness.begin_access(); CHKERRQ(ierr);
  maxvecerr = 0.0; avvecerr = 0.0; maxuberr = 0.0; maxvberr = 0.0;
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; i++) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; j++) {
      if (ice_thickness(i,j) > 0.0) {
        PetscScalar xx = grid.x[i], yy = grid.y[j];
        exactE(xx,yy,&dummy1,&dummy2,&dummy3,&ubexact,&vbexact);
        // compute maximum errors
        const PetscScalar uberr = PetscAbsReal((*vel_adv)(i,j).u - ubexact);
        const PetscScalar vberr = PetscAbsReal((*vel_adv)(i,j).v - vbexact);
        maxuberr = PetscMax(maxuberr,uberr);
        maxvberr = PetscMax(maxvberr,vberr);
        const PetscScalar vecerr = sqrt(uberr*uberr + vberr*vberr);
        maxvecerr = PetscMax(maxvecerr,vecerr);
        avvecerr += vecerr;
      }
    }
  }
  ierr = ice_thickness.end_access(); CHKERRQ(ierr);
  ierr = vel_adv->end_access(); CHKERRQ(ierr);

  ierr = PISMGlobalMax(&maxuberr, &gmaxuberr, grid.com); CHKERRQ(ierr);
  ierr = PISMGlobalMax(&maxvberr, &gmaxvberr, grid.com); CHKERRQ(ierr);

  ierr = PISMGlobalMax(&maxvecerr, &gmaxvecerr, grid.com); CHKERRQ(ierr);
  ierr = PISMGlobalSum(&avvecerr, &gavvecerr, grid.com); CHKERRQ(ierr);
  gavvecerr = gavvecerr/(grid.Mx*grid.My);

  const PetscScalar xpeak = 450e3 * cos(25.0*(M_PI/180.0)),
                    ypeak = 450e3 * sin(25.0*(M_PI/180.0));
  exactE(xpeak,ypeak,&dummy1,&dummy2,&dummy3,&ubexact,&vbexact);
  exactmaxspeed = sqrt(ubexact*ubexact + vbexact*vbexact);
  return 0;
}


PetscErrorCode IceCompModel::additionalAtStartTimestep() {
  PetscErrorCode    ierr;

  if (exactOnly == PETSC_TRUE && testname != 'K')
    dt_force = config.get("maximum_time_step_years", "years", "seconds");

  if (testname == 'F' || testname == 'G') {
    ierr = getCompSourcesTestFG(); CHKERRQ(ierr);
  }

  return 0;
}


PetscErrorCode IceCompModel::additionalAtEndTimestep() {
  PetscErrorCode    ierr;

  if (testname == 'A' || testname == 'E') {
    ierr = reset_thickness_tests_AE(); CHKERRQ(ierr);
  }

  // do nothing at the end of the time step unless the user has asked for the
  // exact solution to overwrite the numerical solution
  if (exactOnly == PETSC_FALSE)
    return 0;

  // because user wants exact solution, fill gridded values from exact formulas;
  // important notes:
  //     (1) the numerical computation *has* already occurred, in run(),
  //           and we just overwrite it with the exact solution here
  //     (2) certain diagnostic quantities like dHdt are computed numerically,
  //           and not overwritten here; while cbar,csurf,cflx,wsurf are diagnostic
  //           quantities recomputed at the end of the run for writing into
  //           NetCDF, in particular dHdt is not recomputed before being written
  //           into the output file, so it is actually numerical
  switch (testname) {
  case 'A':
  case 'B':
  case 'C':
  case 'D':
  case 'H':
    ierr = fillSolnTestABCDH(); CHKERRQ(ierr);
    break;
  case 'E':
    ierr = fillSolnTestE(); CHKERRQ(ierr);
    break;
  case 'F':
  case 'G':
    ierr = fillSolnTestFG(); CHKERRQ(ierr); // see iCMthermo.cc
    break;
  case 'K':
    ierr = fillTemperatureSolnTestsKO(); CHKERRQ(ierr); // see iCMthermo.cc
    break;
  case 'O':
    ierr = fillTemperatureSolnTestsKO(); CHKERRQ(ierr); // see iCMthermo.cc
    ierr = fillBasalMeltRateSolnTestO(); CHKERRQ(ierr); // see iCMthermo.cc
    break;
  case 'L':
    ierr = fillSolnTestL(); CHKERRQ(ierr);
    break;
  default:
    SETERRQ(grid.com, 1,"unknown testname in IceCompModel");
  }

  return 0;
}


PetscErrorCode IceCompModel::summary(bool /* tempAndAge */) {
  //   we always show a summary at every step
  return IceModel::summary(true);
}


PetscErrorCode IceCompModel::reportErrors() {
  // geometry errors to report (for all tests except K and O):
  //    -- max thickness error
  //    -- average (at each grid point on whole grid) thickness error
  //    -- max (thickness)^(2n+2)/n error
  //    -- volume error
  //    -- area error
  // and temperature errors (for tests F & G & K & O):
  //    -- max T error over 3D domain of ice
  //    -- av T error over 3D domain of ice
  // and basal temperature errors (for tests F & G):
  //    -- max basal temp error
  //    -- average (at each grid point on whole grid) basal temp error
  // and bedrock temperature errors (for tests K & O):
  //    -- max Tb error over 3D domain of bedrock
  //    -- av Tb error over 3D domain of bedrock
  // and strain-heating (Sigma) errors (for tests F & G):
  //    -- max Sigma error over 3D domain of ice (in 10^-3 K a^-1)
  //    -- av Sigma error over 3D domain of ice (in 10^-3 K a^-1)
  // and basal melt rate error (for test O):
  //    -- max bmelt error over base of ice
  // and surface velocity errors (for tests F & G):
  //    -- max |<us,vs> - <usex,vsex>| error
  //    -- av |<us,vs> - <usex,vsex>| error
  //    -- max ws error
  //    -- av ws error
  // and basal sliding errors (for test E):
  //    -- max ub error
  //    -- max vb error
  //    -- max |<ub,vb> - <ubexact,vbexact>| error
  //    -- av |<ub,vb> - <ubexact,vbexact>| error

  PetscErrorCode  ierr;

  bool dont_report;
  ierr = PISMOptionsIsSet("-no_report", "Don't report numerical errors",
                          dont_report); CHKERRQ(ierr);

  if (dont_report)
    return 0;

  IceFlowLaw* flow_law = stress_balance->get_ssb_modifier()->get_flow_law();
  if (testname != 'V' &&
      IceFlowLawIsPatersonBuddCold(flow_law, config, EC) == false &&
      ((testname == 'F') || (testname == 'G'))) {
    ierr = verbPrintf(1, grid.com,
                      "pismv WARNING: flow law must be cold part of Paterson-Budd ('-siafd_flow_law arr')\n"
                      "   for reported errors in test %c to be meaningful!\n",
                      testname); CHKERRQ(ierr);
  }

  ierr = verbPrintf(1,grid.com,
     "NUMERICAL ERRORS evaluated at final time (relative to exact solution):\n");
  CHKERRQ(ierr);

  unsigned int start;
  std::string filename;
  bool netcdf_report, append;
  NCTimeseries err("N", "N", grid.get_unit_system());

  err.set_units("1");

  PIO nc(grid.com, "netcdf3", grid.get_unit_system()); // OK to use netcdf3

  ierr = PISMOptionsString("-report_file", "NetCDF error report file",
                           filename, netcdf_report); CHKERRQ(ierr);
  ierr = PISMOptionsIsSet("-append", "Append the NetCDF error report",
                          append); CHKERRQ(ierr);
  if (netcdf_report) {
    ierr = verbPrintf(2,grid.com, "Also writing errors to '%s'...\n", filename.c_str());
    CHKERRQ(ierr);

    // Find the number of records in this file:
    ierr = nc.open(filename, PISM_WRITE, append); CHKERRQ(ierr);
    ierr = nc.inq_dimlen("N", start); CHKERRQ(ierr);

    ierr = nc.write_global_attributes(global_attributes); CHKERRQ(ierr);

    // Write the dimension variable:
    ierr = nc.write_timeseries(err, (size_t)start, (double)(start + 1), PISM_INT); CHKERRQ(ierr);

    // Always write grid parameters:
    err.set_name("dx");
    ierr = err.set_units("meters"); CHKERRQ(ierr);
    ierr = nc.write_timeseries(err, (size_t)start, grid.dx); CHKERRQ(ierr);
    err.set_name("dy");
    ierr = nc.write_timeseries(err, (size_t)start, grid.dy); CHKERRQ(ierr);
    err.set_name("dz");
    ierr = nc.write_timeseries(err, (size_t)start, grid.dzMAX); CHKERRQ(ierr);

    // Always write the test name:
    err.clear_all_strings(); err.clear_all_doubles(); err.set_units("1");
    err.set_name("test");
    ierr = nc.write_timeseries(err, (size_t)start, (double)testname, PISM_BYTE); CHKERRQ(ierr);
  }

  // geometry (thickness, vol) errors if appropriate; reported in m except for relmaxETA
  if ((testname != 'K') && (testname != 'O')) {
    PetscScalar volexact, areaexact, domeHexact, volerr, areaerr, maxHerr, avHerr,
                maxetaerr, centerHerr;
    ierr = computeGeometryErrors(volexact,areaexact,domeHexact,
                                 volerr,areaerr,maxHerr,avHerr,maxetaerr,centerHerr);
            CHKERRQ(ierr);
    ierr = verbPrintf(1,grid.com,
            "geometry  :    prcntVOL        maxH         avH   relmaxETA\n");
            CHKERRQ(ierr);  // no longer reporting centerHerr
    const PetscScalar   m = (2.0 * flow_law->exponent() + 2.0) / flow_law->exponent();
    ierr = verbPrintf(1,grid.com, "           %12.6f%12.6f%12.6f%12.6f\n",
                      100*volerr/volexact, maxHerr, avHerr,
                      maxetaerr/pow(domeHexact,m)); CHKERRQ(ierr);

    if (netcdf_report) {
      err.clear_all_strings(); err.clear_all_doubles(); err.set_units("1");
      err.set_name("relative_volume");
      ierr = err.set_units("percent"); CHKERRQ(ierr);
      err.set_string("long_name", "relative ice volume error");
      ierr = nc.write_timeseries(err, (size_t)start, 100*volerr/volexact); CHKERRQ(ierr);

      err.set_name("relative_max_eta");
      ierr = err.set_units("1"); CHKERRQ(ierr);
      err.set_string("long_name", "relative $\\eta$ error");
      ierr = nc.write_timeseries(err, (size_t)start, maxetaerr/pow(domeHexact,m)); CHKERRQ(ierr);

      err.set_name("maximum_thickness");
      ierr = err.set_units("meters"); CHKERRQ(ierr);
      err.set_string("long_name", "maximum ice thickness error");
      ierr = nc.write_timeseries(err, (size_t)start, maxHerr); CHKERRQ(ierr);

      err.set_name("average_thickness");
      ierr = err.set_units("meters"); CHKERRQ(ierr);
      err.set_string("long_name", "average ice thickness error");
      ierr = nc.write_timeseries(err, (size_t)start, avHerr); CHKERRQ(ierr);
    }
  }

  // temperature errors for F and G
  if ((testname == 'F') || (testname == 'G')) {
    PetscScalar maxTerr, avTerr, basemaxTerr, baseavTerr, basecenterTerr;
    ierr = computeTemperatureErrors(maxTerr, avTerr); CHKERRQ(ierr);
    ierr = computeBasalTemperatureErrors(basemaxTerr, baseavTerr, basecenterTerr);
       CHKERRQ(ierr);
    ierr = verbPrintf(1,grid.com,
       "temp      :        maxT         avT    basemaxT     baseavT\n");
       CHKERRQ(ierr);  // no longer reporting   basecenterT
    ierr = verbPrintf(1,grid.com, "           %12.6f%12.6f%12.6f%12.6f\n",
       maxTerr, avTerr, basemaxTerr, baseavTerr); CHKERRQ(ierr);

    if (netcdf_report) {
      err.clear_all_strings(); err.clear_all_doubles(); err.set_units("1");
      err.set_name("maximum_temperature");
      ierr = err.set_units("Kelvin"); CHKERRQ(ierr);
      err.set_string("long_name", "maximum ice temperature error");
      ierr = nc.write_timeseries(err, (size_t)start, maxTerr); CHKERRQ(ierr);

      err.set_name("average_temperature");
      err.set_string("long_name", "average ice temperature error");
      ierr = nc.write_timeseries(err, (size_t)start, avTerr); CHKERRQ(ierr);

      err.set_name("maximum_basal_temperature");
      err.set_string("long_name", "maximum basal temperature error");
      ierr = nc.write_timeseries(err, (size_t)start, basemaxTerr); CHKERRQ(ierr);
      err.set_name("average_basal_temperature");
      err.set_string("long_name", "average basal temperature error");
      ierr = nc.write_timeseries(err, (size_t)start, baseavTerr); CHKERRQ(ierr);
    }

  } else if ((testname == 'K') || (testname == 'O')) {
    PetscScalar maxTerr, avTerr, maxTberr, avTberr;
    ierr = computeIceBedrockTemperatureErrors(maxTerr, avTerr, maxTberr, avTberr);
       CHKERRQ(ierr);
    ierr = verbPrintf(1,grid.com,
       "temp      :        maxT         avT       maxTb        avTb\n"); CHKERRQ(ierr);
    ierr = verbPrintf(1,grid.com, "           %12.6f%12.6f%12.6f%12.6f\n",
                  maxTerr, avTerr, maxTberr, avTberr); CHKERRQ(ierr);

    if (netcdf_report) {
      err.clear_all_strings(); err.clear_all_doubles(); err.set_units("1");
      err.set_name("maximum_temperature");
      ierr = err.set_units("Kelvin"); CHKERRQ(ierr);
      err.set_string("long_name", "maximum ice temperature error");
      ierr = nc.write_timeseries(err, (size_t)start, maxTerr); CHKERRQ(ierr);

      err.set_name("average_temperature");
      err.set_string("long_name", "average ice temperature error");
      ierr = nc.write_timeseries(err, (size_t)start, avTerr); CHKERRQ(ierr);

      err.set_name("maximum_bedrock_temperature");
      err.set_string("long_name", "maximum bedrock temperature error");
      ierr = nc.write_timeseries(err, (size_t)start, maxTberr); CHKERRQ(ierr);

      err.set_name("average_bedrock_temperature");
      err.set_string("long_name", "average bedrock temperature error");
      ierr = nc.write_timeseries(err, (size_t)start, avTberr); CHKERRQ(ierr);
    }
  }

  // strain_heating errors if appropriate; reported in 10^6 J/(s m^3)
  if ((testname == 'F') || (testname == 'G')) {
    PetscScalar max_strain_heating_error, av_strain_heating_error;
    ierr = compute_strain_heating_errors(max_strain_heating_error, av_strain_heating_error); CHKERRQ(ierr);
    ierr = verbPrintf(1,grid.com,
       "Sigma     :      maxSig       avSig\n"); CHKERRQ(ierr);
    ierr = verbPrintf(1,grid.com, "           %12.6f%12.6f\n",
                  max_strain_heating_error*1.0e6, av_strain_heating_error*1.0e6); CHKERRQ(ierr);

    if (netcdf_report) {
      err.clear_all_strings(); err.clear_all_doubles(); err.set_units("1");
      err.set_name("maximum_sigma");
      ierr = err.set_units("J s-1 m-3"); CHKERRQ(ierr);
      ierr = err.set_glaciological_units("1e6 J s-1 m-3"); CHKERRQ(ierr);
      err.set_string("long_name", "maximum strain heating error");
      ierr = nc.write_timeseries(err, (size_t)start, max_strain_heating_error); CHKERRQ(ierr);

      err.set_name("average_sigma");
      err.set_string("long_name", "average strain heating error");
      ierr = nc.write_timeseries(err, (size_t)start, av_strain_heating_error); CHKERRQ(ierr);
    }
  }

  // surface velocity errors if exact values are available; reported in m/year
  if ((testname == 'F') || (testname == 'G')) {
    PetscScalar maxUerr, avUerr, maxWerr, avWerr;
    ierr = computeSurfaceVelocityErrors(maxUerr, avUerr, maxWerr, avWerr); CHKERRQ(ierr);
    ierr = verbPrintf(1,grid.com,
       "surf vels :     maxUvec      avUvec        maxW         avW\n"); CHKERRQ(ierr);
    ierr = verbPrintf(1,grid.com, "           %12.6f%12.6f%12.6f%12.6f\n",
                  grid.convert(maxUerr, "m/second", "m/year"), grid.convert(avUerr, "m/second", "m/year"), grid.convert(maxWerr, "m/second", "m/year"), grid.convert(avWerr, "m/second", "m/year")); CHKERRQ(ierr);

    if (netcdf_report) {
      err.clear_all_strings(); err.clear_all_doubles(); err.set_units("1");
      err.set_name("maximum_surface_velocity");
      err.set_string("long_name", "maximum ice surface horizontal velocity error");
      ierr = err.set_units("m/s"); CHKERRQ(ierr);
      ierr = err.set_glaciological_units("meters/year"); CHKERRQ(ierr);
      ierr = nc.write_timeseries(err, (size_t)start, maxUerr); CHKERRQ(ierr);

      err.set_name("average_surface_velocity");
      err.set_string("long_name", "average ice surface horizontal velocity error");
      ierr = nc.write_timeseries(err, (size_t)start, avUerr); CHKERRQ(ierr);

      err.set_name("maximum_surface_w");
      err.set_string("long_name", "maximum ice surface vertical velocity error");
      ierr = nc.write_timeseries(err, (size_t)start, maxWerr); CHKERRQ(ierr);

      err.set_name("average_surface_w");
      err.set_string("long_name", "average ice surface vertical velocity error");
      ierr = nc.write_timeseries(err, (size_t)start, avWerr); CHKERRQ(ierr);
    }
  }

  // basal velocity errors if appropriate; reported in m/year except prcntavvec
  if (testname == 'E') {
    PetscScalar exactmaxspeed, maxvecerr, avvecerr, maxuberr, maxvberr;
    ierr = computeBasalVelocityErrors(exactmaxspeed,
                          maxvecerr,avvecerr,maxuberr,maxvberr); CHKERRQ(ierr);
    ierr = verbPrintf(1,grid.com,
       "base vels :  maxvector   avvector  prcntavvec     maxub     maxvb\n");
       CHKERRQ(ierr);
    ierr = verbPrintf(1,grid.com, "           %11.4f%11.5f%12.5f%10.4f%10.4f\n",
                  grid.convert(maxvecerr, "m/second", "m/year"), grid.convert(avvecerr, "m/second", "m/year"),
                  (avvecerr/exactmaxspeed)*100.0,
                  grid.convert(maxuberr, "m/second", "m/year"), grid.convert(maxvberr, "m/second", "m/year")); CHKERRQ(ierr);

    if (netcdf_report) {
      err.clear_all_strings(); err.clear_all_doubles(); err.set_units("1");
      err.set_name("maximum_basal_velocity");
      ierr = err.set_units("m/s"); CHKERRQ(ierr);
      ierr = err.set_glaciological_units("meters/year"); CHKERRQ(ierr);
      ierr = nc.write_timeseries(err, (size_t)start, maxvecerr); CHKERRQ(ierr);

      err.set_name("average_basal_velocity");
      ierr = nc.write_timeseries(err, (size_t)start, avvecerr); CHKERRQ(ierr);
      err.set_name("maximum_basal_u");
      ierr = nc.write_timeseries(err, (size_t)start, maxuberr); CHKERRQ(ierr);
      err.set_name("maximum_basal_v");
      ierr = nc.write_timeseries(err, (size_t)start, maxvberr); CHKERRQ(ierr);

      err.clear_all_strings(); err.clear_all_doubles(); err.set_units("1");
      err.set_name("relative_basal_velocity");
      ierr = err.set_units("percent"); CHKERRQ(ierr);
      ierr = nc.write_timeseries(err, (size_t)start, (avvecerr/exactmaxspeed)*100); CHKERRQ(ierr);
    }
  }

  // basal melt rate errors if appropriate; reported in m/year
  if (testname == 'O') {
    PetscScalar maxbmelterr, minbmelterr;
    ierr = computeBasalMeltRateErrors(maxbmelterr, minbmelterr); CHKERRQ(ierr);
    if (maxbmelterr != minbmelterr) {
       ierr = verbPrintf(1,grid.com,
                         "IceCompModel WARNING: unexpected Test O situation: max and min of bmelt error\n"
                         "  are different: maxbmelterr = %f, minbmelterr = %f\n",
                         grid.convert(maxbmelterr, "m/second", "m/year"),
                         grid.convert(minbmelterr, "m/second", "m/year")); CHKERRQ(ierr);
    }
    ierr = verbPrintf(1,grid.com,
       "basal melt:  max\n"); CHKERRQ(ierr);
    ierr = verbPrintf(1,grid.com, "           %11.5f\n",
                      grid.convert(maxbmelterr, "m/second", "m/year")); CHKERRQ(ierr);

    if (netcdf_report) {
      err.clear_all_strings(); err.clear_all_doubles(); err.set_units("1");
      err.set_name("maximum_basal_melt_rate");
      ierr = err.set_units("m/s"); CHKERRQ(ierr);
      ierr = err.set_glaciological_units("meters/year"); CHKERRQ(ierr);
      ierr = nc.write_timeseries(err, (size_t)start, maxbmelterr); CHKERRQ(ierr);
    }
  }

  if (netcdf_report) {
    ierr = nc.close(); CHKERRQ(ierr);
  }

  ierr = verbPrintf(1,grid.com, "NUM ERRORS DONE\n");  CHKERRQ(ierr);
  return 0;
}

//! \brief Initialize test V.
/*
 Try

 pismv -test V -y 1000 -part_grid -ssa_method fd -cfbc -o fig4-blue.nc
 pismv -test V -y 1000 -part_grid -ssa_method fd -o fig4-green.nc

 to try to reproduce Figure 4.

 Try

 pismv -test V -y 3000 -ssa_method fd -cfbc -o fig5.nc -thickness_calving_threshold 250 -part_grid

 with -Mx 51, -Mx 101, -Mx 201 for figure 5,

 pismv -test V -y 300 -ssa_method fd -o fig6-ab.nc

 for 6a and 6b,

 pismv -test V -y 300 -ssa_method fd -cfbc -part_grid -o fig6-cd.nc

 for 6c and 6d,

 pismv -test V -y 300 -ssa_method fd -cfbc -part_grid -part_redist -o fig6-ef.nc

 for 6e and 6f.

 */
PetscErrorCode IceCompModel::test_V_init() {
  PetscErrorCode ierr;

  // initialize the bed topography
  ierr = bed_topography.set(-1000); CHKERRQ(ierr);

  // set SSA boundary conditions:
  PetscReal upstream_velocity = grid.convert(300.0, "m/year", "m/second"),
    upstream_thk = 600.0;

  ierr = ice_thickness.begin_access(); CHKERRQ(ierr);
  ierr = vBCMask.begin_access(); CHKERRQ(ierr);
  ierr = vBCvel.begin_access(); CHKERRQ(ierr);
  for (PetscInt   i = grid.xs; i < grid.xs+grid.xm; ++i) {
    for (PetscInt j = grid.ys; j < grid.ys+grid.ym; ++j) {
      if (i <= 2) {
        vBCMask(i,j) = 1;
        vBCvel(i,j)  = PISMVector2(upstream_velocity, 0.0);
        ice_thickness(i, j) = upstream_thk;
      } else {
        vBCMask(i,j) = 0;
        vBCvel(i,j)  = PISMVector2(0.0, 0.0);
        ice_thickness(i, j) = 0;
      }
    }
  }
  ierr = vBCvel.end_access(); CHKERRQ(ierr);
  ierr = vBCMask.end_access(); CHKERRQ(ierr);
  ierr = ice_thickness.end_access(); CHKERRQ(ierr);

  ierr = vBCMask.update_ghosts(); CHKERRQ(ierr);

  ierr = vBCvel.update_ghosts(); CHKERRQ(ierr);

  ierr = ice_thickness.update_ghosts(); CHKERRQ(ierr);

  return 0;
}
