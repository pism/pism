// Copyright (C) 2004-2010 Jed Brown, Ed Bueler and Constantine Khroulev
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

#include <cmath>
#include <cstring>
#include <petscda.h>

#include <vector>     // STL vector container; sortable; used in test L
#include <algorithm>  // required by sort(...) in test L

#include "tests/exactTestsABCDE.h"
#include "tests/exactTestsFG.h" 
#include "tests/exactTestH.h" 
#include "tests/exactTestL.h" 

#include "iceCompModel.hh"

const PetscScalar IceCompModel::ablationRateOutside = 0.02; // m/a

IceCompModel::IceCompModel(IceGrid &g, NCConfigVariable &conf, NCConfigVariable &conf_overrides, int mytest)
  : IceModel(g, conf, conf_overrides), tgaIce(NULL) {
  
  // note lots of defaults are set by the IceModel constructor defaults for
  // IceCompModel:
  testname = mytest;
  exactOnly = PETSC_FALSE;
  bedrock_is_ice_forK = PETSC_FALSE;
  
  // Override some defaults from parent class
  config.set("enhancement_factor", 1.0);

  // set values of flags in run() 
  config.set_flag("do_mass_conserve", true);
  config.set_flag("use_ssa_velocity", false);
  config.set_flag("include_bmr_in_continuity", false);
  config.set_flag("do_plastic_till", false);
}


IceCompModel::~IceCompModel() {
}


PetscErrorCode IceCompModel::createVecs() {
  PetscErrorCode ierr;

  ierr = IceModel::createVecs(); CHKERRQ(ierr);

  ierr = vHexactL.create(grid, "HexactL", true); CHKERRQ(ierr);

  ierr = SigmaComp3.create(grid,"SigmaComp", false); CHKERRQ(ierr);
  ierr = SigmaComp3.set_attrs("internal","rate of compensatory strain heating in ice",
			      "W m-3", ""); CHKERRQ(ierr);

  // this ensures that these variables are saved to an output file and are read
  // back in if -i option is used (they are "model_state", in a sense, since
  // PSDummy is used):
  ierr = artm.set_attr("pism_intent", "model_state"); CHKERRQ(ierr);
  ierr = acab.set_attr("pism_intent", "model_state"); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode IceCompModel::set_grid_defaults() {
  PetscErrorCode ierr;

  // This sets the defaults for each test; command-line options can override this.

  // equal spacing is the default for all the tests except K
  grid.ice_vertical_spacing = EQUAL;
  grid.bed_vertical_spacing = EQUAL;

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
    // use 2000km by 2000km by 4000m rectangular domain, but make truely periodic
    grid.Mbz = 2;
    grid.Lx = grid.Ly = 1000e3;
    grid.Lz = 4000;
    grid.Lbz = 1000;		// this is important
    grid.periodicity = XY_PERIODIC;
    grid.ice_vertical_spacing = QUADRATIC;
    grid.bed_vertical_spacing = QUADRATIC;
    break;
  default:
    ierr = PetscPrintf(grid.com, "IceCompModel ERROR : desired test not implemented\n");
    CHKERRQ(ierr);
    PetscEnd();
  }

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

  // These ifs are here (and not in the constructor or init_physics()) because
  // testname actually comes from a command-line *and* because command-line
  // options should be able to override parameter values set here.

  if (testname == 'H') {
    config.set_flag("do_bed_deformation", true);
    config.set_flag("do_bed_iso", true);
  } else
    config.set_flag("do_bed_deformation", false);

  if ((testname == 'F') || (testname == 'G') || (testname == 'K')) {
    config.set_flag("do_temp", true);
    // essentially turn off run-time reporting of extremely low computed
    // temperatures; *they will be reported as errors* anyway
    config.set("global_min_allowed_temp", 0.0);
    config.set("max_low_temp_count", 1000000);
  } else
    config.set_flag("do_temp", false);

  if ((testname == 'A') || (testname == 'E')) {
    config.set_flag("is_dry_simulation", true);
    config.set_flag("ocean_kill", true);
  } else {
    config.set_flag("is_dry_simulation", true);
    config.set_flag("ocean_kill", false);
  }

  // special considerations for K wrt thermal bedrock and pressure-melting
  // (flag thermalBedrock was removed by CK around r783, because it was not used)
  if (testname == 'K') {
    allowAboveMelting = PETSC_FALSE;
    reportPATemps = PETSC_TRUE;
  } else {
    // note temps are generally allowed to go above pressure melting in verify
    allowAboveMelting = PETSC_TRUE;
    reportPATemps = PETSC_FALSE;
  }

  ierr = IceModel::setFromOptions();CHKERRQ(ierr);
  return 0;
}


PetscErrorCode IceCompModel::init_physics() {
  PetscErrorCode ierr;

  // Set the default for IceCompModel:
  ierr = iceFactory.setType(ICE_ARR); CHKERRQ(ierr);

  // Let the base class version read the options (possibly overriding the
  // default set above) and create the IceFlowLaw object.
  ierr = IceModel::init_physics(); CHKERRQ(ierr);
  
  // check on whether the options (already checked) chose the right IceFlowLaw for verification;
  //   need to have a tempFromSoftness() procedure as well as the need for the right
  //   flow law to have the errors make sense
  tgaIce = dynamic_cast<ThermoGlenArrIce*>(ice);
  if (!tgaIce) SETERRQ(1,"IceCompModel requires ThermoGlenArrIce or a derived class");
  if (IceFlowLawIsPatersonBuddCold(ice, config) == PETSC_FALSE) {
    ierr = verbPrintf(1, grid.com, 
       "WARNING: user set -gk; default flow law should be -ice_type arr for IceCompModel\n");
    CHKERRQ(ierr);
  }

  f = ice->rho / config.get("lithosphere_density");  // for simple isostasy

  if (testname != 'K') {
    // now make bedrock have same material properties as ice
    // (note Mbz=1 also, by default, but want ice/rock interface to see
    // pure ice from the point of view of applying geothermal boundary
    // condition, especially in tests F and G)
    config.set("bedrock_thermal_density", tgaIce->rho);
    config.set("bedrock_thermal_conductivity", tgaIce->k);
    config.set("bedrock_thermal_specific_heat_capacity", tgaIce->c_p);
  }

  bool do_bed_deformation = config.get_flag("do_bed_deformation"),
    do_bed_iso = config.get_flag("do_bed_iso");

  if ( (testname == 'H') && ((do_bed_deformation == PETSC_FALSE) || (do_bed_iso == PETSC_FALSE)) ) {
    ierr = verbPrintf(1,grid.com, 
           "IceCompModel WARNING: Test H should be run with option\n"
           "  -bed_def_iso  for the reported errors to be correct.\n"); CHKERRQ(ierr);
  }

  // this switch changes Test K to make material properties for bedrock the same as for ice
  bool biiSet;
  ierr = PISMOptionsIsSet("-bedrock_is_ice", biiSet); CHKERRQ(ierr);
  if (biiSet == PETSC_TRUE) {
    if (testname == 'K') {
      ierr = verbPrintf(1,grid.com,
         "setting material properties of bedrock to those of ice in Test K\n"); CHKERRQ(ierr);
      config.set("bedrock_thermal_density", tgaIce->rho);
      config.set("bedrock_thermal_conductivity", tgaIce->k);
      config.set("bedrock_thermal_specific_heat_capacity", tgaIce->c_p);
      bedrock_is_ice_forK = PETSC_TRUE;
    } else {
      ierr = verbPrintf(1,grid.com,
         "IceCompModel WARNING: option -bedrock_is_ice ignored; only applies to Test K\n");
         CHKERRQ(ierr);
    }
  }

  return 0;
}

PetscErrorCode IceCompModel::set_vars_from_options() {
  PetscErrorCode ierr;

  // -boot_from command-line option is not allowed here.
  ierr = stop_if_set(grid.com, "-boot_from"); CHKERRQ(ierr);

  ierr = SigmaComp3.set(0.0); CHKERRQ(ierr);
  ierr = vtillphi.set(config.get("default_till_phi")); CHKERRQ(ierr);

  ierr = verbPrintf(3,grid.com, "initializing Test %c from formulas ...\n",testname);  CHKERRQ(ierr);

  // all have no uplift or Hmelt
  ierr = vuplift.set(0.0); CHKERRQ(ierr);
  ierr = vHmelt.set(0.0); CHKERRQ(ierr);
  ierr = vbasalMeltRate.set(0.0); CHKERRQ(ierr);

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
    ierr = initTestK(); CHKERRQ(ierr);  // see iCMthermo.cc
    break;
  case 'L':
    ierr = initTestL(); CHKERRQ(ierr);
    break;
  default:  SETERRQ(1,"Desired test not implemented by IceCompModel.\n");
  }

  return 0;
}


void IceCompModel::mapcoords(const PetscInt i, const PetscInt j,
                             PetscScalar &x, PetscScalar &y, PetscScalar &r) {
  // compute x,y,r on grid from i,j
  PetscScalar ifrom0, jfrom0;

  ifrom0=static_cast<PetscScalar>(i)-static_cast<PetscScalar>(grid.Mx - 1)/2.0;
  jfrom0=static_cast<PetscScalar>(j)-static_cast<PetscScalar>(grid.My - 1)/2.0;
  x=grid.dx*ifrom0;
  y=grid.dy*jfrom0;
  r = sqrt(PetscSqr(x) + PetscSqr(y));
}


//! Reimplement IceModel::basalVelocitySIA(), for Test E.  Gives zero SIA-type sliding in other tests.
PetscScalar IceCompModel::basalVelocitySIA(PetscScalar xIN, PetscScalar yIN,
                                           PetscScalar H, PetscScalar /*T*/,
                                           PetscScalar /*alpha*/, PetscScalar /*muIN*/,
					   PetscScalar /*min_T*/) const {

  if (testname == 'E') {
    const PetscScalar r1 = 200e3, r2 = 700e3,   /* define region of sliding */
                      theta1 = 10 * (pi/180), theta2 = 40 * (pi/180);
    const PetscScalar x = fabs(xIN), y = fabs(yIN);
    const PetscScalar r = sqrt(x * x + y * y);
    PetscScalar       theta;
    if (x < 1.0)
      theta = pi / 2.0;
    else
      theta = atan(y / x);
  
    if ((r > r1) && (r < r2) && (theta > theta1) && (theta < theta2)) {
      // now INSIDE sliding region
      const PetscScalar rbot = (r2 - r1) * (r2 - r1),
                        thetabot = (theta2 - theta1) * (theta2 - theta1);
      const PetscScalar mu_max = 2.5e-11; /* Pa^-1 m s^-1; max sliding coeff */
      PetscScalar muE = mu_max * (4.0 * (r - r1) * (r2 - r) / rbot) 
                               * (4.0 * (theta - theta1) * (theta2 - theta) / thetabot);
      return muE * tgaIce->rho * standard_gravity * H;
    } else
      return 0.0;
  } else
    return 0.0;  // zero sliding for other tests
}


PetscErrorCode IceCompModel::initTestABCDEH() {
  PetscErrorCode  ierr;
  PetscScalar     A0, T0, **H, **accum, dummy1, dummy2, dummy3;
  const PetscScalar LforAE = 750e3; // m

  // compute T so that A0 = A(T) = Acold exp(-Qcold/(R T))  (i.e. for ThermoGlenArrIce);
  // set all temps to this constant
  A0 = 1.0e-16/secpera;    // = 3.17e-24  1/(Pa^3 s);  (EISMINT value) flow law parameter
  T0 = tgaIce->tempFromSoftness(A0);
  ierr = artm.set(T0); CHKERRQ(ierr);
  ierr =   T3.set(T0); CHKERRQ(ierr);
  ierr =  Tb3.set(T0); CHKERRQ(ierr);
  ierr = vGhf.set(Ggeo); CHKERRQ(ierr);
  
  ierr = vMask.set(MASK_SHEET); CHKERRQ(ierr);
  if (testname == 'E') { // value is not used by IceCompModel::basalVelocitySIA(),
    config.set("mu_sliding", 1.0); //    but this acts as flag to allow sliding
  } else {
    config.set("mu_sliding", 0.0);
  }

  ierr = acab.get_array(accum); CHKERRQ(ierr);
  ierr = vH.get_array(H); CHKERRQ(ierr);
  if ((testname == 'A') || (testname == 'E')) {
    ierr = vMask.begin_access(); CHKERRQ(ierr);
  }
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      PetscScalar r,xx,yy;
      mapcoords(i,j,xx,yy,r);
      switch (testname) {
        case 'A':
          exactA(r,&H[i][j],&accum[i][j]);
          if (r >= LforAE)
            vMask(i,j) = MASK_OCEAN_AT_TIME_0;
          break;
        case 'B':
          exactB(grid.year*secpera,r,&H[i][j],&accum[i][j]);
          break;
        case 'C':
          exactC(grid.year*secpera,r,&H[i][j],&accum[i][j]);
          break;
        case 'D':
          exactD(grid.year*secpera,r,&H[i][j],&accum[i][j]);
          break;
        case 'E':
          exactE(xx,yy,&H[i][j],&accum[i][j],&dummy1,&dummy2,&dummy3);
          if (r >= LforAE)
            vMask(i,j) = MASK_OCEAN_AT_TIME_0;
          break;
        case 'H':
          exactH(f,grid.year*secpera,r,&H[i][j],&accum[i][j]);
          break;
        default:  SETERRQ(1,"test must be A, B, C, D, E, or H");
      }
    }
  }
  ierr = acab.end_access(); CHKERRQ(ierr);
  ierr = vH.end_access(); CHKERRQ(ierr);
  if ((testname == 'A') || (testname == 'E')) {
    ierr = vMask.end_access(); CHKERRQ(ierr);
  }

  ierr = vH.beginGhostComm(); CHKERRQ(ierr);
  ierr = vH.endGhostComm(); CHKERRQ(ierr);

  if (testname == 'H') {
    ierr = vH.copy_to(vh); CHKERRQ(ierr);
    ierr = vh.scale(1-f); CHKERRQ(ierr);
    ierr = vH.copy_to(vbed); CHKERRQ(ierr);
    ierr = vbed.scale(-f); CHKERRQ(ierr);
  } else {  // flat bed case otherwise
    ierr = vH.copy_to(vh); CHKERRQ(ierr);
    ierr = vbed.set(0.0); CHKERRQ(ierr);
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
  bool operator()(rgrid a, rgrid b) { return (a.r >= b.r); }
};


PetscErrorCode IceCompModel::initTestL() {
  PetscErrorCode  ierr;
  PetscScalar     A0, T0, **H, **accum, **bed;

  if (testname != 'L')  { SETERRQ(1,"test must be 'L'"); }
  
  // compute T so that A0 = A(T) = Acold exp(-Qcold/(R T))  (i.e. for ThermoGlenArrIce);
  // set all temps to this constant
  A0 = 1.0e-16/secpera;    // = 3.17e-24  1/(Pa^3 s);  (EISMINT value) flow law parameter
  T0 = tgaIce->tempFromSoftness(A0);
  ierr = artm.set(T0); CHKERRQ(ierr);
  ierr =   T3.set(T0); CHKERRQ(ierr);
  ierr =  Tb3.set(T0); CHKERRQ(ierr);
  ierr = vGhf.set(Ggeo); CHKERRQ(ierr);
  
  ierr = vMask.set(MASK_SHEET); CHKERRQ(ierr);
  config.set("mu_sliding", 0.0);  // note reimplementation of basalVelocitySIA() in IceCompModel

  // setup to evaluate test L; requires solving an ODE numerically using sorted list
  //   of radii, sorted in decreasing radius order
  const int  MM = grid.xm * grid.ym;

  std::vector<rgrid> rrv(MM);  // destructor at end of scope
  for (PetscInt i = 0; i < grid.xm; i++) {
    for (PetscInt j = 0; j < grid.ym; j++) {
      const PetscInt  k = i * grid.ym + j;
      rrv[k].i = i + grid.xs;  rrv[k].j = j + grid.ys;
      PetscScalar  junkx, junky;
      mapcoords(rrv[k].i, rrv[k].j, junkx, junky, rrv[k].r);
    }
  }
  std::sort(rrv.begin(), rrv.end(), rgridReverseSort()); // so rrv[k].r > rrv[k+1].r

  // get soln to test L at these radii; solves ODE only once (on each processor)
  double *rr, *HH, *bb, *aa;
  rr = new double[MM];
  for (PetscInt k = 0; k < MM; k++)
    rr[k] = rrv[k].r;
  HH = new double[MM];  bb = new double[MM];  aa = new double[MM];
  ierr = exactL_list(rr, MM, HH, bb, aa);
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
  delete [] rr;
  
  ierr = acab.get_array(accum); CHKERRQ(ierr);
  ierr = vH.get_array(H); CHKERRQ(ierr);
  ierr = vbed.get_array(bed); CHKERRQ(ierr);
  for (PetscInt k = 0; k < MM; k++) {
    H    [rrv[k].i][rrv[k].j] = HH[k];
    bed  [rrv[k].i][rrv[k].j] = bb[k];
    accum[rrv[k].i][rrv[k].j] = aa[k];
  }
  ierr = acab.end_access(); CHKERRQ(ierr);
  ierr = vH.end_access(); CHKERRQ(ierr);
  ierr = vbed.end_access(); CHKERRQ(ierr);
  delete [] HH;  delete [] bb;  delete [] aa;

  ierr = vH.beginGhostComm(); CHKERRQ(ierr);
  ierr = vH.endGhostComm(); CHKERRQ(ierr);
  ierr = vbed.beginGhostComm(); CHKERRQ(ierr);
  ierr = vbed.endGhostComm(); CHKERRQ(ierr);

  // store copy of vH for "-eo" runs and for evaluating geometry errors
  ierr = vH.copy_to(vHexactL); CHKERRQ(ierr);

  // set surface to H+b
  ierr = vH.add(1.0, vbed, vh); CHKERRQ(ierr);
  ierr = vh.beginGhostComm(); CHKERRQ(ierr);
  ierr = vh.endGhostComm(); CHKERRQ(ierr);
  return 0;
}


PetscErrorCode IceCompModel::getCompSourcesTestCDH() {
  PetscErrorCode  ierr;
  PetscScalar     **accum, dummy;

  // before flow step, set accumulation from exact values;
  ierr = acab.get_array(accum); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      PetscScalar r,xx,yy;
      mapcoords(i,j,xx,yy,r);
      switch (testname) {
        case 'C':
          exactC(grid.year*secpera,r,&dummy,&accum[i][j]);
          break;
        case 'D':
          exactD(grid.year*secpera,r,&dummy,&accum[i][j]);
          break;
        case 'H':
          exactH(f,grid.year*secpera,r,&dummy,&accum[i][j]);
          break;
        default:  SETERRQ(1,"testname must be C, D, or H");
      }
    }
  }
  ierr = acab.end_access(); CHKERRQ(ierr);
  return 0;
}


PetscErrorCode IceCompModel::fillSolnTestABCDH() {
  PetscErrorCode  ierr;
  PetscScalar     **H, **accum;

  ierr = acab.get_array(accum); CHKERRQ(ierr);
  ierr = vH.get_array(H); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      PetscScalar r,xx,yy;
      mapcoords(i,j,xx,yy,r);
      switch (testname) {
        case 'A':
          exactA(r,&H[i][j],&accum[i][j]);
          break;
        case 'B':
          exactB(grid.year*secpera,r,&H[i][j],&accum[i][j]);
          break;
        case 'C':
          exactC(grid.year*secpera,r,&H[i][j],&accum[i][j]);
          break;
        case 'D':
          exactD(grid.year*secpera,r,&H[i][j],&accum[i][j]);
          break;
        case 'H':
          exactH(f,grid.year*secpera,r,&H[i][j],&accum[i][j]);
          break;
        default:  SETERRQ(1,"test must be A, B, C, D, or H");
      }
    }
  }

  ierr = acab.end_access(); CHKERRQ(ierr);
  ierr = vH.end_access(); CHKERRQ(ierr);

  ierr = vH.beginGhostComm(); CHKERRQ(ierr);
  ierr = vH.endGhostComm(); CHKERRQ(ierr);

  if (testname == 'H') {
    ierr = vH.copy_to(vh); CHKERRQ(ierr);
    ierr = vh.scale(1-f); CHKERRQ(ierr);
    ierr = vH.copy_to(vbed); CHKERRQ(ierr);
    ierr = vbed.scale(-f); CHKERRQ(ierr);
    ierr = vbed.beginGhostComm(); CHKERRQ(ierr);
    ierr = vbed.endGhostComm(); CHKERRQ(ierr);
  } else {
    ierr = vH.copy_to(vh); CHKERRQ(ierr);
  }
  ierr = vh.beginGhostComm(); CHKERRQ(ierr);
  ierr = vh.endGhostComm(); CHKERRQ(ierr);
  return 0;
}


PetscErrorCode IceCompModel::fillSolnTestE() {
  PetscErrorCode  ierr;
  PetscScalar     **H, **accum, **ub, **vb, dummy;

  ierr = acab.get_array(accum); CHKERRQ(ierr);
  ierr = vH.get_array(H); CHKERRQ(ierr);
  ierr = vub.get_array(ub); CHKERRQ(ierr);
  ierr = vvb.get_array(vb); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      PetscScalar r,xx,yy;
      mapcoords(i,j,xx,yy,r);
      exactE(xx,yy,&H[i][j],&accum[i][j],&dummy,&ub[i][j],&vb[i][j]);
    }
  }
  ierr = acab.end_access(); CHKERRQ(ierr);
  ierr = vH.end_access(); CHKERRQ(ierr);
  ierr = vub.end_access(); CHKERRQ(ierr);
  ierr = vvb.end_access(); CHKERRQ(ierr);

  ierr = vH.beginGhostComm(); CHKERRQ(ierr);
  ierr = vH.endGhostComm(); CHKERRQ(ierr);
  ierr = vH.copy_to(vh); CHKERRQ(ierr);

  ierr = vub.beginGhostComm(); CHKERRQ(ierr);
  ierr = vub.endGhostComm(); CHKERRQ(ierr);
  ierr = vvb.beginGhostComm(); CHKERRQ(ierr);
  ierr = vvb.endGhostComm(); CHKERRQ(ierr);

  return 0;
}


PetscErrorCode IceCompModel::fillSolnTestL() {
  PetscErrorCode  ierr;

  ierr = vHexactL.beginGhostComm(); CHKERRQ(ierr);
  ierr = vHexactL.endGhostComm(); CHKERRQ(ierr);
  ierr = vH.copy_from(vHexactL);

  ierr = vbed.add(1.0, vH, vh);	CHKERRQ(ierr); //  h = H + bed = 1 * H + bed
  ierr = vh.beginGhostComm(); CHKERRQ(ierr);
  ierr = vh.endGhostComm(); CHKERRQ(ierr);

  // note bed was filled at initialization and hasn't changed
  return 0;
}


PetscErrorCode IceCompModel::computeGeometryErrors(
      PetscScalar &gvolexact, PetscScalar &gareaexact, PetscScalar &gdomeHexact,
      PetscScalar &volerr, PetscScalar &areaerr,
      PetscScalar &gmaxHerr, PetscScalar &gavHerr, PetscScalar &gmaxetaerr,
      PetscScalar &centerHerr) {
  // compute errors in thickness, eta=thickness^{(2n+2)/n}, volume, area
  
  PetscErrorCode  ierr;
  PetscScalar     **H, **HexactL;
  PetscScalar     Hexact, vol, area, domeH, volexact, areaexact, domeHexact;
  PetscScalar     Herr, avHerr, etaerr;

  PetscScalar     dummy, z, dummy1, dummy2, dummy3, dummy4, dummy5;

  ierr = vH.get_array(H); CHKERRQ(ierr);
  if (testname == 'L') {
    ierr = vHexactL.get_array(HexactL); CHKERRQ(ierr);
  }

  vol = 0; area = 0; domeH = 0;
  volexact = 0; areaexact = 0; domeHexact = 0;
  Herr = 0; avHerr=0; etaerr = 0;

  // area of grid square in square km:
  const PetscScalar   a = grid.dx * grid.dy * 1e-3 * 1e-3;
  const PetscScalar   m = (2.0 * tgaIce->exponent() + 2.0) / tgaIce->exponent();
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      if (H[i][j] > 0) {
        area += a;
        vol += a * H[i][j] * 1e-3;
      }
      PetscScalar r,xx,yy;
      mapcoords(i,j,xx,yy,r);
      switch (testname) {
        case 'A':
          exactA(r,&Hexact,&dummy);
          break;
        case 'B':
          exactB(grid.year*secpera,r,&Hexact,&dummy);
          break;
        case 'C':
          exactC(grid.year*secpera,r,&Hexact,&dummy);
          break;
        case 'D':
          exactD(grid.year*secpera,r,&Hexact,&dummy);
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
            bothexact(grid.year*secpera,r,&z,1,ApforG,
                      &Hexact,&dummy,&dummy5,&dummy1,&dummy2,&dummy3,&dummy4);
          }
          break;
        case 'H':
          exactH(f,grid.year*secpera,r,&Hexact,&dummy);
          break;
        case 'K':
          Hexact = 3000.0;
          break;
        case 'L':
          Hexact = HexactL[i][j];
          break;
        default:  SETERRQ(1,"test must be A, B, C, D, E, F, G, H, K, or L");
      }

      if (Hexact > 0) {
        areaexact += a;
        volexact += a * Hexact * 1e-3;
      }
      if (i == (grid.Mx - 1)/2 && j == (grid.My - 1)/2) {
        domeH = H[i][j];
        domeHexact = Hexact;
      }
      // compute maximum errors
      Herr = PetscMax(Herr,PetscAbsReal(H[i][j] - Hexact));
      etaerr = PetscMax(etaerr,PetscAbsReal(pow(H[i][j],m) - pow(Hexact,m)));
      // add to sums for average errors
      avHerr += PetscAbsReal(H[i][j] - Hexact);
    }
  }

  ierr = vH.end_access(); CHKERRQ(ierr);
  if (testname == 'L') {
    ierr = vHexactL.end_access(); CHKERRQ(ierr);
  }
  
  // globalize (find errors over all processors) 
  PetscScalar gvol, garea, gdomeH;
  ierr = PetscGlobalSum(&volexact, &gvolexact, grid.com); CHKERRQ(ierr);
  ierr = PetscGlobalMax(&domeHexact, &gdomeHexact, grid.com); CHKERRQ(ierr);
  ierr = PetscGlobalSum(&areaexact, &gareaexact, grid.com); CHKERRQ(ierr);

  ierr = PetscGlobalSum(&vol, &gvol, grid.com); CHKERRQ(ierr);
  ierr = PetscGlobalSum(&area, &garea, grid.com); CHKERRQ(ierr);
  volerr = PetscAbsReal(gvol - gvolexact);
  areaerr = PetscAbsReal(garea - gareaexact);

  ierr = PetscGlobalMax(&Herr, &gmaxHerr, grid.com); CHKERRQ(ierr);
  ierr = PetscGlobalSum(&avHerr, &gavHerr, grid.com); CHKERRQ(ierr);
  gavHerr = gavHerr/(grid.Mx*grid.My);
  ierr = PetscGlobalMax(&etaerr, &gmaxetaerr, grid.com); CHKERRQ(ierr);
  
  ierr = PetscGlobalMax(&domeH, &gdomeH, grid.com); CHKERRQ(ierr);
  centerHerr = PetscAbsReal(gdomeH - gdomeHexact);
  
  return 0;
}


PetscErrorCode IceCompModel::computeBasalVelocityErrors(
      PetscScalar &exactmaxspeed,
      PetscScalar &gmaxvecerr, PetscScalar &gavvecerr,
      PetscScalar &gmaxuberr, PetscScalar &gmaxvberr) {

  PetscErrorCode ierr;
  PetscScalar    **H, **ub, **vb;
  PetscScalar    maxvecerr, avvecerr, maxuberr, maxvberr;
  PetscScalar    ubexact,vbexact, dummy1,dummy2,dummy3;
  
  if (testname != 'E')
    SETERRQ(1,"basal velocity errors only computable for test E\n");

  ierr = vub.get_array(ub); CHKERRQ(ierr);
  ierr = vvb.get_array(vb); CHKERRQ(ierr);
  ierr = vH.get_array(H); CHKERRQ(ierr);
  maxvecerr = 0.0; avvecerr = 0.0; maxuberr = 0.0; maxvberr = 0.0;
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; i++) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; j++) {
      if (H[i][j] > 0.0) {
        PetscScalar r,xx,yy;
        mapcoords(i,j,xx,yy,r);
        exactE(xx,yy,&dummy1,&dummy2,&dummy3,&ubexact,&vbexact); 
        // compute maximum errors
        const PetscScalar uberr = PetscAbsReal(ub[i][j] - ubexact);
        const PetscScalar vberr = PetscAbsReal(vb[i][j] - vbexact);
        maxuberr = PetscMax(maxuberr,uberr);
        maxvberr = PetscMax(maxvberr,vberr);
        const PetscScalar vecerr = sqrt(uberr*uberr + vberr*vberr);
        maxvecerr = PetscMax(maxvecerr,vecerr);
        avvecerr += vecerr;      
      }
    }
  }
  ierr = vH.end_access(); CHKERRQ(ierr);
  ierr = vub.end_access(); CHKERRQ(ierr);
  ierr = vvb.end_access(); CHKERRQ(ierr);

  ierr = PetscGlobalMax(&maxuberr, &gmaxuberr, grid.com); CHKERRQ(ierr);
  ierr = PetscGlobalMax(&maxvberr, &gmaxvberr, grid.com); CHKERRQ(ierr);

  ierr = PetscGlobalMax(&maxvecerr, &gmaxvecerr, grid.com); CHKERRQ(ierr);
  ierr = PetscGlobalSum(&avvecerr, &gavvecerr, grid.com); CHKERRQ(ierr);
  gavvecerr = gavvecerr/(grid.Mx*grid.My);
  
  const PetscScalar xpeak = 450e3 * cos(25.0*(pi/180.0)),
                    ypeak = 450e3 * sin(25.0*(pi/180.0));
  exactE(xpeak,ypeak,&dummy1,&dummy2,&dummy3,&ubexact,&vbexact);
  exactmaxspeed = sqrt(ubexact*ubexact + vbexact*vbexact);
  return 0;
}


PetscErrorCode IceCompModel::additionalAtStartTimestep() {
  PetscErrorCode    ierr;
  
  ierr = verbPrintf(5,grid.com,
               "additionalAtStartTimestep() in IceCompModel entered with test %c",
               testname); CHKERRQ(ierr);

  if (exactOnly == PETSC_TRUE)
    dt_force = config.get("maximum_time_step_years") * secpera;

  // these have no changing boundary conditions or comp sources:
  if (strchr("ABEKL",testname) != NULL) 
    return 0;

  switch (testname) {
    case 'C':
    case 'D':
    case 'H':
      ierr = getCompSourcesTestCDH();
      break;
    case 'F':
    case 'G':
      ierr = getCompSourcesTestFG();  // see iCMthermo.cc
      break;
    default:  SETERRQ(1,"only tests CDHFG have comp source update at start time step\n");
  }

  return 0;
}


PetscErrorCode IceCompModel::additionalAtEndTimestep() {
  PetscErrorCode    ierr;
  
  ierr = verbPrintf(5,grid.com,
               "additionalAtEndTimestep() in IceCompModel entered with test %c",testname);
               CHKERRQ(ierr);

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
      ierr = fillSolnTestABCDH();
      break;
    case 'E':
      ierr = fillSolnTestE();
      break;
    case 'F':
    case 'G':
      ierr = fillSolnTestFG();  // see iCMthermo.cc
      break;
    case 'K':
      ierr = fillSolnTestK();  // see iCMthermo.cc
      break;
    case 'L':
      ierr = fillSolnTestL();
      break;
    default:  SETERRQ(1,"unknown testname in IceCompModel");
  }

  return 0;
}


PetscErrorCode IceCompModel::summary(bool /* tempAndAge */, bool useHomoTemp) {
  //   we always show a summary at every step
  return IceModel::summary(true,useHomoTemp);
}


PetscErrorCode IceCompModel::reportErrors() {
  // geometry errors to report (for all tests): 
  //    -- max thickness error
  //    -- average (at each grid point on whole grid) thickness error
  //    -- max (thickness)^(2n+2)/n error
  //    -- volume error
  //    -- area error
  // and temperature errors (for tests F & G):
  //    -- max T error over 3D domain of ice
  //    -- av T error over 3D domain of ice
  // and basal temperature errors (for tests F & G):
  //    -- max basal temp error
  //    -- average (at each grid point on whole grid) basal temp error
  // and strain-heating (Sigma) errors (for tests F & G):
  //    -- max Sigma error over 3D domain of ice (in 10^-3 K a^-1)
  //    -- av Sigma error over 3D domain of ice (in 10^-3 K a^-1)
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
  ierr = verbPrintf(1,grid.com, 
     "NUMERICAL ERRORS evaluated at final time (relative to exact solution):\n");
  CHKERRQ(ierr);

  int start;
  char filename[TEMPORARY_STRING_LENGTH];
  PetscTruth netcdf_report;
  NCTimeseries err;
  ierr = PetscOptionsGetString(PETSC_NULL, "-report_file", filename,
			       TEMPORARY_STRING_LENGTH, &netcdf_report); CHKERRQ(ierr);

  if (netcdf_report) {
    ierr = verbPrintf(2,grid.com, "Also writing errors to '%s'...\n", filename);
    CHKERRQ(ierr);

    // Find the number of records in this file:
    NCTool nc(grid.com, grid.rank);
    // append = true; check_dims = false
    ierr = nc.open_for_writing(filename); CHKERRQ(ierr);
    ierr = nc.get_dim_length("N", &start); CHKERRQ(ierr);
    ierr = nc.close(); CHKERRQ(ierr);

    ierr = global_attributes.write(filename); CHKERRQ(ierr);

    // Write the dimension variable:
    err.init("N", "N", grid.com, grid.rank);
    ierr = err.write(filename, (size_t)start, (double)(start + 1), NC_INT); CHKERRQ(ierr);

    // Always write grid parameters:
    err.short_name = "dx";
    ierr = err.set_units("meters"); CHKERRQ(ierr);
    ierr = err.write(filename, (size_t)start, grid.dx); CHKERRQ(ierr);
    err.short_name = "dy";
    ierr = err.write(filename, (size_t)start, grid.dy); CHKERRQ(ierr);
    err.short_name = "dz";
    ierr = err.write(filename, (size_t)start, grid.dzMAX); CHKERRQ(ierr);
    err.short_name = "dzb";
    ierr = err.write(filename, (size_t)start, grid.dzbMAX); CHKERRQ(ierr);

    // Always write the test name:
    err.reset();
    err.short_name = "test";
    ierr = err.write(filename, (size_t)start, (double)testname, NC_BYTE); CHKERRQ(ierr);
  }
  // geometry (thickness, vol) errors if appropriate; reported in m except for relmaxETA
  if (testname != 'K') {
    PetscScalar volexact, areaexact, domeHexact, volerr, areaerr, maxHerr, avHerr,
                maxetaerr, centerHerr;
    ierr = computeGeometryErrors(volexact,areaexact,domeHexact,
                                 volerr,areaerr,maxHerr,avHerr,maxetaerr,centerHerr);
            CHKERRQ(ierr);
    ierr = verbPrintf(1,grid.com, 
            "geometry  :    prcntVOL        maxH         avH   relmaxETA\n");
            CHKERRQ(ierr);  // no longer reporting centerHerr
    const PetscScalar   m = (2.0 * tgaIce->exponent() + 2.0) / tgaIce->exponent();
    ierr = verbPrintf(1,grid.com, "           %12.6f%12.6f%12.6f%12.6f\n",
                      100*volerr/volexact, maxHerr, avHerr,
                      maxetaerr/pow(domeHexact,m)); CHKERRQ(ierr);

    if (netcdf_report) {
      err.reset();
      err.short_name = "relative_volume";
      ierr = err.set_units("percent"); CHKERRQ(ierr);
      err.set_string("long_name", "relative ice volume error");
      ierr = err.write(filename, (size_t)start, 100*volerr/volexact); CHKERRQ(ierr);

      err.short_name = "relative_max_eta";
      ierr = err.set_units("1"); CHKERRQ(ierr);
      err.set_string("long_name", "relative $\\eta$ error");
      ierr = err.write(filename, (size_t)start, maxetaerr/pow(domeHexact,m)); CHKERRQ(ierr);

      err.short_name = "maximum_thickness";
      ierr = err.set_units("meters"); CHKERRQ(ierr);
      err.set_string("long_name", "maximum ice thickness error");
      ierr = err.write(filename, (size_t)start, maxHerr); CHKERRQ(ierr);

      err.short_name = "average_thickness";
      ierr = err.set_units("meters"); CHKERRQ(ierr);
      err.set_string("long_name", "average ice thickness error");
      ierr = err.write(filename, (size_t)start, avHerr); CHKERRQ(ierr);
    }
  }

  // temperature errors if appropriate; reported in K
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
      err.reset();
      err.short_name = "maximum_temperature";
      ierr = err.set_units("Kelvin"); CHKERRQ(ierr);
      err.set_string("long_name", "maximum ice temperature error");
      ierr = err.write(filename, (size_t)start, maxTerr); CHKERRQ(ierr);

      err.short_name = "average_temperature";
      err.set_string("long_name", "average ice temperature error");
      ierr = err.write(filename, (size_t)start, avTerr); CHKERRQ(ierr);

      err.short_name = "maximum_basal_temperature";
      err.set_string("long_name", "maximum basal temperature error");
      ierr = err.write(filename, (size_t)start, basemaxTerr); CHKERRQ(ierr);
      err.short_name = "average_basal_temperature";
      err.set_string("long_name", "average basal temperature error");
      ierr = err.write(filename, (size_t)start, baseavTerr); CHKERRQ(ierr);
    }

  } else if (testname == 'K') {
    PetscScalar maxTerr, avTerr, maxTberr, avTberr;
    ierr = computeIceBedrockTemperatureErrors(maxTerr, avTerr, maxTberr, avTberr);
       CHKERRQ(ierr);
    ierr = verbPrintf(1,grid.com, 
       "temp      :        maxT         avT       maxTb        avTb\n"); CHKERRQ(ierr);
    ierr = verbPrintf(1,grid.com, "           %12.6f%12.6f%12.6f%12.6f\n", 
                  maxTerr, avTerr, maxTberr, avTberr); CHKERRQ(ierr);

    if (netcdf_report) {
      err.reset();
      err.short_name = "maximum_temperature";
      ierr = err.set_units("Kelvin"); CHKERRQ(ierr);
      err.set_string("long_name", "maximum ice temperature error");
      ierr = err.write(filename, (size_t)start, maxTerr); CHKERRQ(ierr);

      err.short_name = "average_temperature";
      err.set_string("long_name", "average ice temperature error");
      ierr = err.write(filename, (size_t)start, avTerr); CHKERRQ(ierr);

      err.short_name = "maximum_bedrock_temperature";
      err.set_string("long_name", "maximum bedrock temperature error");
      ierr = err.write(filename, (size_t)start, maxTberr); CHKERRQ(ierr);

      err.short_name = "average_bedrock_temperature";
      err.set_string("long_name", "average bedrock temperature error");
      ierr = err.write(filename, (size_t)start, avTberr); CHKERRQ(ierr);
    }
  }

  // Sigma errors if appropriate; reported in 10^6 J/(s m^3)
  if ((testname == 'F') || (testname == 'G')) {
    PetscScalar maxSigerr, avSigerr;
    ierr = computeSigmaErrors(maxSigerr, avSigerr); CHKERRQ(ierr);
    ierr = verbPrintf(1,grid.com, 
       "Sigma     :      maxSig       avSig\n"); CHKERRQ(ierr);
    ierr = verbPrintf(1,grid.com, "           %12.6f%12.6f\n", 
                  maxSigerr*1.0e6, avSigerr*1.0e6); CHKERRQ(ierr);

    if (netcdf_report) {
      err.reset();
      err.short_name = "maximum_sigma";
      ierr = err.set_units("J s-1 m-3"); CHKERRQ(ierr);
      ierr = err.set_glaciological_units("1e6 J s-1 m-3"); CHKERRQ(ierr);
      err.set_string("long_name", "maximum strain heating error");
      ierr = err.write(filename, (size_t)start, maxSigerr); CHKERRQ(ierr);

      err.short_name = "average_sigma";
      err.set_string("long_name", "average strain heating error");
      ierr = err.write(filename, (size_t)start, avSigerr); CHKERRQ(ierr);
    }
  }

  // surface velocity errors if exact values are available; reported in m/a
  if ((testname == 'F') || (testname == 'G')) {
    PetscScalar maxUerr, avUerr, maxWerr, avWerr;
    ierr = computeSurfaceVelocityErrors(maxUerr, avUerr, maxWerr, avWerr); CHKERRQ(ierr);
    ierr = verbPrintf(1,grid.com, 
       "surf vels :     maxUvec      avUvec        maxW         avW\n"); CHKERRQ(ierr);
    ierr = verbPrintf(1,grid.com, "           %12.6f%12.6f%12.6f%12.6f\n", 
                  maxUerr*secpera, avUerr*secpera, maxWerr*secpera, avWerr*secpera); CHKERRQ(ierr);

    if (netcdf_report) {
      err.reset();
      err.short_name = "maximum_surface_velocity";
      err.set_string("long_name", "maximum ice surface horizontal velocity error");
      ierr = err.set_units("m/s"); CHKERRQ(ierr);
      ierr = err.set_glaciological_units("meters/year"); CHKERRQ(ierr);
      ierr = err.write(filename, (size_t)start, maxUerr); CHKERRQ(ierr);

      err.short_name = "average_surface_velocity";
      err.set_string("long_name", "average ice surface horizontal velocity error");
      ierr = err.write(filename, (size_t)start, avUerr); CHKERRQ(ierr);

      err.short_name = "maximum_surface_w";
      err.set_string("long_name", "maximum ice surface vertical velocity error");
      ierr = err.write(filename, (size_t)start, maxWerr); CHKERRQ(ierr);

      err.short_name = "average_surface_w";
      err.set_string("long_name", "average ice surface vertical velocity error");
      ierr = err.write(filename, (size_t)start, avWerr); CHKERRQ(ierr);
    }
  }

  // basal velocity errors if appropriate; reported in m/a except prcntavvec
  if (testname == 'E') {
    PetscScalar exactmaxspeed, maxvecerr, avvecerr, maxuberr, maxvberr;
    ierr = computeBasalVelocityErrors(exactmaxspeed,
                          maxvecerr,avvecerr,maxuberr,maxvberr); CHKERRQ(ierr);
    ierr = verbPrintf(1,grid.com, 
       "base vels :  maxvector   avvector  prcntavvec     maxub     maxvb\n");
       CHKERRQ(ierr);
    ierr = verbPrintf(1,grid.com, "           %11.4f%11.5f%12.5f%10.4f%10.4f\n", 
                  maxvecerr*secpera, avvecerr*secpera, 
                  (avvecerr/exactmaxspeed)*100.0,
                  maxuberr*secpera, maxvberr*secpera); CHKERRQ(ierr);

    if (netcdf_report) {
      err.reset();
      err.short_name = "maximum_basal_velocity";
      ierr = err.set_units("m/s"); CHKERRQ(ierr);
      ierr = err.set_glaciological_units("meters/year"); CHKERRQ(ierr);
      ierr = err.write(filename, (size_t)start, maxvecerr); CHKERRQ(ierr);

      err.short_name = "average_basal_velocity";
      ierr = err.write(filename, (size_t)start, avvecerr); CHKERRQ(ierr);
      err.short_name = "maximum_basal_u";
      ierr = err.write(filename, (size_t)start, maxuberr); CHKERRQ(ierr);
      err.short_name = "maximum_basal_v";
      ierr = err.write(filename, (size_t)start, maxvberr); CHKERRQ(ierr);

      err.reset();
      err.short_name = "relative_basal_velocity";
      ierr = err.set_units("percent"); CHKERRQ(ierr);
      ierr = err.write(filename, (size_t)start, (avvecerr/exactmaxspeed)*100); CHKERRQ(ierr);
    }
  }

  ierr = verbPrintf(1,grid.com, "NUM ERRORS DONE\n");  CHKERRQ(ierr);
  return 0;
}

