// Copyright (C) 2004-2016 Jed Brown, Ed Bueler and Constantine Khroulev
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
#include <cassert>

#include <vector>     // STL vector container; sortable; used in test L
#include <algorithm>  // required by sort(...) in test L

#include "tests/exactTestsABCD.h"
#include "tests/exactTestsFG.h"
#include "tests/exactTestH.h"
#include "tests/exactTestL.h"

#include "iceCompModel.hh"
#include "base/stressbalance/sia/SIAFD.hh"
#include "base/stressbalance/ShallowStressBalance.hh"
#include "base/rheology/PatersonBuddCold.hh"
#include "base/stressbalance/PISMStressBalance.hh"
#include "base/enthalpyConverter.hh"
#include "base/util/io/PIO.hh"
#include "base/util/pism_options.hh"
#include "coupler/ocean/POConstant.hh"
#include "PSVerification.hh"
#include "base/util/Mask.hh"
#include "base/util/error_handling.hh"
#include "earth/PISMBedDef.hh"
#include "base/util/IceGrid.hh"
#include "base/util/PISMTime.hh"
#include "base/util/PISMConfigInterface.hh"
#include "base/util/Context.hh"
#include "base/util/io/io_helpers.hh"
#include "base/util/Logger.hh"
#include "base/util/pism_utilities.hh"

namespace pism {

const double IceCompModel::secpera = 3.15569259747e7;

using units::convert;

IceCompModel::IceCompModel(IceGrid::Ptr g, Context::Ptr context, int mytest)
  : IceModel(g, context) {

  // note lots of defaults are set by the IceModel constructor

  // defaults for IceCompModel:
  testname = mytest;
  exactOnly = false;
  bedrock_is_ice_forK = false;

  // Override some defaults from parent class
  m_config->set_double("sia_enhancement_factor", 1.0);
  // none use bed smoothing & bed roughness parameterization
  m_config->set_double("bed_smoother_range", 0.0);

  // set values of flags in run()
  m_config->set_boolean("do_mass_conserve", true);
  m_config->set_boolean("include_bmr_in_continuity", false);

  if (testname == 'V') {
    m_config->set_string("ssa_flow_law", "isothermal_glen");
    m_config->set_double("ice_softness", pow(1.9e8, -m_config->get_double("sia_Glen_exponent")));
  } else {
    // Set the default for IceCompModel:
    m_config->set_string("sia_flow_law", "arr");
  }
}

void IceCompModel::createVecs() {

  IceModel::createVecs();

  vHexactL.create(m_grid, "HexactL", WITH_GHOSTS, 2);

  strain_heating3_comp.create(m_grid,"strain_heating_comp", WITHOUT_GHOSTS);
  strain_heating3_comp.set_attrs("internal","rate of compensatory strain heating in ice",
                                 "W m-3", "");
}

void IceCompModel::setFromOptions() {

  m_log->message(2, "starting Test %c ...\n", testname);

  /* This switch turns off actual numerical evolution and simply reports the
     exact solution. */
  bool flag = options::Bool("-eo", "exact only");
  if (flag) {
    exactOnly = true;
    m_log->message(1, "!!EXACT SOLUTION ONLY, NO NUMERICAL SOLUTION!!\n");
  }

  // These ifs are here (and not in the constructor or later) because
  // testname actually comes from a command-line *and* because command-line
  // options should be able to override parameter values set here.

  if (testname == 'H') {
    m_config->set_string("bed_deformation_model", "iso");
  } else
    m_config->set_string("bed_deformation_model", "none");

  if ((testname == 'F') || (testname == 'G') || (testname == 'K') || (testname == 'O')) {
    m_config->set_boolean("do_energy", true);
    // essentially turn off run-time reporting of extremely low computed
    // temperatures; *they will be reported as errors* anyway
    m_config->set_double("global_min_allowed_temp", 0.0);
    m_config->set_double("max_low_temp_count", 1000000);
  } else {
    m_config->set_boolean("do_energy", false);
  }

  m_config->set_boolean("is_dry_simulation", true);

  // special considerations for K and O wrt thermal bedrock and pressure-melting
  if ((testname == 'K') || (testname == 'O')) {
    m_config->set_boolean("temperature_allow_above_melting", false);
  } else {
    // note temps are generally allowed to go above pressure melting in verify
    m_config->set_boolean("temperature_allow_above_melting", true);
  }

  if (testname == 'V') {
    // no sub-shelf melting
    m_config->set_boolean("include_bmr_in_continuity", false);

    // this test is isothermal
    m_config->set_boolean("do_energy", false);

    // use the SSA solver
    m_config->set_string("stress_balance_model", "ssa");

    // this certainly is not a "dry simulation"
    m_config->set_boolean("is_dry_simulation", false);

    m_config->set_boolean("ssa_dirichlet_bc", true);
  }

  m_config->set_boolean("do_cold_ice_methods", true);

  IceModel::setFromOptions();
}

void IceCompModel::allocate_bedrock_thermal_unit() {

  if (btu != NULL) {
    return;
  }

  // this switch changes Test K to make material properties for bedrock the same as for ice
  bool biiSet = options::Bool("-bedrock_is_ice", "set bedrock properties to those of ice");
  if (biiSet == true) {
    if (testname == 'K') {
      m_log->message(1,
                 "setting material properties of bedrock to those of ice in Test K\n");
      m_config->set_double("bedrock_thermal_density", m_config->get_double("ice_density"));
      m_config->set_double("bedrock_thermal_conductivity", m_config->get_double("ice_thermal_conductivity"));
      m_config->set_double("bedrock_thermal_specific_heat_capacity", m_config->get_double("ice_specific_heat_capacity"));
      bedrock_is_ice_forK = true;
    } else {
      m_log->message(1,
                 "IceCompModel WARNING: option -bedrock_is_ice ignored; only applies to Test K\n");
    }
  }

  if (testname != 'K') {
    // now make bedrock have same material properties as ice
    // (note Mbz=1 also, by default, but want ice/rock interface to see
    // pure ice from the point of view of applying geothermal boundary
    // condition, especially in tests F and G)
    m_config->set_double("bedrock_thermal_density", m_config->get_double("ice_density"));
    m_config->set_double("bedrock_thermal_conductivity", m_config->get_double("ice_thermal_conductivity"));
    m_config->set_double("bedrock_thermal_specific_heat_capacity", m_config->get_double("ice_specific_heat_capacity"));
  }

  btu = new energy::BTU_Verification(m_grid, testname, bedrock_is_ice_forK);
}

void IceCompModel::allocate_stressbalance() {

  using namespace pism::stressbalance;

  if (m_stress_balance != NULL) {
    return;
  }

  EnthalpyConverter::Ptr EC = m_ctx->enthalpy_converter();

  IceModel::allocate_stressbalance();

  if (testname != 'V') {
    // check on whether the options (already checked) chose the right
    // FlowLaw for verification (we need to have the right flow law for
    // errors to make sense)

    rheology::FlowLaw *ice = m_stress_balance->get_ssb_modifier()->flow_law();

    if (not FlowLawIsPatersonBuddCold(ice, *m_config, EC)) {
      m_log->message(1,
                 "WARNING: SIA flow law should be '-sia_flow_law arr' for the selected pismv test.\n");
    }
  }
}

void IceCompModel::allocate_bed_deformation() {

  IceModel::allocate_bed_deformation();

  f = m_config->get_double("ice_density") / m_config->get_double("lithosphere_density");  // for simple isostasy

  std::string bed_def_model = m_config->get_string("bed_deformation_model");

  if ((testname == 'H') && bed_def_model != "iso") {
    m_log->message(1,
               "IceCompModel WARNING: Test H should be run with option\n"
               "  '-bed_def iso'  for the reported errors to be correct.\n");
  }
}

void IceCompModel::allocate_couplers() {
  EnthalpyConverter::Ptr EC = m_ctx->enthalpy_converter();

  // Climate will always come from verification test formulas.
  m_surface = new surface::Verification(m_grid, EC, testname);
  m_ocean   = new ocean::Constant(m_grid);
}

void IceCompModel::set_vars_from_options() {

  // -bootstrap command-line option is not allowed here.
  options::forbidden("-bootstrap");

  strain_heating3_comp.set(0.0);

  m_log->message(3,
             "initializing Test %c from formulas ...\n",testname);

  // all have no uplift
  IceModelVec2S bed_uplift;
  bed_uplift.create(m_grid, "uplift", WITHOUT_GHOSTS);
  bed_uplift.set(0);
  m_beddef->set_uplift(bed_uplift);

  // this is the correct initialization for Test O (and every other
  // test; they all generate zero basal melt rate)
  m_basal_melt_rate.set(0.0);

  // Test-specific initialization:
  switch (testname) {
  case 'A':
  case 'B':
  case 'C':
  case 'D':
  case 'H':
    initTestABCDH();
    break;
  case 'F':
  case 'G':
    initTestFG();  // see iCMthermo.cc
    break;
  case 'K':
  case 'O':
    initTestsKO();  // see iCMthermo.cc
    break;
  case 'L':
    initTestL();
    break;
  case 'V':
    test_V_init();
    break;
  default:
    throw RuntimeError("Desired test not implemented by IceCompModel.");
  }

  compute_enthalpy_cold(m_ice_temperature, m_ice_enthalpy);
}

void IceCompModel::initTestABCDH() {
  double     A0, T0, H, accum;

  EnthalpyConverter::Ptr EC = m_ctx->enthalpy_converter();

  rheology::PatersonBuddCold tgaIce("sia_", *m_config, EC);

  const double time = m_time->current();

  // compute T so that A0 = A(T) = Acold exp(-Qcold/(R T))  (i.e. for PatersonBuddCold);
  // set all temps to this constant
  A0 = 1.0e-16/secpera;    // = 3.17e-24  1/(Pa^3 s);  (EISMINT value) flow law parameter
  T0 = tgaIce.tempFromSoftness(A0);

  m_ice_temperature.set(T0);
  m_geothermal_flux.set(Ggeo);
  m_cell_type.set(MASK_GROUNDED);

  IceModelVec::AccessList list(m_ice_thickness);

  ParallelSection loop(m_grid->com);
  try {
    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      const double r = radius(*m_grid, i, j);
      switch (testname) {
      case 'A':
        exactA(r, &H, &accum);
        m_ice_thickness(i, j)   = H;
        break;
      case 'B':
        exactB(time, r, &H, &accum);
        m_ice_thickness(i, j)   = H;
        break;
      case 'C':
        exactC(time, r, &H, &accum);
        m_ice_thickness(i, j)   = H;
        break;
      case 'D':
        exactD(time, r, &H, &accum);
        m_ice_thickness(i, j)   = H;
        break;
      case 'H':
        exactH(f, time, r, &H, &accum);
        m_ice_thickness(i, j)   = H;
        break;
      default:
        throw RuntimeError("test must be A, B, C, D, or H");
      }
    }
  } catch (...) {
    loop.failed();
  }
  loop.check();

  m_ice_thickness.update_ghosts();

  {
    IceModelVec2S bed_topography;
    bed_topography.create(m_grid, "topg", WITHOUT_GHOSTS);

    if (testname == 'H') {
      bed_topography.copy_from(m_ice_thickness);
      bed_topography.scale(-f);
    } else {  // flat bed case otherwise
      bed_topography.set(0.0);
    }
    m_beddef->set_elevation(bed_topography);
  }
}

//! Class used initTestL() in generating sorted list for ODE solver.
class rgrid {
public:
  double r;
  int    i,j;
};

//! Comparison used initTestL() in generating sorted list for ODE solver.
struct rgridReverseSort {
  bool operator()(rgrid a, rgrid b) {
    return (a.r > b.r);
  }
};

void IceCompModel::initTestL() {
  int ierr;
  double     A0, T0;

  EnthalpyConverter::Ptr EC = m_ctx->enthalpy_converter();

  assert(testname == 'L');

  rheology::PatersonBuddCold tgaIce("sia_", *m_config, EC);

  // compute T so that A0 = A(T) = Acold exp(-Qcold/(R T))  (i.e. for PatersonBuddCold);
  // set all temps to this constant
  A0 = 1.0e-16/secpera;    // = 3.17e-24  1/(Pa^3 s);  (EISMINT value) flow law parameter
  T0 = tgaIce.tempFromSoftness(A0);

  m_ice_temperature.set(T0);
  m_geothermal_flux.set(Ggeo);

  // setup to evaluate test L; requires solving an ODE numerically
  //   using sorted list of radii, sorted in decreasing radius order
  const int MM = m_grid->xm() * m_grid->ym();

  std::vector<rgrid> rrv(MM);  // destructor at end of scope
  int k = 0;
  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    rrv[k].i = i;
    rrv[k].j = j;
    rrv[k].r = radius(*m_grid, i,j);

    k += 1;
  }

  std::sort(rrv.begin(), rrv.end(), rgridReverseSort()); // so rrv[k].r > rrv[k+1].r

  // get soln to test L at these radii; solves ODE only once (on each processor)
  std::vector<double> rr(MM), HH(MM), bb(MM), aa(MM);

  for (k = 0; k < MM; k++) {
    rr[k] = rrv[k].r;
  }

  ierr = exactL_list(&rr[0], MM, &HH[0], &bb[0], &aa[0]);
  switch (ierr) {
     case TESTL_NOT_DONE:
       m_log->message(1,
          "\n\nTest L ERROR: exactL_list() returns 'NOT_DONE' (%d) ...\n\n\n",ierr);
       break;
     case TESTL_NOT_DECREASING:
       m_log->message(1,
          "\n\nTest L ERROR: exactL_list() returns 'NOT_DECREASING' (%d) ...\n\n\n",ierr);
       break;
     case TESTL_INVALID_METHOD:
       m_log->message(1,
          "\n\nTest L ERROR: exactL_list() returns 'INVALID_METHOD' (%d) ...\n\n\n",ierr);
       break;
     case TESTL_NO_LIST:
       m_log->message(1,
          "\n\nTest L ERROR: exactL_list() returns 'NO_LIST' (%d) ...\n\n\n",ierr);
       break;
     default:
       break;
  }
  if (ierr != 0) {
    throw RuntimeError("test L: exactL_list(..) failed");
  }

  {
    IceModelVec2S bed_topography;
    bed_topography.create(m_grid, "topg", WITHOUT_GHOSTS);

    IceModelVec::AccessList list;
    list.add(m_ice_thickness);
    list.add(bed_topography);

    for (k = 0; k < MM; k++) {
      m_ice_thickness(rrv[k].i, rrv[k].j)  = HH[k];
      bed_topography(rrv[k].i, rrv[k].j) = bb[k];
    }

    m_ice_thickness.update_ghosts();
    m_beddef->set_elevation(bed_topography);
  }

  // store copy of ice_thickness for "-eo" runs and for evaluating geometry errors
  vHexactL.copy_from(m_ice_thickness);
}

//! \brief Tests A and E have a thickness B.C. (ice_thickness == 0 outside a circle of radius 750km).
void IceCompModel::reset_thickness_test_A() {
  const double LforAE = 750e3; // m

  IceModelVec::AccessList list(m_ice_thickness);

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    if (radius(*m_grid, i, j) > LforAE) {
      m_ice_thickness(i, j) = 0;
    }
  }

  m_ice_thickness.update_ghosts();
}



void IceCompModel::fillSolnTestABCDH() {
  double     H, accum;

  const double time = m_time->current();

  IceModelVec::AccessList list(m_ice_thickness);

  ParallelSection loop(m_grid->com);
  try {
    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      double r = radius(*m_grid, i, j);
      switch (testname) {
      case 'A':
        exactA(r, &H, &accum);
        m_ice_thickness(i, j)   = H;
        break;
      case 'B':
        exactB(time, r, &H, &accum);
        m_ice_thickness(i, j)   = H;
        break;
      case 'C':
        exactC(time, r, &H, &accum);
        m_ice_thickness(i, j)   = H;
        break;
      case 'D':
        exactD(time, r, &H, &accum);
        m_ice_thickness(i, j)   = H;
        break;
      case 'H':
        exactH(f, time, r, &H, &accum);
        m_ice_thickness(i, j)   = H;
        break;
      default:
        throw RuntimeError("test must be A, B, C, D, or H");
      }
    }
  } catch (...) {
    loop.failed();
  }
  loop.check();

  m_ice_thickness.update_ghosts();

  {
    IceModelVec2S bed_topography;
    bed_topography.create(m_grid, "topg", WITHOUT_GHOSTS);

    if (testname == 'H') {
      bed_topography.copy_from(m_ice_thickness);
      bed_topography.scale(-f);
    } else {
      bed_topography.set(0.0);
    }
    m_beddef->set_elevation(bed_topography);
  }
}


void IceCompModel::fillSolnTestL() {

  vHexactL.update_ghosts();
  m_ice_thickness.copy_from(vHexactL);

  // note bed was filled at initialization and hasn't changed
}


void IceCompModel::computeGeometryErrors(double &gvolexact, double &gareaexact,
                                         double &gdomeHexact, double &volerr,
                                         double &areaerr, double &gmaxHerr,
                                         double &gavHerr, double &gmaxetaerr,
                                         double &centerHerr) {
  // compute errors in thickness, eta=thickness^{(2n+2)/n}, volume, area

  const double time = m_time->current();
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

  double     dummy, z, dummy1, dummy2, dummy3, dummy4, dummy5;

  IceModelVec::AccessList list(m_ice_thickness);
  if (testname == 'L') {
    list.add(vHexactL);
  }

  double
    seawater_density = m_config->get_double("sea_water_density"),
    ice_density      = m_config->get_double("ice_density"),
    Glen_n           = m_config->get_double("sia_Glen_exponent"),
    standard_gravity = m_config->get_double("standard_gravity");

  // area of grid square in square km:
  const double   a = m_grid->dx() * m_grid->dy() * 1e-3 * 1e-3;
  const double   m = (2.0 * Glen_n + 2.0) / Glen_n;

  ParallelSection loop(m_grid->com);
  try {
    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      if (m_ice_thickness(i,j) > 0) {
        area += a;
        vol += a * m_ice_thickness(i,j) * 1e-3;
      }
      double xx = m_grid->x(i), r = radius(*m_grid, i,j);
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
      case 'F':
        if (r > LforFG - 1.0) {  // outside of sheet
          Hexact=0.0;
        } else {
          r=std::max(r,1.0);
          z=0.0;
          bothexact(0.0,r,&z,1,0.0,
                    &Hexact,&dummy,&dummy5,&dummy1,&dummy2,&dummy3,&dummy4);
        }
        break;
      case 'G':
        if (r > LforFG -1.0) {  // outside of sheet
          Hexact=0.0;
        } else {
          r=std::max(r,1.0);
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
          double
            H0 = 600.0,
            v0 = convert(m_sys, 300.0, "m year-1", "m second-1"),
            Q0 = H0 * v0,
            B0 = m_stress_balance->get_stressbalance()->flow_law()->hardness(0, 0),
            C  = pow(ice_density * standard_gravity * (1.0 - ice_density/seawater_density) / (4 * B0), 3);

          Hexact = pow(4 * C / Q0 * xx + 1/pow(H0, 4), -0.25);
        }
        break;
      default:
        throw RuntimeError("test must be A, B, C, D, F, G, H, K, L, or O");
      }

      if (Hexact > 0) {
        areaexact += a;
        volexact += a * Hexact * 1e-3;
      }
      if (i == ((int)m_grid->Mx() - 1)/2 and
          j == ((int)m_grid->My() - 1)/2) {
        domeH = m_ice_thickness(i,j);
        domeHexact = Hexact;
      }
      // compute maximum errors
      Herr = std::max(Herr,fabs(m_ice_thickness(i,j) - Hexact));
      etaerr = std::max(etaerr,fabs(pow(m_ice_thickness(i,j),m) - pow(Hexact,m)));
      // add to sums for average errors
      avHerr += fabs(m_ice_thickness(i,j) - Hexact);
    }
  } catch (...) {
    loop.failed();
  }
  loop.check();

  // globalize (find errors over all processors)
  double gvol, garea, gdomeH;
  gvolexact = GlobalSum(m_grid->com, volexact);
  gdomeHexact = GlobalMax(m_grid->com, domeHexact);
  gareaexact = GlobalSum(m_grid->com, areaexact);

  gvol = GlobalSum(m_grid->com, vol);
  garea = GlobalSum(m_grid->com, area);
  volerr = fabs(gvol - gvolexact);
  areaerr = fabs(garea - gareaexact);

  gmaxHerr = GlobalMax(m_grid->com, Herr);
  gavHerr = GlobalSum(m_grid->com, avHerr);
  gavHerr = gavHerr/(m_grid->Mx()*m_grid->My());
  gmaxetaerr = GlobalMax(m_grid->com, etaerr);

  gdomeH = GlobalMax(m_grid->com, domeH);
  centerHerr = fabs(gdomeH - gdomeHexact);
}


void IceCompModel::additionalAtStartTimestep() {
  if (testname == 'F' || testname == 'G') {
    getCompSourcesTestFG();
  }
}


void IceCompModel::additionalAtEndTimestep() {

  if (testname == 'A') {
    reset_thickness_test_A();
  }

  // do nothing at the end of the time step unless the user has asked for the
  // exact solution to overwrite the numerical solution
  if (not exactOnly) {
    return;
  }

  // because user wants exact solution, fill gridded values from exact formulas;
  // important notes:
  //     (1) the numerical computation *has* already occurred, in run(),
  //           and we just overwrite it with the exact solution here
  //     (2) certain diagnostic quantities like dHdt are computed numerically,
  //           and not overwritten here; while velbar_mag,velsurf_mag,flux_mag,wsurf are diagnostic
  //           quantities recomputed at the end of the run for writing into
  //           NetCDF, in particular dHdt is not recomputed before being written
  //           into the output file, so it is actually numerical
  switch (testname) {
  case 'A':
  case 'B':
  case 'C':
  case 'D':
  case 'H':
    fillSolnTestABCDH();
    break;
  case 'F':
  case 'G':
    fillSolnTestFG(); // see iCMthermo.cc
    break;
  case 'K':
    fillTemperatureSolnTestsKO(); // see iCMthermo.cc
    break;
  case 'O':
    fillTemperatureSolnTestsKO(); // see iCMthermo.cc
    fillBasalMeltRateSolnTestO(); // see iCMthermo.cc
    break;
  case 'L':
    fillSolnTestL();
    break;
  default:
    throw RuntimeError::formatted("unknown testname %c in IceCompModel", testname);
  }
}


void IceCompModel::summary(bool /* tempAndAge */) {
  //   we always show a summary at every step
  IceModel::summary(true);
}


void IceCompModel::reportErrors() {
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

  bool no_report = options::Bool("-no_report", "Don't report numerical errors");

  if (no_report) {
    return;
  }

  EnthalpyConverter::Ptr EC = m_ctx->enthalpy_converter();

  rheology::FlowLaw* flow_law = m_stress_balance->get_ssb_modifier()->flow_law();
  if ((testname == 'F' or testname == 'G') and
      testname != 'V' and
      not FlowLawIsPatersonBuddCold(flow_law, *m_config, EC)) {
    m_log->message(1,
               "pismv WARNING: flow law must be cold part of Paterson-Budd ('-siafd_flow_law arr')\n"
               "   for reported errors in test %c to be meaningful!\n",
               testname);
  }

  m_log->message(1,
             "NUMERICAL ERRORS evaluated at final time (relative to exact solution):\n");

  unsigned int start;
  TimeseriesMetadata err("N", "N", m_sys);

  err.set_string("units", "1");

  PIO nc(m_grid->com, "netcdf3"); // OK to use netcdf3

  options::String report_file("-report_file", "NetCDF error report file");
  bool append = options::Bool("-append", "Append the NetCDF error report");

  IO_Mode mode = append ? PISM_READWRITE : PISM_READWRITE_MOVE;

  if (report_file.is_set()) {
    m_log->message(2, "Also writing errors to '%s'...\n", report_file->c_str());

    // Find the number of records in this file:
    nc.open(report_file, mode);
    start = nc.inq_dimlen("N");

    io::write_global_attributes(nc, m_output_global_attributes);

    // Write the dimension variable:
    io::write_timeseries(nc, err, (size_t)start, (double)(start + 1), PISM_INT);

    // Always write grid parameters:
    err.set_name("dx");
    err.set_string("units", "meters");
    io::write_timeseries(nc, err, (size_t)start, m_grid->dx());
    err.set_name("dy");
    io::write_timeseries(nc, err, (size_t)start, m_grid->dy());
    err.set_name("dz");
    io::write_timeseries(nc, err, (size_t)start, m_grid->dz_max());

    // Always write the test name:
    err.clear_all_strings(); err.clear_all_doubles(); err.set_string("units", "1");
    err.set_name("test");
    io::write_timeseries(nc, err, (size_t)start, (double)testname, PISM_BYTE);
  }

  // geometry (thickness, vol) errors if appropriate; reported in m except for relmaxETA
  if ((testname != 'K') && (testname != 'O')) {
    double volexact, areaexact, domeHexact, volerr, areaerr, maxHerr, avHerr,
                maxetaerr, centerHerr;
    computeGeometryErrors(volexact,areaexact,domeHexact,
                          volerr,areaerr,maxHerr,avHerr,maxetaerr,centerHerr);
    m_log->message(1,
               "geometry  :    prcntVOL        maxH         avH   relmaxETA\n");  // no longer reporting centerHerr
    const double   m = (2.0 * flow_law->exponent() + 2.0) / flow_law->exponent();
    m_log->message(1, "           %12.6f%12.6f%12.6f%12.6f\n",
               100*volerr/volexact, maxHerr, avHerr,
               maxetaerr/pow(domeHexact,m));

    if (report_file.is_set()) {
      err.clear_all_strings(); err.clear_all_doubles(); err.set_string("units", "1");
      err.set_name("relative_volume");
      err.set_string("units", "percent");
      err.set_string("long_name", "relative ice volume error");
      io::write_timeseries(nc, err, (size_t)start, 100*volerr/volexact);

      err.set_name("relative_max_eta");
      err.set_string("units", "1");
      err.set_string("long_name", "relative $\\eta$ error");
      io::write_timeseries(nc, err, (size_t)start, maxetaerr/pow(domeHexact,m));

      err.set_name("maximum_thickness");
      err.set_string("units", "meters");
      err.set_string("long_name", "maximum ice thickness error");
      io::write_timeseries(nc, err, (size_t)start, maxHerr);

      err.set_name("average_thickness");
      err.set_string("units", "meters");
      err.set_string("long_name", "average ice thickness error");
      io::write_timeseries(nc, err, (size_t)start, avHerr);
    }
  }

  // temperature errors for F and G
  if ((testname == 'F') || (testname == 'G')) {
    double maxTerr, avTerr, basemaxTerr, baseavTerr, basecenterTerr;
    computeTemperatureErrors(maxTerr, avTerr);
    computeBasalTemperatureErrors(basemaxTerr, baseavTerr, basecenterTerr);
    m_log->message(1,
               "temp      :        maxT         avT    basemaxT     baseavT\n");  // no longer reporting   basecenterT
    m_log->message(1, "           %12.6f%12.6f%12.6f%12.6f\n",
               maxTerr, avTerr, basemaxTerr, baseavTerr);

    if (report_file.is_set()) {
      err.clear_all_strings(); err.clear_all_doubles(); err.set_string("units", "1");
      err.set_name("maximum_temperature");
      err.set_string("units", "Kelvin");
      err.set_string("long_name", "maximum ice temperature error");
      io::write_timeseries(nc, err, (size_t)start, maxTerr);

      err.set_name("average_temperature");
      err.set_string("long_name", "average ice temperature error");
      io::write_timeseries(nc, err, (size_t)start, avTerr);

      err.set_name("maximum_basal_temperature");
      err.set_string("long_name", "maximum basal temperature error");
      io::write_timeseries(nc, err, (size_t)start, basemaxTerr);
      err.set_name("average_basal_temperature");
      err.set_string("long_name", "average basal temperature error");
      io::write_timeseries(nc, err, (size_t)start, baseavTerr);
    }

  } else if ((testname == 'K') || (testname == 'O')) {
    double maxTerr, avTerr, maxTberr, avTberr;
    computeIceBedrockTemperatureErrors(maxTerr, avTerr, maxTberr, avTberr);
    m_log->message(1,
               "temp      :        maxT         avT       maxTb        avTb\n");
    m_log->message(1, "           %12.6f%12.6f%12.6f%12.6f\n",
               maxTerr, avTerr, maxTberr, avTberr);

    if (report_file.is_set()) {
      err.clear_all_strings(); err.clear_all_doubles(); err.set_string("units", "1");
      err.set_name("maximum_temperature");
      err.set_string("units", "Kelvin");
      err.set_string("long_name", "maximum ice temperature error");
      io::write_timeseries(nc, err, (size_t)start, maxTerr);

      err.set_name("average_temperature");
      err.set_string("long_name", "average ice temperature error");
      io::write_timeseries(nc, err, (size_t)start, avTerr);

      err.set_name("maximum_bedrock_temperature");
      err.set_string("long_name", "maximum bedrock temperature error");
      io::write_timeseries(nc, err, (size_t)start, maxTberr);

      err.set_name("average_bedrock_temperature");
      err.set_string("long_name", "average bedrock temperature error");
      io::write_timeseries(nc, err, (size_t)start, avTberr);
    }
  }

  // strain_heating errors if appropriate; reported in 10^6 J/(s m^3)
  if ((testname == 'F') || (testname == 'G')) {
    double max_strain_heating_error, av_strain_heating_error;
    compute_strain_heating_errors(max_strain_heating_error, av_strain_heating_error);
    m_log->message(1,
               "Sigma     :      maxSig       avSig\n");
    m_log->message(1, "           %12.6f%12.6f\n",
               max_strain_heating_error*1.0e6, av_strain_heating_error*1.0e6);

    if (report_file.is_set()) {
      err.clear_all_strings(); err.clear_all_doubles(); err.set_string("units", "1");
      err.set_name("maximum_sigma");
      err.set_string("units", "J s-1 m-3");
      err.set_string("glaciological_units", "1e6 J s-1 m-3");
      err.set_string("long_name", "maximum strain heating error");
      io::write_timeseries(nc, err, (size_t)start, max_strain_heating_error);

      err.set_name("average_sigma");
      err.set_string("long_name", "average strain heating error");
      io::write_timeseries(nc, err, (size_t)start, av_strain_heating_error);
    }
  }

  // surface velocity errors if exact values are available; reported in m year-1
  if ((testname == 'F') || (testname == 'G')) {
    double maxUerr, avUerr, maxWerr, avWerr;
    computeSurfaceVelocityErrors(maxUerr, avUerr, maxWerr, avWerr);
    m_log->message(1,
               "surf vels :     maxUvec      avUvec        maxW         avW\n");
    m_log->message(1,
               "           %12.6f%12.6f%12.6f%12.6f\n",
               convert(m_sys, maxUerr, "m second-1", "m year-1"),
               convert(m_sys, avUerr, "m second-1", "m year-1"),
               convert(m_sys, maxWerr, "m second-1", "m year-1"),
               convert(m_sys, avWerr, "m second-1", "m year-1"));

    if (report_file.is_set()) {
      err.clear_all_strings(); err.clear_all_doubles(); err.set_string("units", "1");
      err.set_name("maximum_surface_velocity");
      err.set_string("long_name", "maximum ice surface horizontal velocity error");
      err.set_string("units", "m second-1");
      err.set_string("glaciological_units", "meters year-1");
      io::write_timeseries(nc, err, (size_t)start, maxUerr);

      err.set_name("average_surface_velocity");
      err.set_string("long_name", "average ice surface horizontal velocity error");
      io::write_timeseries(nc, err, (size_t)start, avUerr);

      err.set_name("maximum_surface_w");
      err.set_string("long_name", "maximum ice surface vertical velocity error");
      io::write_timeseries(nc, err, (size_t)start, maxWerr);

      err.set_name("average_surface_w");
      err.set_string("long_name", "average ice surface vertical velocity error");
      io::write_timeseries(nc, err, (size_t)start, avWerr);
    }
  }

  // basal melt rate errors if appropriate; reported in m year-1
  if (testname == 'O') {
    double maxbmelterr, minbmelterr;
    computeBasalMeltRateErrors(maxbmelterr, minbmelterr);
    if (maxbmelterr != minbmelterr) {
      m_log->message(1,
                 "IceCompModel WARNING: unexpected Test O situation: max and min of bmelt error\n"
                 "  are different: maxbmelterr = %f, minbmelterr = %f\n",
                 convert(m_sys, maxbmelterr, "m second-1", "m year-1"),
                 convert(m_sys, minbmelterr, "m second-1", "m year-1"));
    }
    m_log->message(1,
               "basal melt:  max\n");
    m_log->message(1, "           %11.5f\n",
               convert(m_sys, maxbmelterr, "m second-1", "m year-1"));

    if (report_file.is_set()) {
      err.clear_all_strings(); err.clear_all_doubles(); err.set_string("units", "1");
      err.set_name("maximum_basal_melt_rate");
      err.set_string("units", "m second-1");
      err.set_string("glaciological_units", "meters year-1");
      io::write_timeseries(nc, err, (size_t)start, maxbmelterr);
    }
  }

  if (report_file.is_set()) {
    nc.close();
  }

  m_log->message(1, "NUM ERRORS DONE\n");
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
void IceCompModel::test_V_init() {

  {
    // initialize the bed topography
    IceModelVec2S bed_topography;
    bed_topography.create(m_grid, "topg", WITHOUT_GHOSTS);
    bed_topography.set(-1000);
    m_beddef->set_elevation(bed_topography);
  }

  // set SSA boundary conditions:
  double upstream_velocity = convert(m_sys, 300.0, "m year-1", "m second-1"),
    upstream_thk = 600.0;

  IceModelVec::AccessList list;
  list.add(m_ice_thickness);
  list.add(m_ssa_dirichlet_bc_mask);
  list.add(m_ssa_dirichlet_bc_values);

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    if (i <= 2) {
      m_ssa_dirichlet_bc_mask(i,j) = 1;
      m_ssa_dirichlet_bc_values(i,j)  = Vector2(upstream_velocity, 0.0);
      m_ice_thickness(i, j) = upstream_thk;
    } else {
      m_ssa_dirichlet_bc_mask(i,j) = 0;
      m_ssa_dirichlet_bc_values(i,j)  = Vector2(0.0, 0.0);
      m_ice_thickness(i, j) = 0;
    }
  }

  m_ssa_dirichlet_bc_mask.update_ghosts();

  m_ssa_dirichlet_bc_values.update_ghosts();

  m_ice_thickness.update_ghosts();
}

} // end of namespace pism
