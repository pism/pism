/* Copyright (C) 2014, 2015 PISM Authors
 *
 * This file is part of PISM.
 *
 * PISM is free software; you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation; either version 3 of the License, or (at your option) any later
 * version.
 *
 * PISM is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License
 * along with PISM; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 */

#include <algorithm>
#include <vector>

#include "PSVerification.hh"
#include "PISMAtmosphere.hh"
#include "flowlaws.hh"
#include "enthalpyConverter.hh"
#include "PISMTime.hh"
#include "PISMConfig.hh"

#include "tests/exactTestsABCDE.h"
#include "tests/exactTestsFG.h"
#include "tests/exactTestH.h"
#include "tests/exactTestL.h"

#include "error_handling.hh"

namespace pism {

const double PSVerification::ablationRateOutside = 0.02; // m/year
const double PSVerification::secpera = 3.15569259747e7;

const double PSVerification::ST = 1.67e-5;
const double PSVerification::Tmin = 223.15;  // K
const double PSVerification::LforFG = 750000; // m
const double PSVerification::ApforG = 200; // m

PSVerification::PSVerification(const IceGrid &g,
                               EnthalpyConverter *EC, int test)
  : PSFormulas(g), m_testname(test), m_EC(EC) {
  // empty
}

PSVerification::~PSVerification() {
  // empty
}

void PSVerification::init() {
  // Make sure that ice surface temperature and climatic mass balance
  // get initialized at the beginning of the run (as far as I can tell
  // this affects zero-length runs only).
  update(m_grid.time->current(), 0);
}

/** Initialize climate inputs of tests K and O.
 * 
 * @return 0 on success
 */
void PSVerification::update_KO() {

  m_climatic_mass_balance.set(0.0);
  m_ice_surface_temp.set(223.15);
}

/** Update the test L climate input (once).
 *
 * Unlike other `update_...()` methods, this one uses [kg m-2 s-1]
 * as units of the climatic_mass_balance.
 *
 * @return 0 on success
 */
void PSVerification::update_L() {
  double     A0, T0;

  ThermoGlenArrIce tgaIce(m_grid.com, "sia_", m_config, m_EC);

  // compute T so that A0 = A(T) = Acold exp(-Qcold/(R T))  (i.e. for ThermoGlenArrIce);
  // set all temps to this constant
  A0 = 1.0e-16/secpera;    // = 3.17e-24  1/(Pa^3 s);  (EISMINT value) flow law parameter
  T0 = tgaIce.tempFromSoftness(A0);

  m_ice_surface_temp.set(T0);

  const double
    ice_density = m_config.get("ice_density"),
    a0          = m_grid.convert(0.3, "m/year", "m/s"),
    L           = 750e3,
    Lsqr        = L * L;

  IceModelVec::AccessList list(m_climatic_mass_balance);
  for (Points p(m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    double r = radius(m_grid, i, j);
    m_climatic_mass_balance(i, j) = a0 * (1.0 - (2.0 * r * r / Lsqr));

    m_climatic_mass_balance(i, j) *= ice_density; // convert to [kg m-2 s-1]
  }
}

void PSVerification::update_V() {

  // initialize temperature; the value used does not matter
  m_ice_surface_temp.set(273.15);

  // initialize mass balance:
  m_climatic_mass_balance.set(0.0);
}

void PSVerification::update(PetscReal t, PetscReal dt) {

  (void) dt;

  switch (m_testname) {
  case 'A':
  case 'B':
  case 'C':
  case 'D':
  case 'E':
  case 'H':
    update_ABCDEH(t);
    break;
  case 'F':
  case 'G':
    update_FG(t);
    break;
  case 'K':
  case 'O':
    update_KO();
    break;
  case 'L':
    {
      update_L();
      // return here; note update_L() uses correct units
      return;
    }
  case 'V':
    update_V();
    break;
  default:
    throw RuntimeError::formatted("Test %c is not implemented.", m_testname);
  }

  // convert from [m/s] to [kg m-2 s-1]
  m_climatic_mass_balance.scale(m_config.get("ice_density"));
}

/** Update climate inputs for tests A, B, C, D, E, H.
 *
 * @return 0 on success
 */
void PSVerification::update_ABCDEH(double time) {
  double         A0, T0, H, accum, dummy1, dummy2, dummy3;

  double f = m_config.get("ice_density") / m_config.get("lithosphere_density");

  ThermoGlenArrIce tgaIce(m_grid.com, "sia_", m_config, m_EC);

  // compute T so that A0 = A(T) = Acold exp(-Qcold/(R T))  (i.e. for ThermoGlenArrIce);
  // set all temps to this constant
  A0 = 1.0e-16/secpera;    // = 3.17e-24  1/(Pa^3 s);  (EISMINT value) flow law parameter
  T0 = tgaIce.tempFromSoftness(A0);

  m_ice_surface_temp.set(T0);

  IceModelVec::AccessList list(m_climatic_mass_balance);
  for (Points p(m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    double xx = m_grid.x(i), yy = m_grid.y(j),
      r = radius(m_grid, i, j);
    switch (m_testname) {
    case 'A':
      exactA(r, &H, &accum);
      break;
    case 'B':
      exactB(time, r, &H, &accum);
      break;
    case 'C':
      exactC(time, r, &H, &accum);
      break;
    case 'D':
      exactD(time, r, &H, &accum);
      break;
    case 'E':
      exactE(xx, yy, &H, &accum, &dummy1, &dummy2, &dummy3);
      break;
    case 'H':
      exactH(f, time, r, &H, &accum);
      break;
    default:
      throw RuntimeError::formatted("test must be A, B, C, D, E, or H, got %c",
                                    m_testname);
    }
    m_climatic_mass_balance(i, j) = accum;
  }
}

void PSVerification::update_FG(double time) {
  unsigned int   Mz = m_grid.Mz();
  double         H, accum;

  std::vector<double> dummy1(Mz), dummy2(Mz), dummy3(Mz), dummy4(Mz), dummy5(Mz);

  IceModelVec::AccessList list(m_climatic_mass_balance);
  list.add(m_ice_surface_temp);

  for (Points p(m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    double r = std::max(radius(m_grid, i, j), 1.0); // avoid singularity at origin

    m_ice_surface_temp(i, j) = Tmin + ST * r;

    if (r > LforFG - 1.0) { // if (essentially) outside of sheet
      accum = - ablationRateOutside / secpera;
    } else {
      if (m_testname == 'F') {
        bothexact(0.0, r, &(m_grid.z()[0]), Mz, 0.0,
                  &H, &accum, &dummy5[0], &dummy1[0], &dummy2[0], &dummy3[0], &dummy4[0]);
      } else {
        bothexact(time, r, &(m_grid.z()[0]), Mz, ApforG,
                  &H, &accum, &dummy5[0], &dummy1[0], &dummy2[0], &dummy3[0], &dummy4[0]);
      }
    }
    m_climatic_mass_balance(i, j) = accum;

  }
}

} // end of namespace pism
