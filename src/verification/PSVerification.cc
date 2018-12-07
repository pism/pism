/* Copyright (C) 2014, 2015, 2016, 2017, 2018 PISM Authors
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
#include "pism/coupler/AtmosphereModel.hh"
#include "pism/rheology/PatersonBuddCold.hh"
#include "pism/util/EnthalpyConverter.hh"
#include "pism/util/Time.hh"
#include "pism/util/ConfigInterface.hh"

#include "tests/exactTestsABCD.h"
#include "tests/exactTestsFG.hh"
#include "tests/exactTestH.h"
#include "tests/exactTestL.hh"

#include "pism/util/error_handling.hh"
#include "pism/util/IceGrid.hh"
#include "pism/util/MaxTimestep.hh"

namespace pism {
namespace surface {

const double Verification::ablationRateOutside = 0.02; // m year-1
const double Verification::secpera = 3.15569259747e7;

const double Verification::ST = 1.67e-5;
const double Verification::Tmin = 223.15;  // K
const double Verification::LforFG = 750000; // m
const double Verification::ApforG = 200; // m

Verification::Verification(IceGrid::ConstPtr g,
                           EnthalpyConverter::Ptr EC, int test)
  : PSFormulas(g), m_testname(test), m_EC(EC) {
  // empty
}

Verification::~Verification() {
  // empty
}

void Verification::init_impl(const Geometry &geometry) {
  // Make sure that ice surface temperature and climatic mass balance
  // get initialized at the beginning of the run (as far as I can tell
  // this affects zero-length runs only).
  update(geometry, m_grid->ctx()->time()->current(), 0);
}

void Verification::define_model_state_impl(const PIO &output) const {
  m_mass_flux->define(output);
  m_temperature->define(output);
}

void Verification::write_model_state_impl(const PIO &output) const {
  m_mass_flux->write(output);
  m_temperature->write(output);
}

MaxTimestep Verification::max_timestep_impl(double t) const {
  (void) t;
  return MaxTimestep("verification surface model");
}

/** Initialize climate inputs of tests K and O.
 * 
 * @return 0 on success
 */
void Verification::update_KO() {

  m_mass_flux->set(0.0);
  m_temperature->set(223.15);
}

/** Update the test L climate input (once).
 *
 * Unlike other `update_...()` methods, this one uses [kg m-2 s-1]
 * as units of the climatic_mass_balance.
 *
 * @return 0 on success
 */
void Verification::update_L() {
  double     A0, T0;

  rheology::PatersonBuddCold tgaIce("stress_balance.sia.", *m_config, m_EC);

  // compute T so that A0 = A(T) = Acold exp(-Qcold/(R T))  (i.e. for PatersonBuddCold);
  // set all temps to this constant
  A0 = 1.0e-16/secpera;    // = 3.17e-24  1/(Pa^3 s);  (EISMINT value) flow law parameter
  T0 = tgaIce.tempFromSoftness(A0);

  m_temperature->set(T0);

  const double
    ice_density = m_config->get_double("constants.ice.density"),
    a0          = units::convert(m_sys, 0.3, "m year-1", "m second-1"),
    L           = 750e3,
    Lsqr        = L * L;

  IceModelVec::AccessList list(*m_mass_flux);
  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    double r = radius(*m_grid, i, j);
    (*m_mass_flux)(i, j) = a0 * (1.0 - (2.0 * r * r / Lsqr));

    (*m_mass_flux)(i, j) *= ice_density; // convert to [kg m-2 s-1]
  }
}

void Verification::update_V() {

  // initialize temperature; the value used does not matter
  m_temperature->set(273.15);

  // initialize mass balance:
  m_mass_flux->set(0.0);
}

void Verification::update_impl(const Geometry &geometry, double t, double dt) {
  (void) geometry;
  (void) dt;

  switch (m_testname) {
  case 'A':
  case 'B':
  case 'C':
  case 'D':
  case 'H':
    update_ABCDH(t);
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
    throw RuntimeError::formatted(PISM_ERROR_LOCATION, "Test %c is not implemented.", m_testname);
  }

  // convert from [m second-1] to [kg m-2 s-1]
  m_mass_flux->scale(m_config->get_double("constants.ice.density"));
}

/** Update climate inputs for tests A, B, C, D, E, H.
 *
 * @return 0 on success
 */
void Verification::update_ABCDH(double time) {
  double A0, T0, accum;

  double f = m_config->get_double("constants.ice.density") / m_config->get_double("bed_deformation.mantle_density");

  rheology::PatersonBuddCold tgaIce("stress_balance.sia.", *m_config, m_EC);

  // compute T so that A0 = A(T) = Acold exp(-Qcold/(R T))  (i.e. for PatersonBuddCold);
  // set all temps to this constant
  A0 = 1.0e-16/secpera;    // = 3.17e-24  1/(Pa^3 s);  (EISMINT value) flow law parameter
  T0 = tgaIce.tempFromSoftness(A0);

  m_temperature->set(T0);

  IceModelVec::AccessList list(*m_mass_flux);
  ParallelSection loop(m_grid->com);
  try {
    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      const double r = radius(*m_grid, i, j);
      switch (m_testname) {
      case 'A':
        accum = exactA(r).M;
        break;
      case 'B':
        accum = exactB(time, r).M;
        break;
      case 'C':
        accum = exactC(time, r).M;
        break;
      case 'D':
        accum = exactD(time, r).M;
        break;
      case 'H':
        accum = exactH(f, time, r).M;
        break;
      default:
        throw RuntimeError::formatted(PISM_ERROR_LOCATION, "test must be A, B, C, D, or H, got %c",
                                      m_testname);
      }
      (*m_mass_flux)(i, j) = accum;
    }
  } catch (...) {
    loop.failed();
  }
  loop.check();
}

void Verification::update_FG(double time) {

  const double t = m_testname == 'F' ? 0.0 : time;
  const double A = m_testname == 'F' ? 0.0 : ApforG;

  IceModelVec::AccessList list{m_mass_flux.get(), m_temperature.get()};

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    // avoid singularity at origin
    const double r = std::max(radius(*m_grid, i, j), 1.0);

    (*m_temperature)(i, j) = Tmin + ST * r;

    if (r > LforFG - 1.0) {
      // if (essentially) outside of sheet
      (*m_mass_flux)(i, j) = - ablationRateOutside / secpera;
    } else {
      (*m_mass_flux)(i, j) = exactFG(t, r, m_grid->z(), A).M;
    }
  }
}

} // end of namespace surface
} // end of namespace pism
