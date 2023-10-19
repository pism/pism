/* Copyright (C) 2016, 2017, 2023 PISM Authors
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

#include "pism/verification/TemperatureModel_Verification.hh"

#include "pism/util/error_handling.hh"
#include "pism/verification/tests/exactTestsFG.hh"
#include "pism/verification/tests/exactTestK.h"
#include "pism/verification/tests/exactTestO.h"
#include "pism/energy/utilities.hh"
#include "pism/util/Time.hh"
#include "pism/util/Context.hh"

namespace pism {
namespace energy {

static const double ApforG = 200; // m
static const double LforFG = 750000; // m
static const double ST     = 1.67e-5;
static const double Tmin   = 223.15; // K

TemperatureModel_Verification::TemperatureModel_Verification(
    std::shared_ptr<const Grid> grid,
    std::shared_ptr<const stressbalance::StressBalance> stress_balance, int testname,
    bool bedrock_is_ice)
    : TemperatureModel(grid, stress_balance),
      m_testname(testname),
      m_bedrock_is_ice(bedrock_is_ice) {
  // empty
}

void TemperatureModel_Verification::initialize_impl(const array::Scalar &basal_melt_rate,
                                                    const array::Scalar &ice_thickness,
                                                    const array::Scalar &surface_temperature,
                                                    const array::Scalar &climatic_mass_balance,
                                                    const array::Scalar &basal_heat_flux) {

  // ignore provided basal melt rate
  (void) basal_melt_rate;

  m_basal_melt_rate.set(0.0);

  switch (m_testname) {
  case 'F':
  case 'G':
    initTestFG();
    break;
  case 'K':
  case 'O':
    initTestsKO();
    break;
  default:
    TemperatureModel::initialize_impl(m_basal_melt_rate, ice_thickness, surface_temperature,
                                      climatic_mass_balance, basal_heat_flux);
  }

  m_ice_temperature.update_ghosts();
  m_basal_melt_rate.update_ghosts();

  // this will update ghosts of m_ice_enthalpy
  compute_enthalpy_cold(m_ice_temperature, ice_thickness, m_ice_enthalpy);
}

void TemperatureModel_Verification::initTestFG() {

  array::AccessScope list{&m_ice_temperature};

  const double time = m_testname == 'F' ? 0.0 : this->time().current();
  const double A    = m_testname == 'F' ? 0.0 : ApforG;

  for (auto p = m_grid->points(); p; p.next()) {
    const int i = p.i(), j = p.j();

    // avoid singularity at origin
    const double r = std::max(grid::radius(*m_grid, i, j), 1.0);

    if (r > LforFG - 1.0) { // if (essentially) outside of sheet
      m_ice_temperature.set_column(i, j, Tmin + ST * r);
    } else {
      TestFGParameters P = exactFG(time, r, m_grid->z(), A);
      m_ice_temperature.set_column(i, j, P.T.data());
    }
  }
}

void TemperatureModel_Verification::initTestsKO() {

  const unsigned int Mz = m_grid->Mz();

  std::vector<double> temperature(Mz);

  const double time = this->time().current();

  // evaluate exact solution in a column; all columns are the same
  for (unsigned int k = 0; k < Mz; k++) {
    if (m_testname == 'K') {
      temperature[k] = exactK(time, m_grid->z(k), m_bedrock_is_ice ? 1 : 0).T;
    } else {
      temperature[k] = exactO(m_grid->z(k)).TT;
    }
  }

  // fill m_ice_temperature
  array::AccessScope list(m_ice_temperature);

  ParallelSection loop(m_grid->com);
  try {
    for (auto p = m_grid->points(); p; p.next()) {
      m_ice_temperature.set_column(p.i(), p.j(), temperature.data());
    }
  } catch (...) {
    loop.failed();
  }
  loop.check();
}


} // end of namespace energy
} // end of namespace pism
