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

#include "EISMINTII.hh"
#include "pism/coupler/AtmosphereModel.hh"
#include "pism/util/ConfigInterface.hh"
#include "pism/util/pism_options.hh"
#include "pism/util/error_handling.hh"
#include "pism/util/IceGrid.hh"
#include "pism/util/MaxTimestep.hh"

namespace pism {
namespace surface {

EISMINTII::EISMINTII(IceGrid::ConstPtr g, int experiment)
  : PSFormulas(g), m_experiment(experiment) {
  // empty
}

EISMINTII::~EISMINTII() {
  // empty
}

void EISMINTII::init_impl(const Geometry &geometry) {
  (void) geometry;

  using units::convert;

  m_log->message(2,
             "setting parameters for surface mass balance"
             " and temperature in EISMINT II experiment %c ... \n",
             m_experiment);

  // EISMINT II specified values for parameters
  m_S_b = convert(m_sys, 1.0e-2 * 1e-3, "year-1", "second-1"); // Grad of accum rate change
  m_S_T = 1.67e-2 * 1e-3;       // K/m  Temp gradient

  // these are for A,E,G,H,I,K:
  m_M_max = convert(m_sys, 0.5, "m year-1", "m second-1"); // Max accumulation
  m_R_el  = 450.0e3;            // Distance to equil line (SMB=0)
  m_T_min = 238.15;

  switch (m_experiment) {
  case 'B':                     // supposed to start from end of experiment A and:
    m_T_min = 243.15;
    break;
  case 'C':
  case 'J':
  case 'L':                     // supposed to start from end of experiment A (for C;
    //   resp I and K for J and L) and:
    m_M_max = convert(m_sys, 0.25, "m year-1", "m second-1");
    m_R_el  = 425.0e3;
    break;
  case 'D':                     // supposed to start from end of experiment A and:
    m_R_el  = 425.0e3;
    break;
  case 'F':                     // start with zero ice and:
    m_T_min = 223.15;
    break;
  }

  // if user specifies Tmin, Tmax, Mmax, Sb, ST, Rel, then use that (override above)
  m_T_min = options::Real("-Tmin", "T min, Kelvin", m_T_min);

  options::Real Mmax("-Mmax", "Maximum accumulation, m year-1",
                     convert(m_sys, m_M_max, "m second-1", "m year-1"));
  if (Mmax.is_set()) {
    m_M_max = convert(m_sys, Mmax, "m year-1", "m second-1");
  }

  options::Real Sb("-Sb", "radial gradient of accumulation rate, (m year-1)/km",
                   convert(m_sys, m_S_b,   "m second-1/m", "m year-1/km"));
  if (Sb.is_set()) {
    m_S_b = convert(m_sys, Sb, "m year-1/km", "m second-1/m");
  }

  options::Real ST("-ST", "radial gradient of surface temperature, K/km",
                   convert(m_sys, m_S_T, "K/m", "K/km"));
  if (ST.is_set()) {
    m_S_T = convert(m_sys, ST, "K/km", "K/m");
  }

  options::Real Rel("-Rel", "radial distance to equilibrium line, km",
                    convert(m_sys, m_R_el, "m", "km"));
  if (Rel.is_set()) {
    m_R_el = convert(m_sys, Rel, "km", "m");
  }

  initialize_using_formulas();
}

MaxTimestep EISMINTII::max_timestep_impl(double t) const {
  (void) t;
  return MaxTimestep("surface EISMINT-II");
}

void EISMINTII::initialize_using_formulas() {

  // center of the accumulation and surface temperature patterns
  double cx = 0.0, cy = 0.0;
  if (m_experiment == 'E') {
    cx += 100.0e3;
    cy += 100.0e3;
  }

  IceModelVec::AccessList list{m_temperature.get(), m_mass_flux.get()};

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    const double r = sqrt(PetscSqr(m_grid->x(i) - cx) + PetscSqr(m_grid->y(j) - cy));

    // accumulation (formula (7) in [Payne et al 2000])
    (*m_mass_flux)(i,j) = std::min(m_M_max, m_S_b * (m_R_el-r));

    // surface temperature (formula (8) in [Payne et al 2000])
    (*m_temperature)(i,j) = m_T_min + m_S_T * r;
  }

  // convert from "m second-1" to "kg m-2 s-1"
  m_mass_flux->scale(m_config->get_double("constants.ice.density"));
}

void EISMINTII::update_impl(const Geometry &geometry, double t, double dt) {
  (void) t;
  (void) dt;
  (void) geometry;

  // do nothing (but an implementation is required)
}

} // end of namespace surface
} // end of namespace pism
