/* Copyright (C) 2014 PISM Authors
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

#include "PS_EISMINTII.hh"
#include "PISMAtmosphere.hh"
#include "PISMConfig.hh"
#include "pism_options.hh"
#include <stdexcept>
#include "error_handling.hh"

namespace pism {

PS_EISMINTII::PS_EISMINTII(IceGrid &g, int experiment)
  : PSFormulas(g), m_experiment(experiment) {
  // empty
}

PS_EISMINTII::~PS_EISMINTII() {
  // empty
}

void PS_EISMINTII::init() {

  verbPrintf(2, m_grid.com, 
             "setting parameters for surface mass balance"
             " and temperature in EISMINT II experiment %c ... \n", 
             m_experiment);

  // EISMINT II specified values for parameters
  m_S_b = m_grid.convert(1.0e-2 * 1e-3, "1/year", "1/s"); // Grad of accum rate change
  m_S_T = 1.67e-2 * 1e-3;       // K/m  Temp gradient

  // these are for A,E,G,H,I,K:
  m_M_max = m_grid.convert(0.5, "m/year", "m/s"); // Max accumulation
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
    m_M_max = m_grid.convert(0.25, "m/year", "m/s");
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
  bool option_set = false;
  OptionsReal("-Tmin", "T min, Kelvin", m_T_min, option_set);

  double
    myMmax = m_grid.convert(m_M_max, "m/s",        "m/year"),
    mySb   = m_grid.convert(m_S_b,   "m/second/m", "m/year/km"),
    myST   = m_grid.convert(m_S_T,   "K/m",        "K/km"),
    myRel  = m_grid.convert(m_R_el,  "m",          "km");

  OptionsReal("-Mmax", "Maximum accumulation, m/year",
              myMmax, option_set);
  if (option_set) {
    m_M_max = m_grid.convert(myMmax, "m/year", "m/second");
  }

  OptionsReal("-Sb", "radial gradient of accumulation rate, (m/year)/km",
              mySb, option_set);
  if (option_set) {
    m_S_b = m_grid.convert(mySb, "m/year/km", "m/second/m");
  }

  OptionsReal("-ST", "radial gradient of surface temperature, K/km",
              myST, option_set);
  if (option_set) {
    m_S_T = m_grid.convert(myST, "K/km", "K/m");
  }

  OptionsReal("-Rel", "radial distance to equilibrium line, km",
              myRel, option_set);
  if (option_set) {
    m_R_el = m_grid.convert(myRel, "km", "m");
  }

  initialize_using_formulas();
}

void PS_EISMINTII::initialize_using_formulas() {

  PetscScalar cx = m_grid.Lx(), cy = m_grid.Ly();
  if (m_experiment == 'E') {
    // shift center
    cx += 100.0e3;
    cy += 100.0e3;
  }

  // now fill in accum and surface temp
  IceModelVec::AccessList list;
  list.add(m_ice_surface_temp);
  list.add(m_climatic_mass_balance);

  for (Points p(m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    // r is distance from center of grid; if E then center is shifted (above)
    const double r = sqrt(PetscSqr(-cx + m_grid.dx()*i)
                          + PetscSqr(-cy + m_grid.dy()*j));
    // set accumulation from formula (7) in (Payne et al 2000)
    m_climatic_mass_balance(i,j) = std::min(m_M_max, m_S_b * (m_R_el-r));
    // set surface temperature
    m_ice_surface_temp(i,j) = m_T_min + m_S_T * r;  // formula (8) in (Payne et al 2000)
  }

  // convert from [m/s] to [kg m-2 s-1]
  m_climatic_mass_balance.scale(m_config.get("ice_density"));
}

void PS_EISMINTII::update(PetscReal t, PetscReal dt) {
  (void) t;
  (void) dt;

  // do nothing (but an implementation is required)
}

} // end of namespace pism
