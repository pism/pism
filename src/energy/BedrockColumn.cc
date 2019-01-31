/* Copyright (C) 2019 PISM Authors
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

#include <cassert>

#include "BedrockColumn.hh"

#include "pism/util/ConfigInterface.hh"

namespace pism {
namespace energy {

BedrockColumn::BedrockColumn(const std::string& prefix,
                             const Config& config, double dz, unsigned int M)
  : m_dz(dz), m_M(M), m_system(M, prefix) {

  assert(M > 1);

  const double
    rho = config.get_double("energy.bedrock_thermal_density"),
    c   = config.get_double("energy.bedrock_thermal_specific_heat_capacity");

  m_k   = config.get_double("energy.bedrock_thermal_conductivity");
  m_D   = m_k / (rho * c);
}

BedrockColumn::~BedrockColumn() {
  // empty
}

/*!
 * Advance the heat equation in time.
 *
 * @param[in] dt time step length
 * @param[in] Q_bottom heat flux into the column through the bottom surface
 * @param[in] T_top temperature at the top surface
 * @param[in] T_old current temperature in the column
 * @param[out] T_new output
 *
 * Note: T_old and T_new may point to the same location.
 */
void BedrockColumn::solve(double dt, double Q_bottom, double T_top,
                          const double *T_old, double *T_new) {

  double R = m_D * dt / (m_dz * m_dz);
  double G = -Q_bottom / m_k;

  m_system.L(0)   = 0.0;                 // not used
  m_system.D(0)   = 1.0 + 2.0 * R;
  m_system.U(0)   = -2.0 * R;
  m_system.RHS(0) = T_old[0] - 2.0 * G * m_dz * R;

  unsigned int N = m_M - 1;

  for (unsigned int k = 1; k < N; ++k) {
    m_system.L(k)   = -R;
    m_system.D(k)   = 1.0 + 2.0 * R;
    m_system.U(k)   = -R;
    m_system.RHS(k) = T_old[k];
  }

  m_system.L(N)   = 0.0;
  m_system.D(N)   = 1.0;
  m_system.U(N)   = 0.0;                 // not used
  m_system.RHS(N) = T_top;

  m_system.solve(m_M, T_new);
}

/*!
 * This version of `solve()` is easier to use in Python.
 */
void BedrockColumn::solve(double dt, double Q_bottom, double T_top,
                          const std::vector<double> &T_old,
                          std::vector<double> &result) {
  result.resize(m_M);
  solve(dt, Q_bottom, T_top, T_old.data(), result.data());
}


} // end of namespace energy
} // end of namespace pism
