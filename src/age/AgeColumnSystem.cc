/* Copyright (C) 2016, 2017 PISM Authors
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

#include "AgeColumnSystem.hh"

#include "pism/util/error_handling.hh"

namespace pism {

AgeColumnSystem::AgeColumnSystem(const std::vector<double>& storage_grid,
                                 const std::string &my_prefix,
                                 double dx, double dy, double dt,
                                 const array::Array3D &age,
                                 const array::Array3D &u3,
                                 const array::Array3D &v3,
                                 const array::Array3D &w3)
  : columnSystemCtx(storage_grid, my_prefix, dx, dy, dt, u3, v3, w3),
    m_age3(age) {

  size_t Mz = m_z.size();
  m_A.resize(Mz);
  m_A_n.resize(Mz);
  m_A_e.resize(Mz);
  m_A_s.resize(Mz);
  m_A_w.resize(Mz);

  m_nu = m_dt / m_dz; // derived constant
}

void AgeColumnSystem::init(int i, int j, double thickness) {
  init_column(i, j, thickness);

  if (m_ks == 0) {
    return;
  }

  coarse_to_fine(m_u3, i, j, &m_u[0]);
  coarse_to_fine(m_v3, i, j, &m_v[0]);
  coarse_to_fine(m_w3, i, j, &m_w[0]);

  coarse_to_fine(m_age3, m_i, m_j,   &m_A[0]);
  coarse_to_fine(m_age3, m_i, m_j+1, &m_A_n[0]);
  coarse_to_fine(m_age3, m_i+1, m_j, &m_A_e[0]);
  coarse_to_fine(m_age3, m_i, m_j-1, &m_A_s[0]);
  coarse_to_fine(m_age3, m_i-1, m_j, &m_A_w[0]);
}

//! First-order upwind scheme with implicit in the vertical: one column solve.
/*!
  The PDE being solved is
  \f[ \frac{\partial \tau}{\partial t} + \frac{\partial}{\partial x}\left(u \tau\right) + \frac{\partial}{\partial y}\left(v \tau\right) + \frac{\partial}{\partial z}\left(w \tau\right) = 1. \f]
 */
void AgeColumnSystem::solve(std::vector<double> &x) {

  TridiagonalSystem &S = *m_solver;

  // set up system: 0 <= k < m_ks
  for (unsigned int k = 0; k < m_ks; k++) {
    // do lowest-order upwinding, explicitly for horizontal
    S.RHS(k) =  (m_u[k] < 0 ?
                 m_u[k] * (m_A_e[k] -  m_A[k]) / m_dx :
                 m_u[k] * (m_A[k]  - m_A_w[k]) / m_dx);
    S.RHS(k) += (m_v[k] < 0 ?
                 m_v[k] * (m_A_n[k] -  m_A[k]) / m_dy :
                 m_v[k] * (m_A[k]  - m_A_s[k]) / m_dy);
    // note it is the age eqn: dage/dt = 1.0 and we have moved the hor.
    //   advection terms over to right:
    S.RHS(k) = m_A[k] + m_dt * (1.0 - S.RHS(k));

    // do lowest-order upwinding, *implicitly* for vertical
    double AA = m_nu * m_w[k];
    if (k > 0) {
      if (AA >= 0) { // upward velocity
        S.L(k) = - AA;
        S.D(k) = 1.0 + AA;
        S.U(k) = 0.0;
      } else { // downward velocity; note  -AA >= 0
        S.L(k) = 0.0;
        S.D(k) = 1.0 - AA;
        S.U(k) = + AA;
      }
    } else { // k == 0 case
      // note L[0] is not used
      if (AA > 0) { // if strictly upward velocity apply boundary condition:
                    // age = 0 because ice is being added to base
        S.D(0) = 1.0;
        S.U(0) = 0.0;
        S.RHS(0) = 0.0;
      } else { // downward velocity; note  -AA >= 0
        S.D(0) = 1.0 - AA;
        S.U(0) = + AA;
        // keep rhs[0] as is
      }
    }
  }  // done "set up system: 0 <= k < m_ks"

  // surface b.c. at m_ks
  if (m_ks > 0) {
    S.L(m_ks) = 0;
    S.D(m_ks) = 1.0;   // ignore U[m_ks]
    S.RHS(m_ks) = 0.0;  // age zero at surface
  }

  // solve it
  try {
    S.solve(m_ks + 1, x);
  }
  catch (RuntimeError &e) {
    e.add_context("solving the tri-diagonal system (AgeColumnSystem) at (%d, %d)\n"
                  "saving system to m-file... ", m_i, m_j);
    reportColumnZeroPivotErrorMFile(m_ks + 1);
    throw;
  }

  // x[k] contains age for k=0,...,ks, but set age of ice above (and
  // at) surface to zero years
  for (unsigned int k = m_ks + 1; k < x.size(); k++) {
    x[k] = 0.0;
  }
}

} // end of namespace pism
