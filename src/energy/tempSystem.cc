// Copyright (C) 2004-2011, 2013, 2014, 2015, 2016, 2017, 2018 Jed Brown, Ed Bueler and Constantine Khroulev
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

#include <cassert>

#include "pism/util/pism_utilities.hh"
#include "pism/util/iceModelVec.hh"
#include "tempSystem.hh"
#include "pism/util/Mask.hh"

#include "pism/util/error_handling.hh"

namespace pism {
namespace energy {

tempSystemCtx::tempSystemCtx(const std::vector<double>& storage_grid,
                             const std::string &prefix,
                             double dx, double dy, double dt,
                             const Config &config,
                             const IceModelVec3 &T3,
                             const IceModelVec3 &u3,
                             const IceModelVec3 &v3,
                             const IceModelVec3 &w3,
                             const IceModelVec3 &strain_heating3)
  : columnSystemCtx(storage_grid, prefix, dx, dy, dt, u3, v3, w3),
    m_T3(T3),
    m_strain_heating3(strain_heating3) {

  // set flags to indicate nothing yet set
  m_surfBCsValid      = false;
  m_basalBCsValid     = false;

  size_t Mz = m_z.size();
  m_T.resize(Mz);
  m_strain_heating.resize(Mz);

  m_T_n.resize(Mz);
  m_T_e.resize(Mz);
  m_T_s.resize(Mz);
  m_T_w.resize(Mz);

  // set physical constants
  m_ice_density = config.get_double("constants.ice.density");
  m_ice_c       = config.get_double("constants.ice.specific_heat_capacity");
  m_ice_k       = config.get_double("constants.ice.thermal_conductivity");

  // set derived constants
  m_nu    = m_dt / m_dz;
  m_rho_c_I = m_ice_density * m_ice_c;
  m_iceK    = m_ice_k / m_rho_c_I;
  m_iceR    = m_iceK * m_dt / (m_dz*m_dz);
}

void tempSystemCtx::initThisColumn(int i, int j, bool is_marginal, MaskValue mask,
                                   double ice_thickness) {

  m_is_marginal = is_marginal;
  m_mask = mask;

  init_column(i, j, ice_thickness);

  if (m_ks == 0) {
    return;
  }

  coarse_to_fine(m_u3, m_i, m_j, &m_u[0]);
  coarse_to_fine(m_v3, m_i, m_j, &m_v[0]);
  coarse_to_fine(m_w3, m_i, m_j, &m_w[0]);
  coarse_to_fine(m_strain_heating3, m_i, m_j, &m_strain_heating[0]);
  coarse_to_fine(m_T3, m_i, m_j, &m_T[0]);

  coarse_to_fine(m_T3, m_i, m_j+1, &m_T_n[0]);
  coarse_to_fine(m_T3, m_i+1, m_j, &m_T_e[0]);
  coarse_to_fine(m_T3, m_i, m_j-1, &m_T_s[0]);
  coarse_to_fine(m_T3, m_i-1, m_j, &m_T_w[0]);

  m_lambda = compute_lambda();
}

void tempSystemCtx::setSurfaceBoundaryValuesThisColumn(double my_Ts) {
  // allow setting surface BCs only once:
  assert(not m_surfBCsValid);

  m_Ts           = my_Ts;
  m_surfBCsValid = true;
}


void tempSystemCtx::setBasalBoundaryValuesThisColumn(double my_G0,
                                                     double my_Tshelfbase, double my_Rb) {
  // allow setting basal BCs only once:
  assert(not m_basalBCsValid);

  m_G0            = my_G0;
  m_Tshelfbase    = my_Tshelfbase;
  m_Rb            = my_Rb;
  m_basalBCsValid = true;
}

double tempSystemCtx::compute_lambda() {
  double result = 1.0; // start with centered implicit for more accuracy
  const double epsilon = 1e-6 / 3.15569259747e7;

  for (unsigned int k = 0; k <= m_ks; k++) {
    const double denom = (fabs(m_w[k]) + epsilon) * m_ice_density * m_ice_c * m_dz;
    result = std::min(result, 2.0 * m_ice_k / denom);
  }
  return result;
}

void tempSystemCtx::solveThisColumn(std::vector<double> &x) {

  TridiagonalSystem &S = *m_solver;

  assert(m_surfBCsValid == true);
  assert(m_basalBCsValid == true);

  // bottom of ice; k=0 eqn
  if (m_ks == 0) { // no ice; set m_T[0] to surface temp if grounded
    // note L[0] not allocated
    S.D(0) = 1.0;
    S.U(0) = 0.0;
    // if floating and no ice then worry only about bedrock temps
    if (mask::ocean(m_mask)) {
      // essentially no ice but floating ... ask OceanCoupler
      S.RHS(0) = m_Tshelfbase;
    } else { // top of bedrock sees atmosphere
      S.RHS(0) = m_Ts;
    }
  } else { // m_ks > 0; there is ice
    // for w, always difference *up* from base, but make it implicit
    if (mask::ocean(m_mask)) {
      // just apply Dirichlet condition to base of column of ice in an ice shelf
      // note that L[0] is not used
      S.D(0) = 1.0;
      S.U(0) = 0.0;
      S.RHS(0) = m_Tshelfbase; // set by OceanCoupler
    } else {
      // there is *grounded* ice; from FV across interface
      S.RHS(0) = m_T[0] + m_dt * (m_Rb / (m_rho_c_I * m_dz));

      double Sigma = 0.0;
      if (not m_is_marginal) {
        Sigma = m_strain_heating[0];
      }
      S.RHS(0) += m_dt * 0.5 * Sigma / m_rho_c_I;

      double UpTu = 0.0;
      double UpTv = 0.0;
      if (not m_is_marginal) {
        UpTu = (m_u[0] < 0 ?
                m_u[0] * (m_T_e[0] -  m_T[0]) / m_dx :
                m_u[0] * (m_T[0]  - m_T_w[0]) / m_dx);
        UpTv = (m_v[0] < 0 ?
                m_v[0] * (m_T_n[0] -  m_T[0]) / m_dy :
                m_v[0] * (m_T[0]  - m_T_s[0]) / m_dy);
      }
      S.RHS(0) -= m_dt  * (0.5 * (UpTu + UpTv));

      // vertical upwinding
      // L[0] = 0.0;  (is not used)
      S.D(0) = 1.0 + 2.0 * m_iceR;
      S.U(0) = - 2.0 * m_iceR;
      if (m_w[0] < 0.0) { // velocity downward: add velocity contribution
        const double AA = m_dt * m_w[0] / (2.0 * m_dz);
        S.D(0) -= AA;
        S.U(0) += AA;
      }
      // apply geothermal flux G0 here
      S.RHS(0) += 2.0 * m_dt * m_G0 / (m_rho_c_I * m_dz);
    }
  }

  // generic ice segment; build 1:m_ks-1 eqns
  for (unsigned int k = 1; k < m_ks; k++) {
    const double AA = m_nu * m_w[k];
    if (m_w[k] >= 0.0) {  // velocity upward
      S.L(k) = - m_iceR - AA * (1.0 - m_lambda/2.0);
      S.D(k) = 1.0 + 2.0 * m_iceR + AA * (1.0 - m_lambda);
      S.U(k) = - m_iceR + AA * (m_lambda/2.0);
    } else {  // velocity downward
      S.L(k) = - m_iceR - AA * (m_lambda/2.0);
      S.D(k) = 1.0 + 2.0 * m_iceR - AA * (1.0 - m_lambda);
      S.U(k) = - m_iceR + AA * (1.0 - m_lambda/2.0);
    }
    S.RHS(k) = m_T[k];

    double Sigma = 0.0;
    if (not m_is_marginal) {
      Sigma = m_strain_heating[k];
    }

    double UpTu = 0.0;
    double UpTv = 0.0;
    if (not m_is_marginal) {
      UpTu = (m_u[k] < 0 ?
              m_u[k] * (m_T_e[k] -  m_T[k]) / m_dx :
              m_u[k] * (m_T[k]  - m_T_w[k]) / m_dx);
      UpTv = (m_v[k] < 0 ?
              m_v[k] * (m_T_n[k] -  m_T[k]) / m_dy :
              m_v[k] * (m_T[k]  - m_T_s[k]) / m_dy);
    }

    S.RHS(k) += m_dt * (Sigma / m_rho_c_I - UpTu - UpTv);
  }

  // surface b.c.
  if (m_ks>0) {
    S.L(m_ks) = 0.0;
    S.D(m_ks) = 1.0;
    // ignore U[m_ks]
    S.RHS(m_ks) = m_Ts;
  }

  // mark column as done
  m_surfBCsValid = false;
  m_basalBCsValid = false;

  // solve it; note melting not addressed yet
  try {
    S.solve(m_ks + 1, x);
  }
  catch (RuntimeError &e) {
    e.add_context("solving the tri-diagonal system (tempSystemCtx) at (%d,%d)\n"
                  "saving system to m-file... ", m_i, m_j);
    reportColumnZeroPivotErrorMFile(m_ks + 1);
    throw;
  }
}


} // end of namespace energy
} // end of namespace pism
