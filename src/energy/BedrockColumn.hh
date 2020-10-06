// Copyright (C) 2019 PISM Authors
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

#ifndef BEDROCK_COLUMN_HH
#define BEDROCK_COLUMN_HH

#include "pism/util/ColumnSystem.hh"

namespace pism {

class Config;

namespace energy {

/*! Tridiagonal linear system for conservation of energy in vertical columns of the
 *  bedrock thermal layer.
 *
 * The top surface uses the Dirichlet boundary condition, the bottom surface a heat flux
 * (Neumann) BC.
 *
 * The implementation uses a second-order discretization in space and the backward-Euler
 * (first-order, fully implicit) time-discretization.
 */
class BedrockColumn {
public:
  BedrockColumn(const std::string &prefix, const Config &config,
                double dz, unsigned int M);
  ~BedrockColumn();

  void solve(double dt, double Q_bottom, double T_top,
             const double *T_old, double *result);

  void solve(double dt, double Q_bottom, double T_top,
             const std::vector<double> &T_old,
             std::vector<double> &result);

private:
  // temperature diffusivity coefficient
  double m_D;
  // thermal conductivity
  double m_k;
  // vertical spacing
  double m_dz;
  // system size
  unsigned int m_M;

  TridiagonalSystem m_system;
};

} // end of namespace energy
} // end of namespace pism

#endif   //  ifndef BEDROCK_COLUMN_HH
