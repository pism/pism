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

#ifndef AGECOLUMNSYSTEM_H
#define AGECOLUMNSYSTEM_H

#include "pism/util/ColumnSystem.hh"

namespace pism {

//! Tridiagonal linear system for vertical column of age (pure advection) problem.
class AgeColumnSystem : public columnSystemCtx {
public:
  AgeColumnSystem(const std::vector<double>& storage_grid,
                  const std::string &my_prefix,
                  double dx, double dy, double dt,
                  const array::Array3D &age,
                  const array::Array3D &u3,
                  const array::Array3D &v3,
                  const array::Array3D &w3);

  void init(int i, int j, double thickness);

  void solve(std::vector<double> &x);
protected:
  const array::Array3D &m_age3;
  double m_nu;
  std::vector<double> m_A, m_A_n, m_A_e, m_A_s, m_A_w;
};

} // end of namespace pism


#endif /* AGECOLUMNSYSTEM_H */
