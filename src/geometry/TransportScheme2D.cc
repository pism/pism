/* Copyright (C) 2026 PISM Authors
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

#include "pism/geometry/TransportScheme.hh"

#include "pism/geometry/MPDATA2.hh"
#include "pism/geometry/UNO.hh"
#include "pism/util/array/CellType.hh"
#include "pism/util/error_handling.hh"

namespace pism {

namespace {

//! First-order upwinding (one MPDATA pass) and MPDATA with `N` passes.
class MPDATA2Scheme : public TransportScheme2D {
public:
  MPDATA2Scheme(std::shared_ptr<const Grid> grid, int N, bool nonoscillatory)
    : m_scheme(grid, N), m_nonoscillatory(nonoscillatory) {
    // empty
  }

  void update(double dt, const array::CellType1 &cell_type, const array::Scalar &x,
              const array::Vector &velocity) {
    m_scheme.update(dt, cell_type, x, velocity, m_nonoscillatory);
  }

  const array::Scalar &x() const {
    return m_scheme.x();
  }

private:
  MPDATA2 m_scheme;
  bool m_nonoscillatory;
};

//! Upstream non-oscillatory schemes UNO2 and UNO3 (Li, 2008).
class UNOScheme : public TransportScheme2D {
public:
  UNOScheme(std::shared_ptr<const Grid> grid, UNOType type) : m_scheme(grid, type) {
    // empty
  }

  void update(double dt, const array::CellType1 &cell_type, const array::Scalar &x,
              const array::Vector &velocity) {
    m_scheme.update(dt, cell_type, x, velocity, true);
  }

  const array::Scalar &x() const {
    return m_scheme.x();
  }

private:
  UNO m_scheme;
};

} // end of anonymous namespace

std::shared_ptr<TransportScheme2D> TransportScheme2D::create(std::shared_ptr<const Grid> grid,
                                                             const std::string &kind, int N,
                                                             bool nonoscillatory) {
  if (kind == "upwind") {
    return std::make_shared<MPDATA2Scheme>(grid, 1, false);
  }

  if (kind == "mpdata") {
    if (N < 1) {
      throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                    "MPDATA requires at least one pass (got %d)", N);
    }
    return std::make_shared<MPDATA2Scheme>(grid, N, nonoscillatory);
  }

  if (kind == "uno2") {
    return std::make_shared<UNOScheme>(grid, PISM_UNO_2);
  }

  if (kind == "uno3") {
    return std::make_shared<UNOScheme>(grid, PISM_UNO_3);
  }

  throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                "unknown 2D transport scheme '%s'"
                                " (allowed: upwind, mpdata, uno2, uno3)",
                                kind.c_str());
}

} // end of namespace pism
