/* Copyright (C) 2022 PISM Authors
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

#include "IceModelVec2Stag.hh"

#include "pism/util/Mask.hh"

#include "pism/util/pism_utilities.hh"
#include "pism/util/array/CellType.hh"
#include "pism/util/IceModelVec2V.hh"

namespace pism {

IceModelVec2Stag::IceModelVec2Stag(IceGrid::ConstPtr grid, const std::string &name,
                                   IceModelVecKind ghostedp,
                                   unsigned int stencil_width)
  : IceModelVec3(grid, name, ghostedp, 2, stencil_width) {
  set_begin_access_use_dof(true);
}

std::array<double,2> absmax(const IceModelVec2Stag &input) {

  double z[2] = {0.0, 0.0};

  IceModelVec::AccessList list(input);
  for (Points p(*input.grid()); p; p.next()) {
    const int i = p.i(), j = p.j();

    z[0] = std::max(z[0], std::abs(input(i, j, 0)));
    z[1] = std::max(z[1], std::abs(input(i, j, 1)));
  }

  double result[2];
  GlobalMax(input.grid()->com, z, result, 2);

  return {result[0], result[1]};
}

void staggered_to_regular(const array::CellType1 &cell_type,
                          const IceModelVec2Stag &input,
                          bool include_floating_ice,
                          IceModelVec2S &result) {

  using mask::grounded_ice;
  using mask::icy;

  assert(cell_type.stencil_width() > 0);
  assert(input.stencil_width() > 0);

  IceGrid::ConstPtr grid = result.grid();

  IceModelVec::AccessList list{&cell_type, &input, &result};

  for (Points p(*grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    if (cell_type.grounded_ice(i, j) or
        (include_floating_ice and cell_type.icy(i, j))) {
      auto M = cell_type.star(i, j);
      auto F = input.star(i, j);

      int n = 0, e = 0, s = 0, w = 0;
      if (include_floating_ice) {
        n = static_cast<int>(icy(M.n));
        e = static_cast<int>(icy(M.e));
        s = static_cast<int>(icy(M.s));
        w = static_cast<int>(icy(M.w));
      } else {
        n = static_cast<int>(grounded_ice(M.n));
        e = static_cast<int>(grounded_ice(M.e));
        s = static_cast<int>(grounded_ice(M.s));
        w = static_cast<int>(grounded_ice(M.w));
      }

      if (n + e + s + w > 0) {
        result(i, j) = (n * F.n + e * F.e + s * F.s + w * F.w) / (n + e + s + w);
      } else {
        result(i, j) = 0.0;
      }
    } else {
      result(i, j) = 0.0;
    }
  }
}

void staggered_to_regular(const array::CellType1 &cell_type,
                          const IceModelVec2Stag &input,
                          bool include_floating_ice,
                          IceModelVec2V &result) {

  using mask::grounded_ice;
  using mask::icy;

  assert(cell_type.stencil_width() > 0);
  assert(input.stencil_width() > 0);

  IceGrid::ConstPtr grid = result.grid();

  IceModelVec::AccessList list{&cell_type, &input, &result};

  for (Points p(*grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    auto M = cell_type.star(i, j);
    auto F = input.star(i, j);

    int n = 0, e = 0, s = 0, w = 0;
    if (include_floating_ice) {
      n = static_cast<int>(icy(M.n));
      e = static_cast<int>(icy(M.e));
      s = static_cast<int>(icy(M.s));
      w = static_cast<int>(icy(M.w));
    } else {
      n = static_cast<int>(grounded_ice(M.n));
      e = static_cast<int>(grounded_ice(M.e));
      s = static_cast<int>(grounded_ice(M.s));
      w = static_cast<int>(grounded_ice(M.w));
    }

    if (e + w > 0) {
      result(i, j).u = (e * F.e + w * F.w) / (e + w);
    } else {
      result(i, j).u = 0.0;
    }

    if (n + s > 0) {
      result(i, j).v = (n * F.n + s * F.s) / (n + s);
    } else {
      result(i, j).v = 0.0;
    }
  }
}

} // end of namespace pism
