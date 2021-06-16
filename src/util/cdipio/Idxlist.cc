/* Copyright (C) 2021 PISM Authors
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

#include "Idxlist.hh"

#include <vector>

namespace pism {
namespace yaxt {

Idxlist::Idxlist(int Mx, int My, int xs, int ys, int xm, int ym, int dof) {
  // domain decomposition for transposed data
  std::vector<int> gdimlen{dof, My, Mx};
  std::vector<int> start{0, ys, xs};
  std::vector<int> count{dof, ym, xm};
  int idxlen = count[0] * count[1] * count[2];
  std::vector<Xt_int> idx(idxlen);

  for (int z = 0; z < count[0]; z++) {
    for (int y = 0; y < count[1]; y++) {
      for (int x = 0; x < count[2]; x++) {
        idx[(z * count[1] + y) * count[2] + x] =
          ((z + start[0]) * gdimlen[1] + (y + start[1])) * gdimlen[2] + (x + start[2]);
      }
    }
  }

  m_value = xt_idxvec_new(idx.data(), idxlen);
}

Idxlist::~Idxlist() {
  if (m_value != NULL) {
    xt_idxlist_delete(m_value);
  }
}

} // end of namespace yaxt
} // end of namespace pism
