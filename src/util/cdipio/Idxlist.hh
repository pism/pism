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
#ifndef PISM_YAXT_IDXLIST
#define PISM_YAXT_IDXLIST

#include "pism/util/Wrapper.hh"

// We have to include mpi.h *outside* of 'extern "C"'.
#include <mpi.h>
extern "C" {
#include "yaxt.h"
}

namespace pism {
namespace yaxt {

/** Wrapper around YAXT Xt_idxlist. Simplifies memory management.
 *
 * The constructor creates a domain decomposition given the size of the global grid and
 * the size and location of the local sub-domain.
 *
 * The destructor calls xt_idxlist_delete.
 */
class Idxlist : public Wrapper< ::Xt_idxlist > {
public:
  Idxlist(int Mx, int My, int xs, int ys, int xm, int ym, int dof);
  Idxlist(int N);
  ~Idxlist();
};

} // end of namespace yaxt
} // end of namespace pism

#endif /* PISM_YAXT_IDXLIST */
