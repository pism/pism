/* Copyright (C) 2015 PISM Authors
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

#ifndef _VEC_H_
#define _VEC_H_

#include <petscvec.h>

#include "Wrapper.hh"

namespace pism {
namespace petsc {
/** Wrapper around PETSc's Vec. Simplifies memory management.
 *
 * The constructor takes ownership of the dm argument passed to it.
 *
 * The destructor call VecDestroy().
 */
class Vec : public petsc::Wrapper< ::Vec > {
public:
  Vec(::Vec dm);
  ~Vec();
};
} // end of namespace petsc
} // end of namespace pism


#endif /* _VEC_H_ */
