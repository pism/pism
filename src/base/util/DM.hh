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

#ifndef _DM_H_
#define _DM_H_

#include <petscdmda.h>

#ifdef PISM_USE_TR1
#include <tr1/memory>
#else
#include <memory>
#endif

namespace pism {
/** Wrapper around PETSc's DM. Simplifies memory management.
 *
 * The constructor takes ownership of the dm argument passed to it.
 *
 * The destructor call DMDestroy().
 */
class PISMDM {
public:
#ifdef PISM_USE_TR1
  typedef std::tr1::shared_ptr<PISMDM> Ptr;
  typedef std::tr1::weak_ptr<PISMDM> WeakPtr;
#else
  typedef std::shared_ptr<PISMDM> Ptr;
  typedef std::weak_ptr<PISMDM> WeakPtr;
#endif
  PISMDM(DM dm);
  ~PISMDM();
  DM get() const;
  operator DM() const;
private:
  DM m_dm;
};

} // end of namespace pism

#endif /* _DM_H_ */
