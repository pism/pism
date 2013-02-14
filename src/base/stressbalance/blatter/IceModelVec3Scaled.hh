/* Copyright (C) 2013 PISM Authors
 *
 * This file is part of PISM.
 *
 * PISM is free software; you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation; either version 2 of the License, or (at your option) any later
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

#ifndef _ICEMODELVEC3SCALED_H_
#define _ICEMODELVEC3SCALED_H_

#include "iceModelVec.hh"

//! Class to reading and writing 3D horizontal velocity components on
//! the "sigma-coordinate" vertical grid.
class IceModelVec3Scaled : public IceModelVec3D {
public:
  IceModelVec3Scaled() {}
  virtual ~IceModelVec3Scaled() {}

  virtual PetscErrorCode create(IceGrid &mygrid, const char my_short_name[], bool local,
                                int myMbz, int stencil_width = 1);
};


#endif /* _ICEMODELVEC3SCALED_H_ */
