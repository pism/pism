/* Copyright (C) 2014 PISM Authors
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

#ifndef _SIAFD_FEVOR_H_
#define _SIAFD_FEVOR_H_

#include "SIAFD.hh"

namespace pism {

class SIAFD_FEvoR : public SIAFD {
public:
  SIAFD_FEvoR(IceGrid &g, EnthalpyConverter &e, const Config &c);
  virtual ~SIAFD_FEvoR();
  virtual PetscErrorCode init(Vars &vars);
protected:
  virtual PetscErrorCode compute_diffusive_flux(IceModelVec2Stag &h_x, IceModelVec2Stag &h_y,
                                                IceModelVec2Stag &result, bool fast);
private:
  IceModelVec3 *enhancement_factor;
};

} // end of namespace pism

#endif /* _SIAFD_FEVOR_H_ */
