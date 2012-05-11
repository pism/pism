// Copyright (C) 2012 PISM Authors
//
// This file is part of PISM.
//
// PISM is free software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the Free Software
// Foundation; either version 2 of the License, or (at your option) any later
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

#ifndef _MISMIPBASALRESISTANCELAW_H_
#define _MISMIPBASALRESISTANCELAW_H_

#include "basal_resistance.hh"

class MISMIPBasalResistanceLaw : public IceBasalResistancePlasticLaw
{
public:
  MISMIPBasalResistanceLaw(PetscReal m, PetscReal c, PetscReal r)
    : IceBasalResistancePlasticLaw(1, false, 1, 0)
  {
    m_MISMIP = m; // power
    C_MISMIP = c; // Pa m^(−1/3) s^(1/3)
    regularize_MISMIP = r;
  }
  virtual ~MISMIPBasalResistanceLaw() {}
  virtual PetscErrorCode printInfo(int verbthresh, MPI_Comm com);
  virtual PetscScalar drag(PetscScalar tauc,
                           PetscScalar vx, PetscScalar vy);
  PetscReal m_MISMIP, C_MISMIP, regularize_MISMIP;
};

#endif /* _MISMIPBASALRESISTANCELAW_H_ */
