// Copyright (C) 2013  David Maxwell
//
// This file is part of PISM.
//
// PISM is free software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the Free Software
// Foundation; either version 3 of the License, or (at your option) any later
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

#ifndef IPLOGRATIOFUNCTIONAL_HH_HSEWI0Q8
#define IPLOGRATIOFUNCTIONAL_HH_HSEWI0Q8

#include "IPFunctional.hh"

class IPLogRatioFunctional : public IPFunctional<IceModelVec2V> {
public:
  IPLogRatioFunctional(IceGrid &grid, IceModelVec2V &u_observed) :
  IPFunctional<IceModelVec2V>(grid), m_u_observed(u_observed), m_normalization(1.) {};
  virtual ~IPLogRatioFunctional() {};

  virtual PetscErrorCode normalize(PetscReal scale);

  virtual PetscErrorCode valueAt(IceModelVec2V &x, PetscReal *OUTPUT);
  virtual PetscErrorCode gradientAt(IceModelVec2V &x, IceModelVec2V &gradient);

protected:
  IceModelVec2V &m_u_observed;
  PetscReal m_normalization;

};

#endif /* end of include guard: IPLOGRATIOFUNCTIONAL_HH_HSEWI0Q8 */
