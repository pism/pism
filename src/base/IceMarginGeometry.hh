// Copyright (C) 2004-2011 Jed Brown, Ed Bueler, Constantine Khroulev
//                         and David Maxwell
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

#ifndef __ICEMARGINGEOMETRY_HH
#define __ICEMARGINGEOMETRY_HH

#include "flowlaws.hh"
#include "NCVariable.hh"
#include "pism_const.hh"

class IceMarginGeometry
{
public:
  IceMarginGeometry( PetscReal seaLevel, const IceFlowLaw &ice, const NCConfigVariable &config);

  inline void computeGeometry(PetscReal bed, PetscReal thickness,
    PetscInt *mask, PetscReal *surface);

protected:
  PetscReal m_alpha;
  PetscReal m_seaLevel;
};

inline void IceMarginGeometry::computeGeometry(PetscReal bed, PetscReal thickness,
                                        PetscInt *mask, PetscReal *surface) {
  const PetscReal  hgrounded = bed + thickness; // FIXME task #7297
  const PetscReal  hfloating = m_seaLevel + m_alpha*thickness;

  const bool is_floating = hfloating>hgrounded + 1.0,
    ice_free = thickness < 0.01;

  if (is_floating) {
    *surface = hfloating;
    if (ice_free) {
      *mask = MASK_ICE_FREE_OCEAN;
        // added just for clarity, needs to be tested for interference in run
    } else {
      *mask = MASK_FLOATING; // to enable for floating front propagation
    }
  } else {  // Grounded
    *surface = hgrounded;
    if (ice_free) {
      *mask = MASK_ICE_FREE_BEDROCK;
    } else {
      *mask = MASK_GROUNDED;
    }
  }
}


#endif