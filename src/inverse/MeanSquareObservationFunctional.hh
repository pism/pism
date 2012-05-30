// Copyright (C) 2012  David Maxwell
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

#ifndef MEANSQUAREOBSERVATIONFUNCTIONAL_HH_N4PR9LVP
#define MEANSQUAREOBSERVATIONFUNCTIONAL_HH_N4PR9LVP

#include "Functional.hh"

class MeanSquareObservationFunctional2V : public Functional<IceModelVec2V> {
public:
  MeanSquareObservationFunctional2V(IceGrid &grid, IceModelVec2S *weights=NULL) :
  Functional<IceModelVec2V>(grid), m_weights(weights), m_normalization(1.) {};
  virtual ~MeanSquareObservationFunctional2V() {};

  virtual PetscErrorCode normalize();

  virtual PetscErrorCode valueAt(IceModelVec2V &x, PetscReal *OUTPUT);
  virtual PetscErrorCode gradientAt(IceModelVec2V &x, IceModelVec2V &gradient);

protected:
  IceModelVec2S *m_weights;
  PetscReal m_normalization;

private:
  MeanSquareObservationFunctional2V(MeanSquareObservationFunctional2V const &);
  MeanSquareObservationFunctional2V & operator=(MeanSquareObservationFunctional2V const &);
};


#endif /* end of include guard: MEANSQUAREOBSERVATIONFUNCTIONAL_HH_N4PR9LVP */
