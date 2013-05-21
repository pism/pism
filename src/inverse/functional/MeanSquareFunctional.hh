// Copyright (C) 2012  David Maxwell
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

#ifndef MEANSQUAREFUNCTIONAL_HH_DZ18EO5C
#define MEANSQUAREFUNCTIONAL_HH_DZ18EO5C


#include "Functional.hh"

class MeanSquareFunctional2S : public IPFunctional<IceModelVec2S> {
public:
  MeanSquareFunctional2S(IceGrid &grid, IceModelVec2S *weights=NULL) :
  IPFunctional<IceModelVec2S>(grid), m_weights(weights), m_normalization(1.) {};
  virtual ~MeanSquareFunctional2S() {};

  virtual PetscErrorCode normalize(PetscReal scale);

  virtual PetscErrorCode valueAt(IceModelVec2S &x, PetscReal *OUTPUT);
  virtual PetscErrorCode dot(IceModelVec2S &a, IceModelVec2S &b, PetscReal *OUTPUT);
  virtual PetscErrorCode gradientAt(IceModelVec2S &x, IceModelVec2S &gradient);

protected:
  IceModelVec2S *m_weights;
  PetscReal m_normalization;

private:
  MeanSquareFunctional2S(MeanSquareFunctional2S const &);
  MeanSquareFunctional2S & operator=(MeanSquareFunctional2S const &);
};


class MeanSquareFunctional2V : public IPFunctional<IceModelVec2V> {
public:
  MeanSquareFunctional2V(IceGrid &grid, IceModelVec2S *weights=NULL) :
  IPFunctional<IceModelVec2V>(grid), m_weights(weights), m_normalization(1.) {};
  virtual ~MeanSquareFunctional2V() {};

  virtual PetscErrorCode normalize(PetscReal scale);

  virtual PetscErrorCode valueAt(IceModelVec2V &x, PetscReal *OUTPUT);
  virtual PetscErrorCode dot(IceModelVec2V &a, IceModelVec2V &b, PetscReal *OUTPUT);
  virtual PetscErrorCode gradientAt(IceModelVec2V &x, IceModelVec2V &gradient);

protected:
  IceModelVec2S *m_weights;
  PetscReal m_normalization;

private:
  MeanSquareFunctional2V(MeanSquareFunctional2V const &);
  MeanSquareFunctional2V & operator=(MeanSquareFunctional2V const &);
};


#endif /* end of include guard: MEANSQUAREFUNCTIONAL_HH_DZ18EO5C */
