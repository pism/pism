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

#ifndef IPMEANSQUAREFUNCTIONAL_HH_DZ18EO5C
#define IPMEANSQUAREFUNCTIONAL_HH_DZ18EO5C


#include "IPFunctional.hh"

class IPMeanSquareFunctional2S : public IPInnerProductFunctional<IceModelVec2S> {
public:
  IPMeanSquareFunctional2S(IceGrid &grid, IceModelVec2S *weights=NULL) :
  IPInnerProductFunctional<IceModelVec2S>(grid), m_weights(weights), m_normalization(1.) {};
  virtual ~IPMeanSquareFunctional2S() {};

  virtual PetscErrorCode normalize(PetscReal scale);

  virtual PetscErrorCode valueAt(IceModelVec2S &x, PetscReal *OUTPUT);
  virtual PetscErrorCode dot(IceModelVec2S &a, IceModelVec2S &b, PetscReal *OUTPUT);
  virtual PetscErrorCode gradientAt(IceModelVec2S &x, IceModelVec2S &gradient);

protected:
  IceModelVec2S *m_weights;
  PetscReal m_normalization;

private:
  IPMeanSquareFunctional2S(IPMeanSquareFunctional2S const &);
  IPMeanSquareFunctional2S & operator=(IPMeanSquareFunctional2S const &);
};


class IPMeanSquareFunctional2V : public IPInnerProductFunctional<IceModelVec2V> {
public:
  IPMeanSquareFunctional2V(IceGrid &grid, IceModelVec2S *weights=NULL) :
  IPInnerProductFunctional<IceModelVec2V>(grid), m_weights(weights), m_normalization(1.) {};
  virtual ~IPMeanSquareFunctional2V() {};

  virtual PetscErrorCode normalize(PetscReal scale);

  virtual PetscErrorCode valueAt(IceModelVec2V &x, PetscReal *OUTPUT);
  virtual PetscErrorCode dot(IceModelVec2V &a, IceModelVec2V &b, PetscReal *OUTPUT);
  virtual PetscErrorCode gradientAt(IceModelVec2V &x, IceModelVec2V &gradient);

protected:
  IceModelVec2S *m_weights;
  PetscReal m_normalization;

private:
  IPMeanSquareFunctional2V(IPMeanSquareFunctional2V const &);
  IPMeanSquareFunctional2V & operator=(IPMeanSquareFunctional2V const &);
};


#endif /* end of include guard: IPMEANSQUAREFUNCTIONAL_HH_DZ18EO5C */
