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

#ifndef FUNCTIONAL_HH_1E2DIXE6
#define FUNCTIONAL_HH_1E2DIXE6

#include "iceModelVec.hh"
#include "FETools.hh"

template<class IMVecType>
class Functional {
public:
  
  Functional(IceGrid &grid) : m_grid(grid), m_element_index(m_grid) { 
    m_quadrature.init(m_grid);
  }

  virtual ~Functional() {};

  virtual PetscErrorCode valueAt(IMVecType &x, PetscReal *OUTPUT) = 0;
  virtual PetscErrorCode gradientAt(IMVecType &x, IMVecType &gradient) = 0;

protected:
  IceGrid &m_grid;

  FEElementMap m_element_index;
  FEQuadrature m_quadrature;
  FEDOFMap     m_dofmap;
  
private:
  // Hide copy/assignment operations
  Functional(Functional const &);
  Functional & operator=(Functional const &);

};

template<class IMVecType>
class IPFunctional : public Functional<IMVecType>{
public:
  IPFunctional(IceGrid &grid) : Functional<IMVecType>(grid) {};
  virtual PetscErrorCode dot(IMVecType &a, IMVecType &b, PetscReal *OUTPUT) = 0;
};

PetscErrorCode gradientFD(Functional<IceModelVec2S> &f, IceModelVec2S &x, IceModelVec2S &gradient);
PetscErrorCode gradientFD(Functional<IceModelVec2V> &f, IceModelVec2V &x, IceModelVec2V &gradient);

#endif /* end of include guard: FUNCTIONAL_HH_1E2DIXE6 */
