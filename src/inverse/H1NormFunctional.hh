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

#ifndef H1NORMFUNCTIONAL_HH_TF8AKRNQ
#define H1NORMFUNCTIONAL_HH_TF8AKRNQ

#include "Functional.hh"

class H1NormFunctional2S : public Functional<IceModelVec2S> {
public:
  H1NormFunctional2S(IceGrid &grid, PetscReal cL2, 
      PetscReal cH1, IceModelVec2Int *dirichletLocations=NULL) :
      Functional<IceModelVec2S>(grid),
      m_cL2(cL2), m_cH1(cH1), m_dirichletIndices(dirichletLocations) {};
  virtual ~H1NormFunctional2S() {};
  
  virtual PetscErrorCode valueAt(IceModelVec2S &x, PetscReal *OUTPUT);
  // virtual PetscErrorCode ip(IceModelVec2S &a, IceModelVec2S &b, PetscReal *v);
  virtual PetscErrorCode gradientAt(IceModelVec2S &x, IceModelVec2S &gradient);

protected:

  PetscReal m_cL2, m_cH1;
  IceModelVec2Int *m_dirichletIndices;

private:
  H1NormFunctional2S(H1NormFunctional2S const &);
  H1NormFunctional2S & operator=(H1NormFunctional2S const &);  
};

#endif /* end of include guard: H1NORMFUNCTIONAL_HH_TF8AKRNQ */
