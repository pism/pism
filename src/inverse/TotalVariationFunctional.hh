// Copyright (C) 2013  David Maxwell
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

#ifndef TOTALVARIATIONFUNCTIONAL_HH_HKBL1T7I
#define TOTALVARIATIONFUNCTIONAL_HH_HKBL1T7I

#include "Functional.hh"

//! Pseduo total variation functional
/*! \f[
J(u) = c\int_\Omega |\nabla u|^q 
\f]
*/

class TotalVariationFunctional2S : public Functional<IceModelVec2S> {
public:
  TotalVariationFunctional2S(IceGrid &grid, PetscReal c, PetscReal q, IceModelVec2Int *dirichletLocations=NULL);

  virtual PetscErrorCode valueAt(IceModelVec2S &x, PetscReal *OUTPUT);
  virtual PetscErrorCode gradientAt(IceModelVec2S &x, IceModelVec2S &gradient);

protected:

  IceModelVec2Int *m_dirichletIndices;
  PetscReal m_schoofReg; // Regularization parameter.
  PetscReal m_c; // scale parameter.
  PetscReal m_lebesgue_exp;

private:
  // Hide copy/assignment operations
  TotalVariationFunctional2S(TotalVariationFunctional2S const &);
  TotalVariationFunctional2S & operator=(TotalVariationFunctional2S const &);
};

#endif /* end of include guard: TOTALVARIATIONFUNCTIONAL_HH_HKBL1T7I */
