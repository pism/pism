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

#ifndef TOTALVARIATIONFUNCTIONAL_HH_HKBL1T7I
#define TOTALVARIATIONFUNCTIONAL_HH_HKBL1T7I

#include "IPFunctional.hh"

//! Pseduo total variation functional
/*! \f[
J(u) = c\int_\Omega (\epsilon^2+|\nabla u|^2)^{q/2} 
\f]
The parameters \f$c\f$ and \f$q\f$ are provided at construction.  Taking \f$q\f$=1 would
yield a total variation functional, save for the regularizing parameter \f$\epsilon\f$
which is determined from config parameters:
\f[
\epsilon= \frac{\tt Schoof\_regularizing\_velocity}{\tt Schoof\_regularizing\_length}.
\f]
*/
class IPTotalVariationFunctional2S : public IPFunctional<IceModelVec2S> {
public:
  IPTotalVariationFunctional2S(IceGrid &grid, PetscReal c, PetscReal q, IceModelVec2Int *dirichletLocations=NULL);

  virtual PetscErrorCode valueAt(IceModelVec2S &x, PetscReal *OUTPUT);
  virtual PetscErrorCode gradientAt(IceModelVec2S &x, IceModelVec2S &gradient);

protected:

  IceModelVec2Int *m_dirichletIndices;
  PetscReal m_schoofReg; // Regularization parameter.
  PetscReal m_c; // scale parameter.
  PetscReal m_lebesgue_exp;

private:
  // Hide copy/assignment operations
  IPTotalVariationFunctional2S(IPTotalVariationFunctional2S const &);
  IPTotalVariationFunctional2S & operator=(IPTotalVariationFunctional2S const &);
};

#endif /* end of include guard: TOTALVARIATIONFUNCTIONAL_HH_HKBL1T7I */
