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

//! Implements a functional for log-ratio errors.
/*!  This type of functional appears in [\ref Morlighemetal2010].
Specifically, given a reference function \f$u_{obs}=[U_i]\f$, and an
IceModelVec2V \f$x=[X_i]\f$,
\f[
J(x) = c_N \sum_i \log\left(\frac{|X_i+U_i|+\epsilon}{|U_{i}|+\epsilon}\right)
\f]
where \f$\epsilon=10^{-4}{\tt inv_ssa_velocity_scale}\f$.  The term \f$X_i+U_i\f$
appears because the argument is expected to already be in the form \f$X_i-U_i\f$.

The normalization constant \f$c_N\f$ is determined implicitly by ::normalize.
*/
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
