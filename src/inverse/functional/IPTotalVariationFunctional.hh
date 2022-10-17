// Copyright (C) 2013, 2014, 2015  David Maxwell and Constantine Khroulev
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

namespace pism {
namespace inverse {

//! Pseduo total variation functional
/*! \f[
  J(u) = c\int_\Omega (\epsilon^2+|\nabla u|^2)^{q/2} 
  \f]
  The parameters \f$c\f$, \f$q\f$ and \f$\epsilon\f$ are provided at construction.  Taking \f$q\f$=1 would
  yield a total variation functional, save for the regularizing parameter \f$\epsilon\f$.
*/
class IPTotalVariationFunctional2S : public IPFunctional<array::Scalar> {
public:
  IPTotalVariationFunctional2S(IceGrid::ConstPtr grid, double c, double q, double eps, array::Scalar *dirichletLocations=NULL);

  virtual void valueAt(array::Scalar &x, double *OUTPUT);
  virtual void gradientAt(array::Scalar &x, array::Scalar &gradient);

protected:

  array::Scalar *m_dirichletIndices;
  double m_c; // scale parameter.
  double m_lebesgue_exp;
  double m_epsilon_sq; // Regularization parameter.

private:
  // Hide copy/assignment operations
  IPTotalVariationFunctional2S(IPTotalVariationFunctional2S const &);
  IPTotalVariationFunctional2S & operator=(IPTotalVariationFunctional2S const &);
};

} // end of namespace inverse
} // end of namespace pism

#endif /* end of include guard: TOTALVARIATIONFUNCTIONAL_HH_HKBL1T7I */
