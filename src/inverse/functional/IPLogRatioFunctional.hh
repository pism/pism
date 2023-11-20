// Copyright (C) 2013, 2014, 2015, 2022  David Maxwell and Constantine Khroulev
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

#include "pism/inverse/functional/IPFunctional.hh"

namespace pism {
namespace inverse {

//! Implements a functional for log-ratio errors.
/*!  This type of functional appears in [\ref Morlighemetal2010].
  Specifically, given a reference function \f$u_{obs}=[U_i]\f$, and an
  array::Vector \f$x=[X_i]\f$,
  \f[
  J(x) = c_N \sum_i W_i\left[\log\left(\frac{|X_i+U_i|^2+\epsilon^2}{|U_{i}|^2+\epsilon^2}\right)\right]^2
  \f]
  where \f$\epsilon\f$ is a regularizing constant and \f$[W_i]\f$ is a vector of weights.
  The term \f$X_i+U_i\f$ appears because the argument is expected to already be in the form
  \f$V_i-U_i\f$, where \f$v=[V_i]\f$ is some approximation of \f$[U_i]\f$ and hence the
  integrand has the form \f$\log(|V_i|/|U_i|)\f$.

  The normalization constant \f$c_N\f$ is determined implicitly by normalize().
*/
class IPLogRatioFunctional : public IPFunctional<array::Vector> {
public:
  IPLogRatioFunctional(std::shared_ptr<const Grid> grid, array::Vector &u_observed, double eps,
                       array::Scalar *weights=NULL) :
    IPFunctional<array::Vector>(grid), m_u_observed(u_observed), m_weights(weights),
    m_normalization(1.), m_eps(eps) {};
  virtual ~IPLogRatioFunctional() {};

  virtual void normalize(double scale);

  virtual void valueAt(array::Vector &x, double *OUTPUT);
  virtual void gradientAt(array::Vector &x, array::Vector &gradient);

protected:
  array::Vector &m_u_observed;
  array::Scalar *m_weights;
  double m_normalization;
  double m_eps;

};

} // end of namespace inverse
} // end of namespace pism

#endif /* end of include guard: IPLOGRATIOFUNCTIONAL_HH_HSEWI0Q8 */
