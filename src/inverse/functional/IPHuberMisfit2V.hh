// Copyright (C) 2026 Andy Aschwanden
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

#ifndef PISM_IPHUBERMISFIT2V_H
#define PISM_IPHUBERMISFIT2V_H

#include "pism/inverse/functional/IPFunctional.hh"
#include "pism/util/array/Scalar.hh"
#include "pism/util/array/Vector.hh"

namespace pism {
namespace inverse {

//! Implements a (possibly weighted) Huber-loss misfit functional on a 2D
//! velocity residual.
/*!
  Given a residual field \f$x_i = u_i - u_{\rm obs,i}\f$ with magnitudes
  \f$|x_i|=\sqrt{u^2+v^2}\f$, the functional is
  \f[
  J(x) = c_N \sum_i w_i \,\rho_\delta(|x_i|)
  \f]
  with the Huber kernel
  \f[
  \rho_\delta(r) =
  \begin{cases}
    \tfrac12 r^2                  & r \le \delta \\
    \delta\,(r - \tfrac12\delta)  & r > \delta
  \end{cases}
  \f]
  In the quadratic regime (\f$r\le\delta\f$) the gradient is the residual itself;
  in the linear regime (\f$r>\delta\f$) the gradient saturates at magnitude
  \f$\delta\f$ along the direction of the residual. This caps the influence of
  individual outliers on the inversion's descent direction.

  The normalization \f$c_N\f$ matches IPMeanSquareFunctional2V: a uniform
  residual of magnitude `scale` everywhere (with all entries in the quadratic
  regime) gives \f$J=\tfrac12\f$. As \f$\delta\to\infty\f$ this functional
  reduces to half of IPMeanSquareFunctional2V — the factor of 1/2 is baked
  into the Huber definition, not a bug.

  Note: this is NOT an inner-product functional (Huber is not bilinear). It
  works with TAO Tikhonov solvers (which only need valueAt + gradientAt). It
  does not work with the SSA Gauss-Newton solver.
*/
class IPHuberMisfit2V : public IPFunctional<array::Vector> {
public:
  /*!
   * @param[in] grid the computational grid
   * @param[in] delta the Huber transition threshold, in m s^-1
   * @param[in] weights vector of weights (NULL implies all weights are 1)
   */
  IPHuberMisfit2V(std::shared_ptr<const Grid> grid,
                  double delta,
                  array::Scalar *weights = NULL)
    : IPFunctional<array::Vector>(grid),
      m_delta(delta),
      m_weights(weights),
      m_normalization(1.) {};
  virtual ~IPHuberMisfit2V() {};

  virtual void normalize(double scale);

  virtual void valueAt(array::Vector &x, double *OUTPUT);
  virtual void gradientAt(array::Vector &x, array::Vector &gradient);

protected:
  double         m_delta;
  array::Scalar *m_weights;
  double         m_normalization;

private:
  IPHuberMisfit2V(IPHuberMisfit2V const &);
  IPHuberMisfit2V & operator=(IPHuberMisfit2V const &);
};

} // end of namespace inverse
} // end of namespace pism

#endif /* PISM_IPHUBERMISFIT2V_H */
