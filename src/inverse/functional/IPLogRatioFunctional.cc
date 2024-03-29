// Copyright (C) 2013, 2014, 2015, 2016, 2017, 2020, 2022, 2023  David Maxwell
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

#include "pism/inverse/functional/IPLogRatioFunctional.hh"
#include "pism/util/Grid.hh"
#include "pism/util/array/Vector.hh"
#include "pism/util/array/Scalar.hh"
#include "pism/util/pism_utilities.hh"

namespace pism {
namespace inverse {

//! Determine the normalization constant for the functional.
/*! Sets the normalization constant \f$c_N\f$ so that
\f[
J(x)=1
\f]
if  \f$|x| = \mathtt{scale}|u_{\rm obs}| \f$ everywhere.
*/
void IPLogRatioFunctional::normalize(double scale) {

  double value = 0;

  double w = 1.0;

  array::AccessScope list(m_u_observed);

  if (m_weights) {
    list.add(*m_weights);
  }

  for (auto p = m_grid->points(); p; p.next()) {
    const int i = p.i(), j = p.j();

    if (m_weights) {
      w = (*m_weights)(i, j);
    }

    Vector2d &u_obs_ij = m_u_observed(i, j);
    double obsMagSq = u_obs_ij.u*u_obs_ij.u + u_obs_ij.v*u_obs_ij.v + m_eps*m_eps;

    double modelMagSq = scale*scale*(u_obs_ij.u*u_obs_ij.u + u_obs_ij.v*u_obs_ij.v) + m_eps*m_eps;

    double v = log(modelMagSq/obsMagSq);
    value += w*v*v;
  }

  m_normalization = GlobalSum(m_grid->com, value);
}

void IPLogRatioFunctional::valueAt(array::Vector &x, double *OUTPUT)  {

  // The value of the objective
  double value = 0;

  double w = 1.;

  array::AccessScope list{&x, &m_u_observed};

  if (m_weights) {
    list.add(*m_weights);
  }

  for (auto p = m_grid->points(); p; p.next()) {
    const int i = p.i(), j = p.j();

    if (m_weights) {
      w = (*m_weights)(i, j);
    }
    Vector2d &x_ij = x(i, j);
    Vector2d &u_obs_ij = m_u_observed(i, j);
    Vector2d u_model_ij = x_ij+u_obs_ij;
    double obsMagSq = u_obs_ij.u*u_obs_ij.u + u_obs_ij.v*u_obs_ij.v + m_eps*m_eps;

    double modelMagSq = (u_model_ij.u*u_model_ij.u + u_model_ij.v*u_model_ij.v)+m_eps*m_eps;
    double v = log(modelMagSq/obsMagSq);
    value += w*v*v;
  }

  value /= m_normalization;

  GlobalSum(m_grid->com, &value, OUTPUT, 1);
}

void IPLogRatioFunctional::gradientAt(array::Vector &x, array::Vector &gradient)  {
  gradient.set(0);

  double w = 1.;

  array::AccessScope list{&x, &gradient, &m_u_observed};

  if (m_weights) {
    list.add(*m_weights);
  }

  for (auto p = m_grid->points(); p; p.next()) {
    const int i = p.i(), j = p.j();

    if (m_weights) {
      w = (*m_weights)(i, j);
    }
    Vector2d &x_ij = x(i, j);
    Vector2d &u_obs_ij = m_u_observed(i, j);
    Vector2d u_model_ij = x_ij+u_obs_ij;

    double obsMagSq = u_obs_ij.u*u_obs_ij.u + u_obs_ij.v*u_obs_ij.v + m_eps*m_eps;
    double modelMagSq = (u_model_ij.u*u_model_ij.u + u_model_ij.v*u_model_ij.v)+m_eps*m_eps;
    double v = log(modelMagSq/obsMagSq);
    double dJdw =  2*w*v/modelMagSq;

    gradient(i, j).u = dJdw*2*u_model_ij.u/m_normalization;
    gradient(i, j).v = dJdw*2*u_model_ij.v/m_normalization;
  }
}

} // end of namespace inverse
} // end of namespace pism
