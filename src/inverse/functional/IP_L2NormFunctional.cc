// Copyright (C) 2012, 2014, 2015, 2016, 2017  David Maxwell and Constantine Khroulev
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

#include "IP_L2NormFunctional.hh"
#include "pism/util/IceGrid.hh"
#include "pism/util/pism_utilities.hh"

namespace pism {
namespace inverse {

void IP_L2NormFunctional2S::valueAt(IceModelVec2S &x, double *OUTPUT) {

  const unsigned int Nq     = m_quadrature.n();
  const unsigned int Nq_max = fem::MAX_QUADRATURE_SIZE;

  // The value of the objective
  double value = 0;

  double x_q[Nq_max];

  IceModelVec::AccessList list(x);

  // Jacobian times weights for quadrature.
  const double* W = m_quadrature.weights();

  // Loop through all LOCAL elements.
  const int
    xs = m_element_index.lxs,
    xm = m_element_index.lxm,
    ys = m_element_index.lys,
    ym = m_element_index.lym;

  for (int j = ys; j < ys + ym; j++) {
    for (int i = xs; i < xs + xm; i++) {
      m_element.reset(i, j);

      // Obtain values of x at the quadrature points for the element.
      double tmp[fem::q1::n_chi];
      m_element.nodal_values(x, tmp);
      quadrature_point_values(m_quadrature, tmp, x_q);

      for (unsigned int q = 0; q < Nq; q++) {
        const double x_qq = x_q[q];
        value += W[q]*x_qq*x_qq;
      } // q
    } // j
  } // i

  GlobalSum(m_grid->com, &value, OUTPUT, 1);
}

void IP_L2NormFunctional2S::dot(IceModelVec2S &a, IceModelVec2S &b, double *OUTPUT) {

  const unsigned int Nq     = m_quadrature.n();
  const unsigned int Nq_max = fem::MAX_QUADRATURE_SIZE;

  // The value of the objective
  double value = 0;

  double a_q[Nq_max];
  double b_q[Nq_max];

  IceModelVec::AccessList list{&a, &b};

  // Jacobian times weights for quadrature.
  const double* W = m_quadrature.weights();

  // Loop through all LOCAL elements.
  const int
    xs = m_element_index.lxs,
    xm = m_element_index.lxm,
    ys = m_element_index.lys,
    ym = m_element_index.lym;

  for (int j = ys; j < ys + ym; j++) {
    for (int i = xs; i < xs + xm; i++) {
      m_element.reset(i, j);

      double tmp[fem::q1::n_chi];
      m_element.nodal_values(a, tmp);
      quadrature_point_values(m_quadrature, tmp, a_q);

      m_element.nodal_values(b, tmp);
      quadrature_point_values(m_quadrature, tmp, b_q);

      for (unsigned int q = 0; q < Nq; q++) {
        value += W[q]*a_q[q]*b_q[q];
      } // q
    } // j
  } // i

  GlobalSum(m_grid->com, &value, OUTPUT, 1);
}

void IP_L2NormFunctional2S::gradientAt(IceModelVec2S &x, IceModelVec2S &gradient) {

  const unsigned int Nk     = fem::q1::n_chi;
  const unsigned int Nq     = m_quadrature.n();
  const unsigned int Nq_max = fem::MAX_QUADRATURE_SIZE;

  // Clear the gradient before doing anything with it!
  gradient.set(0);

  double x_q[Nq_max];
  double gradient_e[Nk];

  IceModelVec::AccessList list{&x, &gradient};

  // An Nq by Nk array of test function values.
  const fem::Germs *test = m_quadrature.test_function_values();

  // Jacobian times weights for quadrature.
  const double* W = m_quadrature.weights();

  // Loop through all local and ghosted elements.
  const int
    xs = m_element_index.xs,
    xm = m_element_index.xm,
    ys = m_element_index.ys,
    ym = m_element_index.ym;

  for (int j = ys; j < ys + ym; j++) {
    for (int i = xs; i < xs + xm; i++) {

      // Reset the DOF map for this element.
      m_element.reset(i, j);

      // Obtain values of x at the quadrature points for the element.
      double tmp[Nk];
      m_element.nodal_values(x, tmp);
      quadrature_point_values(m_quadrature, tmp, x_q);

      // Zero out the element-local residual in prep for updating it.
      for (unsigned int k = 0; k < Nk; k++) {
        gradient_e[k] = 0;
      }

      for (unsigned int q = 0; q < Nq; q++) {
        const double x_qq = x_q[q];
        for (unsigned int k = 0; k < Nk; k++) {
          gradient_e[k] += 2*W[q]*x_qq*test[q][k].val;
        } // k
      } // q
      m_element.add_contribution(gradient_e, gradient);
    } // j
  } // i
}

void IP_L2NormFunctional2V::valueAt(IceModelVec2V &x, double *OUTPUT) {

  const unsigned int Nq     = m_quadrature.n();
  const unsigned int Nq_max = fem::MAX_QUADRATURE_SIZE;

  // The value of the objective
  double value = 0;

  Vector2 x_q[Nq_max];

  IceModelVec::AccessList list(x);

  // Jacobian times weights for quadrature.
  const double* W = m_quadrature.weights();

  // Loop through all local and ghosted elements.
  const int
    xs = m_element_index.lxs,
    xm = m_element_index.lxm,
    ys = m_element_index.lys,
    ym = m_element_index.lym;

  for (int j = ys; j < ys + ym; j++) {
    for (int i = xs; i < xs + xm; i++) {
      m_element.reset(i, j);

      // Obtain values of x at the quadrature points for the element.
      Vector2 tmp[fem::q1::n_chi];
      m_element.nodal_values(x, tmp);
      quadrature_point_values(m_quadrature, tmp, x_q);

      for (unsigned int q = 0; q < Nq; q++) {
        const Vector2 &x_qq = x_q[q];
        value += W[q]*(x_qq.u*x_qq.u + x_qq.v*x_qq.v);
      } // q
    } // j
  } // i

  GlobalSum(m_grid->com, &value, OUTPUT, 1);
}

void IP_L2NormFunctional2V::dot(IceModelVec2V &a, IceModelVec2V &b, double *OUTPUT) {

  const unsigned int Nq     = m_quadrature.n();
  const unsigned int Nq_max = fem::MAX_QUADRATURE_SIZE;

  // The value of the objective
  double value = 0;

  Vector2 a_q[Nq_max];
  Vector2 b_q[Nq_max];

  IceModelVec::AccessList list{&a, &b};

  // Jacobian times weights for quadrature.
  const double* W = m_quadrature.weights();

  // Loop through all LOCAL elements.
  const int
    xs = m_element_index.lxs,
    xm = m_element_index.lxm,
    ys = m_element_index.lys,
    ym = m_element_index.lym;

  for (int j = ys; j < ys + ym; j++) {
    for (int i = xs; i < xs + xm; i++) {
      m_element.reset(i, j);

      // Obtain values of x at the quadrature points for the element.
      Vector2 tmp[fem::q1::n_chi];
      m_element.nodal_values(a, tmp);
      quadrature_point_values(m_quadrature, tmp, a_q);
      m_element.nodal_values(b, tmp);
      quadrature_point_values(m_quadrature, tmp, b_q);

      for (unsigned int q = 0; q < Nq; q++) {
        value += W[q]*(a_q[q].u*b_q[q].u + a_q[q].v*b_q[q].v);
      } // q
    } // j
  } // i

  GlobalSum(m_grid->com, &value, OUTPUT, 1);
}

void IP_L2NormFunctional2V::gradientAt(IceModelVec2V &x, IceModelVec2V &gradient) {

  const unsigned int Nk     = fem::q1::n_chi;
  const unsigned int Nq     = m_quadrature.n();
  const unsigned int Nq_max = fem::MAX_QUADRATURE_SIZE;

  // Clear the gradient before doing anything with it!
  gradient.set(0);

  Vector2 x_q[Nq_max];
  Vector2 gradient_e[Nk];

  IceModelVec::AccessList list{&x, &gradient};

  // An Nq by Nk array of test function values.
  const fem::Germs *test = m_quadrature.test_function_values();

  // Jacobian times weights for quadrature.
  const double* W = m_quadrature.weights();

  // Loop through all local and ghosted elements.
  const int
    xs = m_element_index.xs,
    xm = m_element_index.xm,
    ys = m_element_index.ys,
    ym = m_element_index.ym;

  for (int j = ys; j < ys + ym; j++) {
    for (int i = xs; i < xs + xm; i++) {

      // Reset the DOF map for this element.
      m_element.reset(i, j);

      // Obtain values of x at the quadrature points for the element.
      Vector2 tmp[Nk];
      m_element.nodal_values(x, tmp);
      quadrature_point_values(m_quadrature, tmp, x_q);

      // Zero out the element-local residual in prep for updating it.
      for (unsigned int k = 0; k < Nk; k++) {
        gradient_e[k].u = 0;
        gradient_e[k].v = 0;
      }

      for (unsigned int q = 0; q < Nq; q++) {
        const Vector2 &x_qq = x_q[q];
        for (unsigned int k = 0; k < Nk; k++) {
          double gcommon = 2*W[q]*test[q][k].val;
          gradient_e[k].u += gcommon*x_qq.u;
          gradient_e[k].v += gcommon*x_qq.v;
        } // k
      } // q
      m_element.add_contribution(gradient_e, gradient);
    } // j
  } // i
}

} // end of namespace inverse
} // end of namespace pism
