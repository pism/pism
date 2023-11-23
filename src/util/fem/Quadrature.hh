/* Copyright (C) 2020 PISM Authors
 *
 * This file is part of PISM.
 *
 * PISM is free software; you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation; either version 3 of the License, or (at your option) any later
 * version.
 *
 * PISM is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License
 * along with PISM; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 */
#ifndef PISM_QUADRATURE_H
#define PISM_QUADRATURE_H

#include "pism/util/fem/FEM.hh"

namespace pism {
namespace fem {

//! Numerical integration of finite element functions.
/*! The core of the finite element method is the computation of integrals over elements.
  For nonlinear problems, or problems with non-constant coefficients (%i.e. any real problem)
  the integration has to be done approximately:
  \f[
  \int_E f(x)\; dx \approx \sum_q f(x_q) w_q
  \f]
  for certain quadrature points \f$x_q\f$ and weights \f$w_q\f$.  A quadrature is used
  to evaluate finite element functions at quadrature points, and to compute weights \f$w_q\f$
  for a given element.

  In this concrete implementation, the reference element \f$R\f$ is the square
  \f$[-1,1]\times[-1,1]\f$.  On a given element, nodes (o) and quadrature points (*)
  are ordered as follows:

  ~~~
  3 o------------------o  2
    |  3             2 |
    |    *        *    |
    |                  |
    |                  |
    |    *        *    |
    |  0            1  |
  0 o------------------o  1
  ~~~

  There are four quad points per element, which occur at \f$x,y=\pm 1/\sqrt{3}\f$. This corresponds
  to the tensor product of Gaussian integration on an interval that is exact for cubic functions on
  the interval.

  Integration on a physical element can be thought of as being done by change of variables. The
  quadrature weights need to be modified, and the Quadrature takes care of this. Because all
  elements in an Grid are congruent, the quadrature weights are the same for each element, and
  are computed upon initialization.
*/
class Quadrature {
public:
  const std::vector<QuadPoint>& points() const;
  const std::vector<double>& weights() const;

  QuadPoint point(int k) const {
    return m_points[k];
  }

  double weight(int k) const {
    return m_weights[k];
  }
protected:
  std::vector<QuadPoint> m_points;
  std::vector<double> m_weights;
};


/*!
 * 2-point Gaussian quadrature on an interval of length D.
 */
class Gaussian2 : public Quadrature {
public:
  Gaussian2(double D);
};

//! The 1-point Gaussian quadrature on the square [-1,1]*[-1,1]
class Q1Quadrature1 : public Quadrature {
public:
  Q1Quadrature1();
};

//! The 4-point Gaussian quadrature on the square [-1,1]*[-1,1]
class Q1Quadrature4 : public Quadrature {
public:
  Q1Quadrature4();
};

//! The 9-point Gaussian quadrature on the square [-1,1]*[-1,1]
class Q1Quadrature9 : public Quadrature {
public:
  Q1Quadrature9();
};

//! The 16-point Gaussian quadrature on the square [-1,1]*[-1,1]
class Q1Quadrature16 : public Quadrature {
public:
  Q1Quadrature16();
};

/*!
 * N*N point (NOT Gaussian) quadrature on the square [-1,1]*[-1,1]
 */
class Q1QuadratureN : public Quadrature {
public:
  Q1QuadratureN(unsigned int n);
};

/*!
 * 3-point Gaussian quadrature on the triangle (0,0)-(1,0)-(0,1)
 */
class P1Quadrature3 : public Quadrature {
public:
  P1Quadrature3();
};

/*!
 * 8-point Gaussian quadrature on the cube [-1,1]*[-1,1]*[-1,1]
 */
class Q13DQuadrature8 : public Quadrature {
public:
  Q13DQuadrature8();
};

/*!
 * 1-point Gaussian quadrature on the cube [-1,1]*[-1,1]*[-1,1]
 */
class Q13DQuadrature1 : public Quadrature {
public:
  Q13DQuadrature1();
};

/*!
 * 64-point (4*4*4) Gaussian quadrature on the cube [-1,1]*[-1,1]*[-1,1]
 */
class Q13DQuadrature64 : public Quadrature {
public:
  Q13DQuadrature64();
};
} // end of namespace fem
} // end of namespace pism

#endif /* PISM_QUADRATURE_H */
