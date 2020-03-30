// Copyright (C) 2012, 2013, 2014, 2015, 2016, 2017, 2020  David Maxwell
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

#ifndef IPFUNCTIONAL_HH_1E2DIXE6
#define IPFUNCTIONAL_HH_1E2DIXE6

#include "pism/util/iceModelVec.hh"
#include "pism/util/fem/FEM.hh"

namespace pism {

//! Inverse modeling code.
namespace inverse {
//! Abstract base class for functions from ice model vectors to \f$\mathbb{R}\f$.
/*! Inverse problems frequently involve minimizing a functional,
  such such as the misfit
  \f[
  J(u) = \int_\Omega |u-u_{\rm obs}|^2
  \f]
  between observed (\f$u_{\rm obs}\f$) and modeled (\f$u\f$)
  surface velocities. Subclasses of IPFunctional define such maps,
  and permit computation of their gradients.
*/
template<class IMVecType>
class IPFunctional {

public:
  IPFunctional(IceGrid::ConstPtr grid)
    : m_grid(grid),
      m_element_index(*m_grid),
      m_element(*m_grid, fem::Q1Quadrature4())
  {
    // empty
  }

  virtual ~IPFunctional() {};

  //! Computes the value of the functional at the vector x.
  virtual void valueAt(IMVecType &x, double *OUTPUT) = 0;

  //! Computes the gradient of the functional at the vector x.
  /*! On an \f$m\times n\f$ IceGrid, an IceModelVec \f$x\f$ with \f$d\f$
    degrees of freedom will be \f$d m n\f$-dimensional with components \f$x_i\f$.
    The gradient computed here is the vector of directional derivatives \f$\nabla J\f$ of the functional
    \f$J\f$ with respect to \f$x\f$. Concretely, the \f$i^{\rm th}\f$ component of \f$\nabla J\f$
    is
    \f[
    \nabla J_i = \frac{\partial}{\partial x_i} J(x).
    \f]
    This vector is returned as `gradient`.
  */
  virtual void gradientAt(IMVecType &x, IMVecType &gradient) = 0;

protected:
  IceGrid::ConstPtr m_grid;

  fem::ElementIterator m_element_index;
  fem::Q1Element2       m_element;

private:
  // Hide copy/assignment operations
  IPFunctional(IPFunctional const &);
  IPFunctional & operator=(IPFunctional const &);

};

//! Abstract base class for IPFunctionals arising from an inner product.
/*!
  Frequently functionals have the structure
  \f[
  J(u) = Q(u, u)
  \f]
  where \f$Q\f$ is a symmetric, non-negative definite, bilinear form. Certain
  minimization algorithms only apply to such functionals, which should subclass
  from IPInnerProductFunctional.
*/
template<class IMVecType>
class IPInnerProductFunctional : public IPFunctional<IMVecType> {

public:
  IPInnerProductFunctional(IceGrid::ConstPtr grid) : IPFunctional<IMVecType>(grid) {};

  //! Computes the inner product \f$Q(a, b)\f$.
  virtual void dot(IMVecType &a, IMVecType &b, double *OUTPUT) = 0;

  //! Computes the interior product of a vector with the IPInnerProductFunctional's underlying bilinear form.
  /*! If \f$Q(x, y)\f$ is a bilinear form, and \f$a\f$ is a vector, then the
    interior product of \f$a\f$ with \f$Q\f$ is the functional
    \f[
    I(z) = Q(a, z).
    \f]
    Such a functional is always linear and hence can be represented by taking
    the standard dot product with some vector \f$y\f$:
    \f[
    I(z) = y^T z.
    \f]
    This method returns the vector \f$y\f$.
  */
  virtual void interior_product(IMVecType &x, IMVecType &y) {
    this->gradientAt(x, y);
    y.scale(0.5);
  }

  // Assembles the matrix \f$Q_{ij}\f$ corresponding to the bilinear form.
  /* If \f$\{x_i\}\f$ is a basis for the vector space IMVecType,
    \f$Q_{ij}= Q(x_i, x_j)\f$. */
  // virtual void assemble_form(Mat Q) = 0;

};

//! Computes finite difference approximations of a IPFunctional<IceModelVec2S> gradient.
/*! Useful for debugging a hand coded gradient. */
void gradientFD(IPFunctional<IceModelVec2S> &f, IceModelVec2S &x, IceModelVec2S &gradient);

//! Computes finite difference approximations of a IPFunctional<IceModelVec2V> gradient.
/*! Useful for debugging a hand coded gradient. */
void gradientFD(IPFunctional<IceModelVec2V> &f, IceModelVec2V &x, IceModelVec2V &gradient);

} // end of namespace inverse
} // end of namespace pism

#endif /* end of include guard: FUNCTIONAL_HH_1E2DIXE6 */
