// Copyright (C) 2012  David Maxwell
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

#ifndef FUNCTIONAL_HH_1E2DIXE6
#define FUNCTIONAL_HH_1E2DIXE6

#include "iceModelVec.hh"
#include "FETools.hh"



//! Abstract base class for functions from ice model vectors to \f$\mathbb{R}\f$.
/*! Inverse problems frequently involve minimizing a functional,
such such as the misfit
\f[
J(u) = \int_\Omega |u-u_{\rm obs}|^2
\f]
between observed (\f$u_{\rm obs}\f$) and modeled (\f$u\f$)
surface velocities. Subclasses of Functional define such maps,
and permit computation of their gradients.
*/
template<class IMVecType>
class Functional {

public:  
  Functional(IceGrid &grid) : m_grid(grid), m_element_index(m_grid) { 
    m_quadrature.init(m_grid);
  }

  virtual ~Functional() {};

  //! Computes the value of the functional at the vector x.
  virtual PetscErrorCode valueAt(IMVecType &x, PetscReal *OUTPUT) = 0;

  //! Computes the gradient of the functional at the vector x.
  /*! On an \f$m\times n\f$ IceGrid, an IceModelVec \f$x\f$ with \f$d\f$
  degrees of freedom will be \f$d m n\f$-dimensional with components \f$x_i\f$.
  The gradient computed here is the directional derivative \f$\nabla J\f$ of the functional 
  \f$J\f$ with respect to \f$x\f$. Concretely, the \f$i^{\rm th}\f$ component of \f$\nabla J\f$
  is
  \f[
  \nabla J_i = \frac{\partial}{\partial x_i} J(x).
  \f]
  This vector is returned as \c gradient.
  */
  virtual PetscErrorCode gradientAt(IMVecType &x, IMVecType &gradient) = 0;

protected:
  IceGrid &m_grid;

  FEElementMap m_element_index;
  FEQuadrature m_quadrature;
  FEDOFMap     m_dofmap;
  
private:
  // Hide copy/assignment operations
  Functional(Functional const &);
  Functional & operator=(Functional const &);

};

//! Abstract base class for Functionals arising from an inner product.
/*!
Frequently functionals have the structure
\f[
J(u) = Q(u,u)
\f]
where \f$Q\f$ is a symmetric, positive definite, bilinear form. Certain
minimization algorithms only apply to such functionals, which should subclass
from IPFunctional. 
*/
template<class IMVecType>
class IPFunctional : public Functional<IMVecType>{

public:
  IPFunctional(IceGrid &grid) : Functional<IMVecType>(grid) {};

  //! Computes the inner product \f$Q(a,b)\f$.
  virtual PetscErrorCode dot(IMVecType &a, IMVecType &b, PetscReal *OUTPUT) = 0;

  //! Computes the interior product of a vector with the IPFunctional's underlying bilinear form.
  /*! If \f$Q(x,y)\f$ is a bilinear form, and \f$a\f$ is a vector, then the 
  interior product of \f$a\f$ with \f$Q\f$ is the functional 
  \f[
  I(z) = Q(a,z).
  \f]
  Such a functional is always linear and hence can be represented by taking
  the standard dot product with some vector \f$y\f$:
  \f[
  I(z) = y^T z.
  \f]
  This method returns the vector \f$y\f$.
  */  
  virtual PetscErrorCode interior_product(IMVecType &x, IMVecType &y) {
    PetscErrorCode ierr;
    ierr = this->gradientAt(x,y); CHKERRQ(ierr);
    ierr = y.scale(0.5); CHKERRQ(ierr);
    return 0;
  }
};

//! Computes finite difference approximations of a Functional<IceModelVec2S> gradient.
/*! Useful for debugging a hand coded gradient. */
PetscErrorCode gradientFD(Functional<IceModelVec2S> &f, IceModelVec2S &x, IceModelVec2S &gradient);

//! Computes finite difference approximations of a Functional<IceModelVec2V> gradient.
/*! Useful for debugging a hand coded gradient. */
PetscErrorCode gradientFD(Functional<IceModelVec2V> &f, IceModelVec2V &x, IceModelVec2V &gradient);

#endif /* end of include guard: FUNCTIONAL_HH_1E2DIXE6 */
