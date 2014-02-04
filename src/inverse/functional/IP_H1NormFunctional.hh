// Copyright (C) 2012, 2013, 2014  David Maxwell
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

#ifndef IP_H1NORMFUNCTIONAL_HH_TF8AKRNQ
#define IP_H1NORMFUNCTIONAL_HH_TF8AKRNQ

#include "IPFunctional.hh"


//! Implements a functional corresponding to (the square of) an \f$H^1\f$ norm of a scalar valued function.
/*! The functional is, in continuous terms 
\f[
J(f) = \int_{\Omega} c_{H^1} \left|\nabla f\right|^2 + c_{L^2}f^2 \; dA
\f]
where \f$\Omega\f$ is the square domain. Numerically it is implemented using 
Q1 finite elements.  Integration can be 'restricted', in a sense, to a subset of the domain
using a projection that forces \f$f\f$ to equal zero at nodes specified
by the constructor argument \a dirichletLocations.
*/
class IP_H1NormFunctional2S : public IPInnerProductFunctional<IceModelVec2S> {
public:
  IP_H1NormFunctional2S(IceGrid &grid, ///< computational grid
      double cL2, ///< The constant \f$c_{L^2}\f$.
      double cH1, ///< The constant \f$c_{H^1}\f$.
      IceModelVec2Int *dirichletLocations=NULL ///< Nodes where the function will be set to zero prior to integration.
      ) :
      IPInnerProductFunctional<IceModelVec2S>(grid),
      m_cL2(cL2), m_cH1(cH1), m_dirichletIndices(dirichletLocations) {};
  virtual ~IP_H1NormFunctional2S() {};
  
  virtual PetscErrorCode valueAt(IceModelVec2S &x, double *OUTPUT);
  virtual PetscErrorCode dot(IceModelVec2S &a, IceModelVec2S &b, double *OUTPUT);
  virtual PetscErrorCode gradientAt(IceModelVec2S &x, IceModelVec2S &gradient);
  virtual PetscErrorCode assemble_form(Mat J);

protected:

  double m_cL2, m_cH1;
  IceModelVec2Int *m_dirichletIndices;

private:
  IP_H1NormFunctional2S(IP_H1NormFunctional2S const &);
  IP_H1NormFunctional2S & operator=(IP_H1NormFunctional2S const &);  
};

#endif /* end of include guard: H1NORMFUNCTIONAL_HH_TF8AKRNQ */
