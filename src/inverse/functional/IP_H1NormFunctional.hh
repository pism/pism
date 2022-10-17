// Copyright (C) 2012, 2013, 2014, 2015, 2020  David Maxwell and Constantine Khroulev
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

namespace pism {
namespace inverse {


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
class IP_H1NormFunctional2S : public IPInnerProductFunctional<array::Scalar> {
public:
  /*!
   * @param[in] grid computational grid
   * @param[in] cL2 The constant \f$c_{L^2}\f$.
   * @param[in] cH1 The constant \f$c_{H^1}\f$.
   * @param[in] dirichletLocations Nodes where the function will be set to zero prior to integration.
   */
  IP_H1NormFunctional2S(IceGrid::ConstPtr grid,
                        double cL2,
                        double cH1,
                        array::Scalar *dirichletLocations=NULL)
    : IPInnerProductFunctional<array::Scalar>(grid),
      m_cL2(cL2),
      m_cH1(cH1),
      m_dirichletIndices(dirichletLocations) {};
  virtual ~IP_H1NormFunctional2S() {};
  
  virtual void valueAt(array::Scalar &x, double *OUTPUT);
  virtual void dot(array::Scalar &a, array::Scalar &b, double *OUTPUT);
  virtual void gradientAt(array::Scalar &x, array::Scalar &gradient);
  virtual void assemble_form(Mat J);

protected:

  double m_cL2, m_cH1;
  array::Scalar *m_dirichletIndices;

private:
  IP_H1NormFunctional2S(IP_H1NormFunctional2S const &);
  IP_H1NormFunctional2S & operator=(IP_H1NormFunctional2S const &);  
};

} // end of namespace inverse
} // end of namespace pism

#endif /* end of include guard: H1NORMFUNCTIONAL_HH_TF8AKRNQ */
