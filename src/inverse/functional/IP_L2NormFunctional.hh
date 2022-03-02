// Copyright (C) 2012, 2014, 2015, 2022  David Maxwell and Constantine Khroulev
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

#ifndef IP_L2NORMFUNCTIONAL_HH_BSVF8BMQ
#define IP_L2NORMFUNCTIONAL_HH_BSVF8BMQ

#include "IPFunctional.hh"

namespace pism {
namespace inverse {

//! Implements a functional corresponding to (the square of) an \f$L^2\f$ norm of a scalar valued function.
/*! The functional is, in continuous terms 
  \f[
  J(f) = \int_{\Omega} f^2 \; dA
  \f]
  where \f$\Omega\f$ is the square domain. Numerically it is implemented using 
  Q1 finite elements.
*/
class IP_L2NormFunctional2S : public IPInnerProductFunctional<array::Scalar> {
public:
  IP_L2NormFunctional2S(IceGrid::ConstPtr grid) : IPInnerProductFunctional<array::Scalar>(grid) {};
  virtual ~IP_L2NormFunctional2S() {};
  
  virtual void valueAt(array::Scalar &x, double *OUTPUT);
  virtual void dot(array::Scalar &a, array::Scalar &b, double *v);
  virtual void gradientAt(array::Scalar &x, array::Scalar &gradient);

private:
  IP_L2NormFunctional2S(IP_L2NormFunctional2S const &);
  IP_L2NormFunctional2S & operator=(IP_L2NormFunctional2S const &);  
};

//! Implements a functional corresponding to (the square of) an \f$L^2\f$ norm of a vector valued function.
/*! The functional is, in continuous terms 
  \f[
  J(f) = \int_{\Omega} f^2 \; dA
  \f]
  where \f$\Omega\f$ is the square domain. Numerically it is implemented using 
  Q1 finite elements.
*/
class IP_L2NormFunctional2V : public IPInnerProductFunctional<array::Vector> {
public:
  IP_L2NormFunctional2V(IceGrid::ConstPtr grid) : IPInnerProductFunctional<array::Vector>(grid) {};
  virtual ~IP_L2NormFunctional2V() {};
  
  virtual void valueAt(array::Vector &x, double *v);
  virtual void dot(array::Vector &a, array::Vector &b, double *v);
  virtual void gradientAt(array::Vector &x, array::Vector &gradient);

private:
  IP_L2NormFunctional2V(IP_L2NormFunctional2V const &);
  IP_L2NormFunctional2V & operator=(IP_L2NormFunctional2V const &);  
};

} // end of namespace inverse
} // end of namespace pism

#endif /* end of include guard: IP_L2NORMFUNCTIONAL_HH_BSVF8BMQ */
