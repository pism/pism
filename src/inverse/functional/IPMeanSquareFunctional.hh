// Copyright (C) 2012, 2013, 2014, 2015, 2020, 2022  David Maxwell
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

#ifndef IPMEANSQUAREFUNCTIONAL_HH_DZ18EO5C
#define IPMEANSQUAREFUNCTIONAL_HH_DZ18EO5C

#include "IPFunctional.hh"

namespace pism {
namespace inverse {

//! Implements a functional corresponding to a (possibly weighted) sum of squares of components of an array::Scalar.
/*! If the vector has components \f$x_i\f$ the functional is
  \f[
  J(x) = c_N \sum_{i} w_i x_i^2
  \f]
  where \f$[w_i]\f$ is a vector of weights and \f$c_N\f$ is a normalization constant. The value
  of the normalization constant is set implicitly by a call to normalize().
*/
class IPMeanSquareFunctional2S : public IPInnerProductFunctional<array::Scalar> {
public:
  /*!
   * @param[in] grid the computational grid
   * @param[in] weights Vector of weights (NULL implies all weights are 1)
   */
  IPMeanSquareFunctional2S(IceGrid::ConstPtr grid, 
                           array::Scalar *weights=NULL)
    : IPInnerProductFunctional<array::Scalar>(grid),
      m_weights(weights),
      m_normalization(1.) {};
  virtual ~IPMeanSquareFunctional2S() {};

  virtual void normalize(double scale);

  virtual void valueAt(array::Scalar &x, double *OUTPUT);
  virtual void dot(array::Scalar &a, array::Scalar &b, double *OUTPUT);
  virtual void gradientAt(array::Scalar &x, array::Scalar &gradient);

protected:
  array::Scalar *m_weights;
  double m_normalization;

private:
  IPMeanSquareFunctional2S(IPMeanSquareFunctional2S const &);
  IPMeanSquareFunctional2S & operator=(IPMeanSquareFunctional2S const &);
};


//! Implements a functional corresponding to a (possibly weighted) sum of squares of components of an array::Scalar.
/*! If the vector has component vectors \f$x_i\f$ the functional is
  \f[
  J(x) = c_N \sum_{i} w_i |x_i|^2
  \f]
  where \f$[w_i]\f$ is a vector of weights and \f$c_N\f$ is a normalization constant. The value
  of the normalization constant is set implicitly by a call to normalize().
*/
class IPMeanSquareFunctional2V : public IPInnerProductFunctional<array::Vector> {
public:
  IPMeanSquareFunctional2V(IceGrid::ConstPtr grid, array::Scalar *weights=NULL) :
    IPInnerProductFunctional<array::Vector>(grid), m_weights(weights), m_normalization(1.) {};
  virtual ~IPMeanSquareFunctional2V() {};

  virtual void normalize(double scale);

  virtual void valueAt(array::Vector &x, double *OUTPUT);
  virtual void dot(array::Vector &a, array::Vector &b, double *OUTPUT);
  virtual void gradientAt(array::Vector &x, array::Vector &gradient);

protected:
  array::Scalar *m_weights;
  double m_normalization;

private:
  IPMeanSquareFunctional2V(IPMeanSquareFunctional2V const &);
  IPMeanSquareFunctional2V & operator=(IPMeanSquareFunctional2V const &);
};


} // end of namespace inverse
} // end of namespace pism

#endif /* end of include guard: IPMEANSQUAREFUNCTIONAL_HH_DZ18EO5C */
