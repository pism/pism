// Copyright (C) 2011, 2012, 2013  David Maxwell
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

#ifndef _INVTAUCPARAM_H_
#define _INVTAUCPARAM_H_

#include <petsc.h>
#include "NCVariable.hh"
#include "iceModelVec.hh"

//! Encapsulates a parameterization of \f$\tau_c\f$ as a function of a parameter \f$\zeta\f$.
/*!
When solving an inverse problem for \f$\tau_c\f$, one frequently does
not work with \f$\tau_c\f$ directly but to work with a different
variable \f$\zeta\f$, and a relationship \f$\tau_c=g(\zeta)\f$.
A common choice in the glaciology literature is \f$g(\zeta)=\zeta^2\f$,
which ensures that \f$\tau_c\f$ is non-negative, but has the 
disadvantage that it is a 2-1 parameterization.  A potentially
more satisfactory choice is \f$g(\zeta)=e^\zeta\f$, which ensures
positivitiy, is 1-1, and respects the wide scale variations of \f$\tau_c\f$.

An IPTaucParameterization encapsulates a parameterization, and
is intended to be used in conjuction with an \ref IP_SSATaucForwardProblem 
to indicate the choice of parameterization.

This method of encoding mathematical expressions is flexible and convenient,
but is also slow; it has the overhead that many virtual function calls are
needed if the expression is being called over and over again.  If this
proves to be a significant source of slowness, we could look at 
using the Expression Template idiom, http://drdobbs.com/184401627.

Internally, the functions use units of \f$\tau_c\f$ such that 
the config parameter \a tauc_param_tauc_scale equals one. I.e.
for a conversion function \f$g(\zeta)=\zeta^2\f$, then
\f[
\frac{\tau_c} = \tau_{\rm scale}g(\zeta^2).
\f]
where \f$\tau_{\rm scale} = \text{\tt tauc\_param\_tauc\_scale}\f$.
*/
class IPTaucParameterization 
{
public:
  
  IPTaucParameterization(){ /*do nothing*/ };
  
  virtual ~IPTaucParameterization() {};

  virtual PetscErrorCode init( const NCConfigVariable &config);

  //! Converts from parameterization value \f$\zeta\f$ to \f$\tau_c=g(\zeta)\f$.
  /*!
  \param[in] zeta The parameter value.
  \param[out] value The value \f$g(p)\f$.
  \param[out] derivative The value \f$f'(p)\f$. */
  virtual PetscErrorCode toTauc( PetscReal zeta, PetscReal *value, PetscReal *derivative) = 0;
  
  //! Converts from \f$\tau_c\f$ to a parameterization value \f$\zeta\f$ such that \f$\tau_c=g(\zeta)\f$.  
  /*! More than one such \f$\zeta\f$ may exist; only one is returned. */
  virtual PetscErrorCode fromTauc( PetscReal tauc, PetscReal *OUTPUT) = 0;

  //! Transforms a vector of \f$\zeta\f$ values to a vector of \f$\tau_c\f$ values.
  virtual PetscErrorCode convertToTauc( IceModelVec2S &zeta, IceModelVec2S &tauc, bool communicate = true);

  //! Transforms a vector of \f$\tau_c\f$ values to a vector of \f$\zeta\f$ values.
  virtual PetscErrorCode convertFromTauc( IceModelVec2S &tauc, IceModelVec2S &zeta,  bool communicate = true);
protected:
  
  PetscReal m_tauc_scale;  ///< Value of \f$\tau_c\f$ in PISM units that equals 1 for IPTaucParameterization's.
};

//! Parameterization \f$\tau_c=\tau_{\rm scale}g(\zeta)\f$ with \f$g(\zeta)=\zeta\f$.
class IPTaucParamIdent: public IPTaucParameterization
{
public:
  IPTaucParamIdent(){ /*do nothing*/ };

  virtual ~IPTaucParamIdent() {};

  virtual PetscErrorCode toTauc( PetscReal p, PetscReal *value, PetscReal *derivative);

  virtual PetscErrorCode fromTauc( PetscReal tauc, PetscReal *OUTPUT);
};

//! Parameterization \f$\tau_c=\tau_{\rm scale}g(\zeta)\f$ with \f$g(\zeta)=\zeta^2\f$.
class IPTaucParamSquare: public IPTaucParameterization
{
public:
  IPTaucParamSquare(){ /*do nothing*/ };

  virtual ~IPTaucParamSquare() {};

  virtual PetscErrorCode toTauc( PetscReal p, PetscReal *value, PetscReal *derivative);

  virtual PetscErrorCode fromTauc( PetscReal tauc, PetscReal *OUTPUT);
};

//! Parameterization \f$\tau_c=\tau_{\rm scale}g(\zeta)\f$ with \f$g(\zeta)=\exp(\zeta)\f$.
class IPTaucParamExp: public IPTaucParameterization
{
public:
  IPTaucParamExp(){ /*do nothing*/ };

  virtual ~IPTaucParamExp() {};

  virtual PetscErrorCode init( const NCConfigVariable &config);

  virtual PetscErrorCode toTauc( PetscReal p, PetscReal *value, PetscReal *derivative);

  virtual PetscErrorCode fromTauc( PetscReal tauc, PetscReal *OUTPUT);

private:
  PetscReal m_tauc_eps;
};


//! A monotone non-negative parameterization \f$\tau_c=\tau_{\rm scale}g(\zeta)\f$ that is approximately the identity away from small values of  \f$\tau_c\f$
/*! More specifically, \f$g(\zeta)\rightarrow 0\f$ as \f$\zeta\rightarrow-\infty\f$ and \f$g(\zeta)\approx p\f$ 
    for large values of \f$\zeta\f$.  The transition from a nonlinear to an approximately linear 
    function occurs in the neighbourhood of the parameter \f$tauc_0\f$. */
class IPTaucParamTruncatedIdent: public IPTaucParameterization
{
public:
  IPTaucParamTruncatedIdent( ) {};

  virtual ~IPTaucParamTruncatedIdent() {};

  virtual PetscErrorCode init( const NCConfigVariable &config);

  virtual PetscErrorCode toTauc( PetscReal p, PetscReal *value, PetscReal *derivative);

  virtual PetscErrorCode fromTauc( PetscReal tauc, PetscReal *OUTPUT);

private:
  PetscReal m_tauc0_sq;
  PetscReal m_tauc_eps;
};

#endif //_INVTAUCPARAM_H_

