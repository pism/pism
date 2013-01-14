// Copyright (C) 2011, 2012  David Maxwell
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

#ifndef _INVTAUCPARAM_H_
#define _INVTAUCPARAM_H_

#include <petsc.h>
#include "NCVariable.hh"
#include "iceModelVec.hh"

/*! Encapsulates a parameterization of \f$\tau_c\f$.
It is frequently when solving an inverse problem for \f$\tau_c\f$
not to work with \f$\tau_c\f$ directly but to work with a different
variable \f$p\f$, and a relationship \f$\tau_c=F(p)\f$.
A common choice in the glaciology literature is \f$F(p)=p^2\f$,
which ensures that \f$\tau_c\f$ is non-negative, but has the 
disadvantage that it is a 2-1 parameterization.  A potentially
more satisfactory choice is \f$F(p)=e^p\f$, which ensures
positivitiy, is 1-1, and respects the wide scale variations of \f$\tau_c\f$.

An \class InvTaucParameterization encapsulates a parameterization, and
is intended to be used in conjuction with an \class InvSSAForwardProblem 
to indicate the choice of parameterization.

This method of encoding mathematical expressions is flexible and convenient,
but is also slow; it has the overhead that many virtual function calls are
needed if the expression is being called over and over again.  If this
proves to be a significant source of slowness, we could look at 
using the Expression Template idiom, http://drdobbs.com/184401627.
*/
class InvTaucParameterization 
{
public:
  
  InvTaucParameterization(){ /*do nothing*/ };
  
  virtual ~InvTaucParameterization() {};

  virtual PetscErrorCode init( const NCConfigVariable &config);

  /*! Converts from parameterization value \f$p\f$ to \f$\tau_c=F(p)\f$.
  \param p The parameter value.
  \param value On output, the value \f$F(p)\f$.
  \param derivative On output, the value \f$F'(p)\f$. */
  virtual PetscErrorCode toTauc( PetscReal p, PetscReal *value, PetscReal *derivative) = 0;
  
  /*! Converts from \f$\tau_c\f$ to a parameterization value \f$p\f$ such that
   \f$\tau_c=F(p)\f$.  More than one such \f$p\f$ may exist; only one is returned. */
  virtual PetscErrorCode fromTauc( PetscReal tauc, PetscReal *OUTPUT) = 0;

  virtual PetscErrorCode convertToTauc( IceModelVec2S &zeta, IceModelVec2S &tauc, bool communicate = true);

  virtual PetscErrorCode convertFromTauc( IceModelVec2S &tauc, IceModelVec2S &zeta,  bool communicate = true);
protected:
  
  PetscReal m_tauc_scale;
};

/*! Parameterization \f$\tau_c=F(p)\f$ with $F(p)=p$. */
class InvTaucParamIdent: public InvTaucParameterization
{
public:
  InvTaucParamIdent(){ /*do nothing*/ };

  virtual ~InvTaucParamIdent() {};

  virtual PetscErrorCode toTauc( PetscReal p, PetscReal *value, PetscReal *derivative);

  virtual PetscErrorCode fromTauc( PetscReal tauc, PetscReal *OUTPUT);
};

/*! Parameterization \f$\tau_c=F(p)\f$ with $F(p)=p^2$. */
class InvTaucParamSquare: public InvTaucParameterization
{
public:
  InvTaucParamSquare(){ /*do nothing*/ };

  virtual ~InvTaucParamSquare() {};

  virtual PetscErrorCode toTauc( PetscReal p, PetscReal *value, PetscReal *derivative);

  virtual PetscErrorCode fromTauc( PetscReal tauc, PetscReal *OUTPUT);
};

/*! Parameterization \f$\tau_c=F(p)\f$ with $F(p)=exp(p)$. */
class InvTaucParamExp: public InvTaucParameterization
{
public:
  InvTaucParamExp(){ /*do nothing*/ };

  virtual ~InvTaucParamExp() {};

  virtual PetscErrorCode init( const NCConfigVariable &config);

  virtual PetscErrorCode toTauc( PetscReal p, PetscReal *value, PetscReal *derivative);

  virtual PetscErrorCode fromTauc( PetscReal tauc, PetscReal *OUTPUT);

private:
  PetscReal m_tauc_eps;
};


/*! Monotone parameterization \f$\tau_c=F(p)\f$ where \f$F(p)\rightarrow 0\f$
as \f$p\rightarrow-\infty\f$ and \f$F(p)\approx p\f$ for large values of \f$p\f$.  The transition from a nonlinear to an approximately linear 
function occurs in the neighbourhood of the parameter \f$tauc_0\f$.*/
class InvTaucParamTruncatedIdent: public InvTaucParameterization
{
public:
  InvTaucParamTruncatedIdent( ) {};

  virtual ~InvTaucParamTruncatedIdent() {};

  virtual PetscErrorCode init( const NCConfigVariable &config);

  virtual PetscErrorCode toTauc( PetscReal p, PetscReal *value, PetscReal *derivative);

  virtual PetscErrorCode fromTauc( PetscReal tauc, PetscReal *OUTPUT);

private:
  PetscReal m_tauc0_sq;
  PetscReal m_tauc_eps;
};

/*! Some handy pre-instantiated parameterizations. */
extern InvTaucParamIdent g_InvTaucParamIdent;
extern InvTaucParamSquare g_InvTaucParamSquare;
extern InvTaucParamExp g_InvTaucParamExp;

#endif //_INVTAUCPARAM_H_

