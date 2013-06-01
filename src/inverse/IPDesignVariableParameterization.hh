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


#ifndef IPTAUCPARAMETERIZATION_HH_7ZRETI1S
#define IPTAUCPARAMETERIZATION_HH_7ZRETI1S


#include <petsc.h>
#include "NCVariable.hh"
#include "iceModelVec.hh"

//! Encapsulates a parameterization of a design variable (e.g. \f$\tau_c\f$ for SSA inversions)
//! as a function of a different parameter \f$\zeta\f$.
/*!
When solving an inverse problem for a design variable \f$d\f$ (think of
\f$\tau_c\f$ or hardness for SSA inversions), one frequently does
not work with \f$\d\f$ directly but with a different
variable \f$\zeta\f$, and a relationship \f$d=g(\zeta)\f$.
A common choice in the glaciology literature for \f$\tau_c\f$
is \f$\tau_c=g(\zeta)=\zeta^2\f$, which ensures that \f$\tau_c\f$ is 
non-negative, but has the disadvantage that it is a 2-1 parameterization.  A potentially
more satisfactory choice is \f$g(\zeta)=e^\zeta\f$, which ensures
positivitiy, is 1-1, and respects the wide scale variations of \f$\tau_c\f$.

An IPDesignVariableParameterization implements such a parameterization map.

This method of encoding mathematical expressions is flexible and convenient,
but is also slow; it has the overhead that many virtual function calls are
needed if the expression is being called over and over again.  If this
proves to be a significant source of slowness, we could look at 
using the Expression Template idiom, http://drdobbs.com/184401627.

For certain Tikhonov inversions, it is important to mainain well-scaled
variables.  If the design parameter name is 'foo', internally the parameterizations 
use units of \f$d\f$ such that  the config parameter \a design_param_foo_scale equals one. I.e.
if 'foo' is 'tauc', then  for a conversion function \f$g(\zeta)=\zeta^2\f$, 
\f[
\frac{d} = d_{\rm scale}g(\zeta^2).
\f]
where \f$d_{\rm scale}={\tt design_param_tauc_scale}\f$.
*/
class IPDesignVariableParameterization 
{
public:
  
  IPDesignVariableParameterization() {};
  
  virtual ~IPDesignVariableParameterization() {};

  //! Initializes the scale parameters of the parameterization.
  /*! Every IPDesignVariableParameterization has an associated scale for the design variable 
  \f$d_{\rm scale}\f that equals 1 in internal units.  The scale for a design variable named \a foo
  is stored in an NCConfigVariable file as design_param_foo_scale.  Subclasses may have additional
  parameters that are follow the naming convention \a design_param_foo_*.
  
  \param config          The config file to read the scale parameters from.
  \param design_var_name The associated name of the design variable, e.g. 'tauc' or 'hardav'
  */  
  virtual PetscErrorCode set_scales( const NCConfigVariable &config, const char *design_var_name);

  //! Converts from parameterization value \f$\zeta\f$ to \f$d=g(\zeta)\f$.
  /*!
  \param[in] zeta The parameter value.
  \param[out] value The value \f$g(\zeta)\f$.
  \param[out] derivative The value \f$g'(\zeta)\f$. */
  virtual PetscErrorCode toDesignVariable( PetscReal zeta, PetscReal *value, PetscReal *derivative) = 0;
  
  //! Converts from \f$d\f$ to a parameterization value \f$\zeta\f$ such that \f$d=g(\zeta)\f$.  
  /*! More than one such \f$\zeta\f$ may exist; only one is returned. */
  virtual PetscErrorCode fromDesignVariable( PetscReal d, PetscReal *OUTPUT) = 0;

  //! Transforms a vector of \f$\zeta\f$ values to a vector of \f$d\f$ values.
  virtual PetscErrorCode convertToDesignVariable( IceModelVec2S &zeta, IceModelVec2S &d, bool communicate = true);

  //! Transforms a vector of \f$d\f$ values to a vector of \f$\zeta\f$ values.
  virtual PetscErrorCode convertFromDesignVariable( IceModelVec2S &d, IceModelVec2S &zeta,  bool communicate = true);
protected:
  
  PetscReal m_d_scale;  ///< Value of \f$\d\f$ in PISM units that equals 1 for IPDesignVariableParameterization's units.
};

//! Parameterization \f$d=d_{\rm scale}g(\zeta)\f$ with \f$g(\zeta)=\zeta\f$.
class IPDesignVariableParamIdent: public IPDesignVariableParameterization
{
public:
  IPDesignVariableParamIdent(){ /*do nothing*/ };

  virtual ~IPDesignVariableParamIdent() {};

  virtual PetscErrorCode toDesignVariable( PetscReal p, PetscReal *value, PetscReal *derivative);

  virtual PetscErrorCode fromDesignVariable( PetscReal tauc, PetscReal *OUTPUT);
};

//! Parameterization \f$\tau_c=\tau_{\rm scale}g(\zeta)\f$ with \f$g(\zeta)=\zeta^2\f$.
class IPDesignVariableParamSquare: public IPDesignVariableParameterization
{
public:
  IPDesignVariableParamSquare(){ /*do nothing*/ };

  virtual ~IPDesignVariableParamSquare() {};

  virtual PetscErrorCode toDesignVariable( PetscReal p, PetscReal *value, PetscReal *derivative);

  virtual PetscErrorCode fromDesignVariable( PetscReal tauc, PetscReal *OUTPUT);
};

//! Parameterization \f$\tau_c=\tau_{\rm scale}g(\zeta)\f$ with \f$g(\zeta)=\exp(\zeta)\f$.
class IPDesignVariableParamExp: public IPDesignVariableParameterization
{
public:
  IPDesignVariableParamExp(){ /*do nothing*/ };

  virtual ~IPDesignVariableParamExp() {};

  virtual PetscErrorCode set_scales( const NCConfigVariable &config, const char *design_var_name);

  virtual PetscErrorCode toDesignVariable( PetscReal p, PetscReal *value, PetscReal *derivative);

  virtual PetscErrorCode fromDesignVariable( PetscReal tauc, PetscReal *OUTPUT);

private:
  PetscReal m_d_eps;
};


//! A monotone non-negative parameterization \f$\tau_c=\tau_{\rm scale}g(\zeta)\f$ that is approximately the identity away from small values of  \f$\tau_c\f$
/*! More specifically, \f$g(\zeta)\rightarrow 0\f$ as \f$\zeta\rightarrow-\infty\f$ and \f$g(\zeta)\approx p\f$ 
    for large values of \f$\zeta\f$.  The transition from a nonlinear to an approximately linear 
    function occurs in the neighbourhood of the parameter \f$d_0\f$. */
class IPDesignVariableParamTruncatedIdent: public IPDesignVariableParameterization
{
public:
  IPDesignVariableParamTruncatedIdent( ) {};

  virtual ~IPDesignVariableParamTruncatedIdent() {};

  virtual PetscErrorCode set_scales( const NCConfigVariable &config, const char *design_var_name);

  virtual PetscErrorCode toDesignVariable( PetscReal p, PetscReal *value, PetscReal *derivative);

  virtual PetscErrorCode fromDesignVariable( PetscReal d, PetscReal *OUTPUT);

private:
  PetscReal m_d0_sq;
  PetscReal m_d_eps;
};

#endif /* end of include guard: IPTAUCPARAMETERIZATION_HH_7ZRETI1S */
