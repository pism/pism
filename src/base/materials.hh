// Copyright (C) 2004-2010 Jed Brown, Ed Bueler, and Constantine Khroulev
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

#ifndef __materials_hh
#define __materials_hh

#include <petsc.h>
#include "enthalpyConverter.hh"
#include "NCVariable.hh"

//! Class containing physical constants and the constitutive relation describing till for SSA.
/*!
This \e pseudo -plastic type can actually describe anything from linearly 
viscous till to purely plastic till.
 */
class IceBasalResistancePlasticLaw {
public:
  IceBasalResistancePlasticLaw(PetscScalar regularizationConstant, bool pseudoPlastic,
                   PetscScalar pseudoExponent, PetscScalar pseudoUThreshold);
  virtual PetscErrorCode printInfo(int verbthresh, MPI_Comm com);
  virtual PetscScalar drag(PetscScalar tauc,
                           PetscScalar vx, PetscScalar vy);
  // Also get the derivative of drag with respect to \f$ alpha=\frac 1 2 \abs{u}^2 \f$.
  virtual void dragWithDerivative(PetscReal tauc, PetscScalar vx, PetscScalar vy,
                                  PetscScalar *drag, PetscScalar *ddrag) const;
  virtual ~IceBasalResistancePlasticLaw() {} // class w virtual methods needs virtual destructor?

  PetscReal   plastic_regularize, pseudo_q, pseudo_u_threshold;
  bool pseudo_plastic;
};


//! Extension coefficient to maintain well-posedness/"ellipticity" of SSA where ice thickness is zero.
/*!
More specifically, the SSA equations are
\latexonly
\def\ddt#1{\ensuremath{\frac{\partial #1}{\partial t}}}
\def\ddx#1{\ensuremath{\frac{\partial #1}{\partial x}}}
\def\ddy#1{\ensuremath{\frac{\partial #1}{\partial y}}}
\begin{equation*}
  - 2 \ddx{}\left[\nu H \left(2 \ddx{u} + \ddy{v}\right)\right]
        - \ddy{}\left[\nu H \left(\ddy{u} + \ddx{v}\right)\right]
        + \tau_{(b)x}  =  - \rho g H \ddx{h},
\end{equation*}
\endlatexonly
and another similar equation.  Schoof \ref SchoofStream shows that these PDEs
are the variational equations for a functional.

The quantity \f$\nu H\f$ is the nonlinear coefficient in this (roughly-speaking)
elliptic pair of PDEs.  Conceptually it is a membrane strength.  Well-posedness of the SSA problem 
requires either a precisely-defined boundary and an appropriate boundary condition
\e or a nonzero value of \f$\nu H\f$ at all points.  This class provides that nonzero value.
 */
class SSAStrengthExtension {
public:
  SSAStrengthExtension();
  virtual ~SSAStrengthExtension();
  //! Set strength with units (viscosity times thickness).
  virtual PetscErrorCode set_notional_strength(PetscReal my_nuH);
  //! Set minimum thickness to trigger use of extension.
  virtual PetscErrorCode set_min_thickness(PetscReal my_min_thickness);
  virtual PetscReal notional_strength() const;           //!< Returns strength with units (viscosity times thickness).
  virtual PetscReal min_thickness_for_extension() const; //!< Returns minimum thickness to trigger use of extension.
private:
  PetscReal  min_thickness, nuH;
};


#endif /* __materials_hh */

