// Copyright (C) 2006-2010 Ed Bueler, Constantine Khroulev, and Jed Brown
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

#ifndef __stressBalanceSSA_hh
#define __stressBalanceSSA_hh

#include <petscksp.h>
#include "pism_const.hh"
#include "grid.hh"
#include "materials.hh"
#include "iceModelVec.hh"


//! Where ice thickness is zero the SSA is no longer "elliptic".  This class provides an extension coefficient to maintain well-posedness/ellipticity.
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
class newSSAStrengthExtension {
public:
  newSSAStrengthExtension();
  //! Set strength with units (viscosity times thickness).
  PetscErrorCode set_notional_strength(PetscReal my_nuH);
  //! Set minimum thickness to trigger use of extension.
  PetscErrorCode set_min_thickness(PetscReal my_min_thickness);
  //! Returns strength with units (viscosity times thickness).
  PetscReal      get_notional_strength() const { return nuH; }
  //! Returns minimum thickness to trigger use of extension.
  PetscReal      get_min_thickness_for_extension() const { return min_thickness; }
private:
  PetscReal  min_thickness, nuH;
};


//! Solve shallow shelf approximation (SSA) ice sheet stress balance for PISM.
class StressBalanceSSA {

public:
  StressBalanceSSA(IceGrid* g, IceFlowLaw* ssa_ice, IceBasalResistancePlasticLaw* ssa_basal,
                   IceModelVec2S *ssa_tauc, IceModelVec2S *ssa_mask, IceModelVec2S *ssa_hardav);
  virtual ~StressBalanceSSA();

  //! Sets initial-guess velocities to zero.
  virtual PetscErrorCode setGuessZero();
  //! Provide initial-guess velocities.
  virtual PetscErrorCode setGuess(IceModelVec2S *ubar_guess, IceModelVec2S *vbar_guess);

  //! Solve stress balance for vertically-integrated horizontal velocity.
  virtual PetscErrorCode solve(PetscInt *numiter);
  //! Solve stress balance for vertically-integrated horizontal velocity, but provide the effective viscosity (two components on the staggered grid).
  virtual PetscErrorCode solve(IceModelVec2S vNuH[2], PetscInt *numiter);

  //! Put components of basal stress, applied by bedrock to base of ice, in provided IceModelVec2.
  virtual PetscErrorCode getBasalStress(IceModelVec2S *vbs_x, IceModelVec2S *vbs_y);

  //! Modify/update the provided field (e.g. by adding SSA solution to current SIA solution).
  virtual PetscErrorCode modifyHorizontalVelocityComponents(IceModelVec3 *u, IceModelVec3 *v);

  //! Modify/update the provided field of strain heating values.
  virtual PetscErrorCode modifySigma(IceModelVec3 *sigma);

  //! Modify/update the provided field of basal friction rates.
  virtual PetscErrorCode modifyBasalFrictionalHeating();

protected:
  virtual PetscErrorCode updateAveragedHardness();
  virtual PetscErrorCode computeEffectiveViscosity(IceModelVec2S vNuH[2], PetscReal epsilon);
  virtual PetscErrorCode testConvergenceOfNu(IceModelVec2S vNuH[2], IceModelVec2S vNuHOld[2],
                                             PetscReal *norm, PetscReal *normChange);
  virtual PetscErrorCode assembleSSAMatrix(bool includeBasalShear, IceModelVec2S vNuH[2], Mat A);
  virtual PetscErrorCode assembleSSARhs(bool surfGradInward, Vec rhs);
  virtual PetscErrorCode moveVelocityToDAVectors(Vec x);

  IceModelVec2V ssavel;
  IceModelVec2S vaveragedhardness;
  
  KSP        SSAKSP;
  Mat        SSAStiffnessMatrix;
  Vec        SSAX,      //!< Global vector for solution of the linear system  Ax=b.
             SSARHS;    //!< Global vector for right-hand side of the linear system  Ax=b.
  Vec        SSAXLocal; //!< Local copy of the solution uses this to map back to a DA-based vector.
  VecScatter SSAScatterGlobalToLocal;

  IceGrid*         grid;
  NCConfigVariable config;

  IceFlowLaw*                   ice;
  IceBasalResistancePlasticLaw* basal;
  SSAStrengthExtension          ssaStrengthExtend;

  IceModelVec2S *vtauc, *vmask, *vhardav;  // pointers to external IMVecs (e.g. from IceModel)

private:
  PetscErrorCode initAndAllocate(IceGrid* g);
  PetscErrorCode deallocate();

  PetscScalar basalDragx(PetscScalar **tauc, PISMVector2 **uv,
                         PetscInt i, PetscInt j) const;
  PetscScalar basalDragy(PetscScalar **tauc, PISMVector2 **uv,
                         PetscInt i, PetscInt j) const;
};

#endif

