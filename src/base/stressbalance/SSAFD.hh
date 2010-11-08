// Copyright (C) 2004--2010 Jed Brown, Ed Bueler and Constantine Khroulev
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

#ifndef _SSAFD_H_
#define _SSAFD_H_

#include "ShallowStressBalance.hh"
#include <petscksp.h>

// //! Where ice thickness is zero the SSA is no longer "elliptic".  This class provides an extension coefficient to maintain well-posedness/ellipticity.
// /*!
// More specifically, the SSA equations are
// \latexonly
// \def\ddt#1{\ensuremath{\frac{\partial #1}{\partial t}}}
// \def\ddx#1{\ensuremath{\frac{\partial #1}{\partial x}}}
// \def\ddy#1{\ensuremath{\frac{\partial #1}{\partial y}}}
// \begin{equation*}
//   - 2 \ddx{}\left[\nu H \left(2 \ddx{u} + \ddy{v}\right)\right]
//         - \ddy{}\left[\nu H \left(\ddy{u} + \ddx{v}\right)\right]
//         + \tau_{(b)x}  =  - \rho g H \ddx{h},
// \end{equation*}
// \endlatexonly
// and another similar equation.  Schoof \ref SchoofStream shows that these PDEs
// are the variational equations for a functional.

// The quantity \f$\nu H\f$ is the nonlinear coefficient in this (roughly-speaking)
// elliptic pair of PDEs.  Conceptually it is a membrane strength.  Well-posedness of the SSA problem 
// requires either a precisely-defined boundary and an appropriate boundary condition
// \e or a nonzero value of \f$\nu H\f$ at all points.  This class provides that nonzero value.
//  */
// class SSAStrengthExtension {
// public:
//   SSAStrengthExtension();
//   //! Set strength with units (viscosity times thickness).
//   PetscErrorCode set_notional_strength(PetscReal my_nuH);
//   //! Set minimum thickness to trigger use of extension.
//   PetscErrorCode set_min_thickness(PetscReal my_min_thickness);
//   //! Returns strength with units (viscosity times thickness).
//   PetscReal      get_notional_strength() const { return nuH; }
//   //! Returns minimum thickness to trigger use of extension.
//   PetscReal      get_min_thickness() const { return min_thickness; }
// private:
//   PetscReal  min_thickness, nuH;
// };


//! PISM's SSA solver implementation
class SSAFD : public ShallowStressBalance
{
public:
  SSAFD(IceGrid &g, IceBasalResistancePlasticLaw &b, IceFlowLaw &i, EnthalpyConverter &e,
        const NCConfigVariable &c)
    : ShallowStressBalance(g, b, i, e, c) {}

  virtual ~SSAFD() { deallocate(); }

  SSAStrengthExtension strength_extension;
  virtual PetscErrorCode init(PISMVars &vars);
  virtual PetscErrorCode update(bool fast);

protected:
  virtual PetscErrorCode allocate_internals();

  virtual PetscErrorCode deallocate();

  virtual PetscErrorCode solve();

  virtual PetscErrorCode compute_nuH_staggered(IceModelVec2Stag &result,
                                               PetscReal epsilon); // done

  virtual PetscErrorCode compute_nuH_norm(PetscReal &norm,
                                          PetscReal &norm_change); // done

  virtual PetscErrorCode assemble_matrix(bool include_basal_shear, Mat A);

  virtual PetscErrorCode assemble_rhs(Vec rhs); // done

  virtual PetscErrorCode compute_driving_stress(IceModelVec2V &taud);

  virtual PetscErrorCode compute_hardav_staggered(IceModelVec2Stag &result); // done

  virtual PetscErrorCode compute_basal_frictional_heating(IceModelVec2S &result); // done

  virtual PetscErrorCode compute_D2(IceModelVec2S &result); // done

  IceModelVec2Mask *mask;
  IceModelVec2S *thickness, *tauc, *surface, *bed;
  IceModelVec2Stag hardness, nuH, nuH_old;
  IceModelVec2V taud;
  IceModelVec3 *enthalpy;

  // objects used by the SSA solver (internally)
  KSP SSAKSP;
  Mat SSAStiffnessMatrix;
  Vec SSAX, SSARHS;  // Global vectors for solution of the linear system and the RHS.
  DA  SSADA;
};


#endif /* _SSAFD_H_ */
