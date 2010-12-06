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
#include "PISMDiagnostic.hh"

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
class SSAStrengthExtension {
public:
  SSAStrengthExtension() {
    // FIXME: use the config database.

    min_thickness = 50.0;   // m
    // minimum thickness (for SSA velocity computation) at which 
    // NuH switches from vertical integral to constant value
    // this value strongly related to calving front
    // force balance, but the geometry itself is not affected by this value
    const PetscReal
      DEFAULT_CONSTANT_HARDNESS_FOR_SSA = 1.9e8,  // Pa s^{1/3}; see p. 49 of MacAyeal et al 1996
      DEFAULT_TYPICAL_STRAIN_RATE = (100.0 / secpera) / (100.0 * 1.0e3);  // typical strain rate is 100 m/yr per 
    nuH = min_thickness * DEFAULT_CONSTANT_HARDNESS_FOR_SSA
      / (2.0 * pow(DEFAULT_TYPICAL_STRAIN_RATE,2./3.)); // Pa s m
    // COMPARE: 30.0 * 1e6 * secpera = 9.45e14 is Ritz et al (2001) value of
    //          30 MPa yr for \bar\nu
  }
  virtual ~SSAStrengthExtension() {}
  //! Set strength with units (viscosity times thickness).
  virtual PetscErrorCode set_notional_strength(PetscReal my_nuH)
  { nuH = my_nuH; return 0; }
  //! Set minimum thickness to trigger use of extension.
  virtual PetscErrorCode set_min_thickness(PetscReal my_min_thickness)
  { min_thickness = my_min_thickness; return 0; }
  //! Returns strength with units (viscosity times thickness).
  virtual PetscReal      get_notional_strength() const { return nuH; }
  //! Returns minimum thickness to trigger use of extension.
  virtual PetscReal      get_min_thickness() const { return min_thickness; }
private:
  PetscReal  min_thickness, nuH;
};


//! PISM's SSA solver implementation
class SSAFD : public ShallowStressBalance
{
  friend class SSAFD_taud;
public:
  SSAFD(IceGrid &g, IceBasalResistancePlasticLaw &b, IceFlowLaw &i, EnthalpyConverter &e,
        const NCConfigVariable &c);

  virtual ~SSAFD() { deallocate(); }

  SSAStrengthExtension strength_extension;

  virtual PetscErrorCode init(PISMVars &vars);

  virtual PetscErrorCode update(bool fast);

  virtual PetscErrorCode set_initial_guess(IceModelVec2V &guess);

  virtual PetscErrorCode write_model_state(string filename);

  //! Add pointers to diagnostic quantities to a dictionary.
  virtual void get_diagnostics(map<string, PISMDiagnostic*> &dict);

protected:
  virtual PetscErrorCode allocate();

  virtual PetscErrorCode deallocate();

  virtual PetscErrorCode solve();

  virtual PetscErrorCode compute_nuH_staggered(IceModelVec2Stag &result,
                                               PetscReal epsilon);

  virtual PetscErrorCode compute_nuH_norm(PetscReal &norm,
                                          PetscReal &norm_change);

  virtual PetscErrorCode assemble_matrix(bool include_basal_shear, Mat A);

  virtual PetscErrorCode assemble_rhs(Vec rhs);

  virtual PetscErrorCode compute_driving_stress(IceModelVec2V &taud);

  virtual PetscErrorCode compute_hardav_staggered(IceModelVec2Stag &result);

  virtual PetscErrorCode compute_basal_frictional_heating(IceModelVec2S &result);

  virtual PetscErrorCode compute_D2(IceModelVec2S &result);

  virtual PetscErrorCode compute_maximum_velocity();

  virtual PetscErrorCode writeSSAsystemMatlab();

  virtual PetscErrorCode update_nuH_viewers();

  virtual PetscErrorCode stdout_report(string &result);

  IceModelVec2Mask *mask;
  IceModelVec2S *thickness, *tauc, *surface, *bed;
  IceModelVec2Stag hardness, nuH, nuH_old;
  IceModelVec2V taud, velocity_old;
  IceModelVec3 *enthalpy;

  string stdout_ssa;

  // objects used by the SSA solver (internally)
  KSP SSAKSP;
  Mat SSAStiffnessMatrix;
  Vec SSAX, SSARHS;  // Global vectors for solution of the linear system and the RHS.
  DA  SSADA;
};

//! \brief Computes the driving stress (taud).
class SSAFD_taud : public PISMDiag<SSAFD>
{
public:
  SSAFD_taud(SSAFD *m, IceGrid &g, PISMVars &my_vars);
  PetscErrorCode compute(IceModelVec* &result);
};

#endif /* _SSAFD_H_ */
