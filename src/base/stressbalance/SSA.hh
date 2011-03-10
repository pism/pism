// Copyright (C) 2004--2011 Jed Brown, Ed Bueler and Constantine Khroulev
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

#ifndef _SSA_H_
#define _SSA_H_

#include "ShallowStressBalance.hh"
#include "PISMDiagnostic.hh"

//! Gives an extension coefficient to maintain ellipticity of SSA where ice is thin.
/*!
The SSA is basically a nonlinear elliptic, but vector-valued, equation which
determines the ice velocity field from the driving stress, the basal shear
stress, the ice hardness, and some boundary conditions.  The problem loses
ellipticity (coercivity) if the thickness actually goes to zero.  This class
provides an extension coefficient to maintain ellipticity.

More specifically, the SSA equations are
\f[
\def\ddx#1{\ensuremath{\frac{\partial #1}{\partial x}}}
\def\ddy#1{\ensuremath{\frac{\partial #1}{\partial y}}}
  - 2 \ddx{}\left[\nu H \left(2 \ddx{u} + \ddy{v}\right)\right]
        - \ddy{}\left[\nu H \left(\ddy{u} + \ddx{v}\right)\right]
        + \tau_{(b)x}  =  - \rho g H \ddx{h},
\f]
and another similar equation for the \f$y\f$-component.  Schoof
\ref SchoofStream shows that these PDEs are the variational equations for a
coercive functional, thus (morally) elliptic.

The quantity \f$\nu H\f$ is the nonlinear coefficient, and conceptually it is a
membrane strength.  This class extends \f$\nu H\f$ to have a minimum value
at all points.  It is a class, and not just a configuration constant, because 
setting both the thickness \f$H\f$ and the value \f$\nu H\f$ are allowed, and
setting each of these does not affect the value of the other.
 */
class SSAStrengthExtension {
public:
  SSAStrengthExtension(const NCConfigVariable &c) : config(c) {
    min_thickness = config.get("min_thickness_strength_extension_ssa");
    constant_nu = config.get("constant_nu_strength_extension_ssa");
  }

  virtual ~SSAStrengthExtension() {}

  //! Set strength = (viscosity times thickness).
  /*! Determines nu by input strength and current min_thickness. */
  virtual PetscErrorCode set_notional_strength(PetscReal my_nuH) {
     if (my_nuH <= 0.0) SETERRQ(1,"nuH must be positive");
     constant_nu = my_nuH / min_thickness;
     return 0;
  }

  //! Set minimum thickness to trigger use of extension.
  /*! Preserves strength (nuH) by also updating using current nu.  */
  virtual PetscErrorCode set_min_thickness(PetscReal my_min_thickness) {
     if (my_min_thickness <= 0.0) SETERRQ(1,"min_thickness must be positive");
     PetscReal nuH = constant_nu * min_thickness;
     min_thickness = my_min_thickness;
     constant_nu = nuH / min_thickness;
     return 0;
  }

  //! Returns strength = (viscosity times thickness).
  virtual PetscReal get_notional_strength() const {
    return constant_nu * min_thickness;
  }

  //! Returns minimum thickness to trigger use of extension.
  virtual PetscReal get_min_thickness() const { return min_thickness; }

private:
  const NCConfigVariable &config;
  PetscReal  min_thickness, constant_nu;
};

//! Callback for constructing a new SSA subclass.  The caller is
//! responsible for deleting the newly constructed SSA.
/*! The factory idiom gives a way to implement runtime polymorphism for the 
choice of SSA algorithm.  The factory is a function pointer that takes 
all the arguments of an SSA constructor and returns a newly constructed instance.
Subclasses of SSA should provide an associated function pointer matching the
SSAFactory typedef */
class SSA;
typedef SSA * (*SSAFactory)(IceGrid &, IceBasalResistancePlasticLaw &, 
              IceFlowLaw &, EnthalpyConverter &, const NCConfigVariable &);


//! PISM's SSA solver
class SSA : public ShallowStressBalance
{
  friend class SSA_taud;
public:
  SSA(IceGrid &g, IceBasalResistancePlasticLaw &b, IceFlowLaw &i, EnthalpyConverter &e,
        const NCConfigVariable &c);

  SSAStrengthExtension *strength_extension;

  virtual ~SSA() { 
    deallocate();
    delete strength_extension;
  }

  virtual PetscErrorCode init(PISMVars &vars);

  virtual PetscErrorCode update(bool fast);

  virtual PetscErrorCode set_initial_guess(IceModelVec2V &guess);

  virtual PetscErrorCode stdout_report(string &result);

  virtual void add_vars_to_output(string keyword, set<string> &result);
  virtual PetscErrorCode define_variables(set<string> vars, const NCTool &nc, nc_type nctype);
  virtual PetscErrorCode write_variables(set<string> vars, string filename);
  virtual PetscErrorCode write_model_state(string filename);

  //! Add pointers to diagnostic quantities to a dictionary.
  virtual void get_diagnostics(map<string, PISMDiagnostic*> &dict);

protected:
  virtual PetscErrorCode allocate();

  virtual PetscErrorCode deallocate();

  virtual PetscErrorCode solve()  = 0;

  virtual PetscErrorCode compute_driving_stress(IceModelVec2V &taud);

  virtual PetscErrorCode compute_basal_frictional_heating(IceModelVec2S &result);

  virtual PetscErrorCode compute_D2(IceModelVec2S &result);

  virtual PetscErrorCode compute_maximum_velocity();

  IceModelVec2Mask *mask;
  IceModelVec2S *thickness, *tauc, *surface, *bed;
  IceModelVec2V taud, velocity_old;
  IceModelVec3 *enthalpy;

  string stdout_ssa;

  // objects used by the SSA solver (internally)
  DA  SSADA;                    // dof=2 DA (grid.da2 has dof=1)
  Vec SSAX;  // global vector for solution

  // profiling
  int event_ssa;
};

//! \brief Computes the driving stress (taud).
class SSA_taud : public PISMDiag<SSA>
{
public:
  SSA_taud(SSA *m, IceGrid &g, PISMVars &my_vars);
  PetscErrorCode compute(IceModelVec* &result);
};

#endif /* _SSA_H_ */

