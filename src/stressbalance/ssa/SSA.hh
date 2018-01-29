// Copyright (C) 2004--2018 Jed Brown, Ed Bueler and Constantine Khroulev
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

#ifndef _SSA_H_
#define _SSA_H_

#include "pism/stressbalance/ShallowStressBalance.hh"
#include "pism/util/IceModelVec2CellType.hh"

namespace pism {

class Geometry;

namespace stressbalance {

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
  SSAStrengthExtension(const Config &c);

  void set_notional_strength(double my_nuH);
  void set_min_thickness(double my_min_thickness);
  double get_notional_strength() const;
  double get_min_thickness() const;
private:
  double m_min_thickness, m_constant_nu;
};

//! Callback for constructing a new SSA subclass.  The caller is
//! responsible for deleting the newly constructed SSA.
/*! The factory idiom gives a way to implement runtime polymorphism for the 
  choice of SSA algorithm.  The factory is a function pointer that takes 
  all the arguments of an SSA constructor and returns a newly constructed instance.
  Subclasses of SSA should provide an associated function pointer matching the
  SSAFactory typedef */
class SSA;
typedef SSA * (*SSAFactory)(IceGrid::ConstPtr);


//! PISM's SSA solver.
/*!
  An object of this type solves equations for the vertically-constant horizontal
  velocity of ice that is sliding over land or is floating.  The equations are, in
  their clearest divergence form
  \f[ - \frac{\partial T_{ij}}{\partial x_j} - \tau_{(b)i} = f_i \f]
  where \f$i,j\f$ range over \f$x,y\f$, \f$T_{ij}\f$ is a depth-integrated viscous
  stress tensor (%i.e. equation (2.6) in [\ref SchoofStream]).
  These equations determine velocity in a more-or-less elliptic manner.
  Here \f$\tau_{(b)i}\f$ are the components of the basal shear stress applied to
  the base of the ice.  The right-hand side \f$f_i\f$ is the driving shear stress,
  \f[ f_i = - \rho g H \frac{\partial h}{\partial x_i}. \f]
  Here \f$H\f$ is the ice thickness and \f$h\f$ is the elevation of the surface of
  the ice.  More concretely, the SSA equations are
  \f{align*}
  - 2 \left[\nu H \left(2 u_x + v_y\right)\right]_x
  - \left[\nu H \left(u_y + v_x\right)\right]_y
  - \tau_{(b)1}  &= - \rho g H h_x, \\
  - \left[\nu H \left(u_y + v_x\right)\right]_x
  - 2 \left[\nu H \left(u_x + 2 v_y\right)\right]_y
  - \tau_{(b)2}  &= - \rho g H h_y, 
  \f}
  where \f$u\f$ is the \f$x\f$-component of the velocity and \f$v\f$ is the
  \f$y\f$-component of the velocity [\ref MacAyeal, \ref Morland, \ref WeisGreveHutter].

  Derived classes actually implement numerical methods to solve these equations.
  This class is virtual, but it actually implements some helper functions believed
  to be common to all implementations (%i.e. regular grid implementations) and it
  provides the basic fields.
*/
class SSA : public ShallowStressBalance {
public:
  SSA(IceGrid::ConstPtr g);
  virtual ~SSA();

  SSAStrengthExtension *strength_extension;

  virtual void update(const Inputs &inputs, bool full_update);

  void set_initial_guess(const IceModelVec2V &guess);

  virtual std::string stdout_report() const;

  const IceModelVec2V& driving_stress() const;
protected:
  virtual void define_model_state_impl(const PIO &output) const;
  virtual void write_model_state_impl(const PIO &output) const;

  virtual void init_impl();

  virtual DiagnosticList diagnostics_impl() const;

  virtual void compute_driving_stress(const Geometry &geometry, IceModelVec2V &result) const;

  virtual void solve(const Inputs &inputs) = 0;

  IceModelVec2CellType m_mask;
  IceModelVec2V m_taud;

  std::string m_stdout_ssa;

  // objects used by the SSA solver (internally)
  petsc::DM::Ptr  m_da;               // dof=2 DA
  IceModelVec2V m_velocity_global; // global vector for solution

  // profiling
  int m_event_ssa;
};

} // end of namespace stressbalance
} // end of namespace pism

#endif /* _SSA_H_ */
