// Copyright (C) 2009--2017, 2020, 2021, 2022 Jed Brown and Ed Bueler and Constantine Khroulev and David Maxwell
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

#ifndef _SSAFEM_H_
#define _SSAFEM_H_

#include "SSA.hh"
#include "pism/util/fem/FEM.hh"
#include "pism/util/petscwrappers/SNES.hh"
#include "pism/util/TerminationReason.hh"
#include "pism/util/Mask.hh"

namespace pism {

namespace stressbalance {

//! Factory function for constructing a new SSAFEM.
SSA * SSAFEMFactory(IceGrid::ConstPtr grid);

//! PISM's SSA solver: the finite element method implementation written by Jed and David
/*!
  Jed's original code is in rev 831: src/base/ssaJed/...
  The SSAFEM duplicates most of the functionality of SSAFD, using the finite element method.
*/
class SSAFEM : public SSA {
public:
  SSAFEM(IceGrid::ConstPtr g);

  virtual ~SSAFEM() = default;

protected:
  virtual void init_impl();
  void cache_inputs(const Inputs &inputs);

  //! Storage for SSA coefficients at element nodes.
  //!
  //! All fields must be "double" or structures containing "double"
  //! for IceModelVec2 to work correctly.
  struct Coefficients {
    //! ice thickness
    double thickness;
    //! bed elevation
    double bed;
    //! sea level
    double sea_level;
    //! basal yield stress
    double tauc;
    //! ice hardness
    double hardness;
    //! prescribed gravitational driving stress
    Vector2 driving_stress;
  };

  Array2SGhosted<1> m_bc_mask;
  IceModelVec2V m_bc_values;

  GeometryCalculator m_gc;
  double m_alpha;
  double m_rho_g;

  IceModelVec2<Coefficients> m_coefficients;

  void quad_point_values(const fem::Element &Q,
                         const Coefficients *x,
                         int *mask,
                         double *thickness,
                         double *tauc,
                         double *hardness) const;

  void explicit_driving_stress(const fem::Element &E,
                               const Coefficients *x,
                               Vector2 *driving_stress) const;

  void driving_stress(const fem::Element &E,
                      const Coefficients *x,
                      Vector2 *driving_stress) const;

  void PointwiseNuHAndBeta(double thickness,
                           double hardness,
                           int mask,
                           double tauc,
                           const Vector2 &U,
                           const Vector2 &U_x,
                           const Vector2 &U_y,
                           double *nuH, double *dnuH,
                           double *beta, double *dbeta);

  void compute_local_function(Vector2 const *const *const velocity,
                              Vector2 **residual);

  void compute_local_jacobian(Vector2 const *const *const velocity, Mat J);

  virtual void solve(const Inputs &inputs);

  TerminationReason::Ptr solve_with_reason(const Inputs &inputs);

  TerminationReason::Ptr solve_nocache();

  //! Adaptor for gluing SNESDAFormFunction callbacks to an SSAFEM.
  /* The callbacks from SNES are mediated via SNESDAFormFunction, which has the
     convention that its context argument is a pointer to a struct
     having a DA as its first entry.  The CallbackData fulfills
     this requirement, and allows for passing the callback on to an honest
     SSAFEM object. */
  struct CallbackData {
    DM da;
    SSAFEM *ssa;
  };

  // objects used internally
  CallbackData m_callback_data;

  petsc::SNES m_snes;

  //! Storage for node types (interior, boundary, exterior).
  Array2SGhosted<1> m_node_type;
  //! Boundary integral (CFBC contribution to the residual).
  IceModelVec2V m_boundary_integral;

  double m_dirichletScale;
  double m_beta_ice_free_bedrock;
  double m_epsilon_ssa;

  fem::ElementIterator m_element_index;
  fem::Q1Element2 m_q1_element;
  // fem::P1Element m_p1_element;

  // Support for direct specification of driving stress to the FEM SSA solver. This helps
  // with certain test cases where the grid is periodic but the driving stress cannot be the
  // gradient of a periodic function. (See commit ffb4be16.)
  const IceModelVec2S *m_driving_stress_x;
  const IceModelVec2S *m_driving_stress_y;
private:
  void cache_residual_cfbc(const Inputs &inputs);
  void monitor_jacobian(Mat Jac);
  void monitor_function(Vector2 const *const *const velocity_global,
                        Vector2 const *const *const residual_global);

  //! SNES callbacks.
  /*! These simply forward the call on to the SSAFEM member of the CallbackData */
  static PetscErrorCode function_callback(DMDALocalInfo *info,
                                          Vector2 const *const *const velocity,
                                          Vector2 **residual,
                                          CallbackData *fe);
  static PetscErrorCode jacobian_callback(DMDALocalInfo *info,
                                          Vector2 const *const *const xg,
                                          Mat A, Mat J, CallbackData *fe);
};


} // end of namespace stressbalance
} // end of namespace pism

#endif /* _SSAFEM_H_ */
