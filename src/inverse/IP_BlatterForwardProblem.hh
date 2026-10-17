// Copyright (C) 2026 Andy Aschwanden and Constantine Khroulev
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

#ifndef IP_BLATTERFORWARDPROBLEM_HH
#define IP_BLATTERFORWARDPROBLEM_HH

#include "pism/stressbalance/blatter/Blatter.hh"
#include "pism/inverse/IPDesignVariableParameterization.hh"
#include "pism/util/TerminationReason.hh"
#include "pism/util/petscwrappers/KSP.hh"
#include "pism/util/petscwrappers/Mat.hh"

namespace pism {
namespace inverse {

//! Design-variable-agnostic part of a Blatter inverse forward problem.
/*!
  A Blatter forward problem for inversion maps a 2D design variable
  (parameterized by \f$\zeta\f$) to the 2D **surface** velocity extracted
  from the 3D Blatter solution. Everything that does not depend on *which*
  design variable is used lives here:

  - construction of the `stressbalance::Inputs` from `Grid::variables()` and
    the forward SNES solve (linearize_at),
  - extraction of the surface velocity (the operator \f$P\f$ in the
    documentation),
  - the reduced linearization \f$DF = -P J_{\rm State}^{-1} J_{\rm Design}\f$
    (apply_linearization) and its transpose (apply_linearization_transpose),
    including the adjoint KSP and the three adjoint methods selected by
    `inverse.adjoint.method`.

  Derived classes supply the design-variable-specific pieces:

  - set_design(): convert \f$\zeta\f$ to the physical design variable,
  - fill_inputs(): hand the design variable to the Blatter solver,
  - apply_jacobian_design_3d() / apply_jacobian_design_transpose_3d(): the
    design Jacobian \f$J_{\rm Design} = \partial \mathcal{R}/\partial \zeta\f$
    of the 3D residual and its transpose.

  See IP_BlatterTaucForwardProblem (basal yield stress; the design Jacobian
  is a basal-face integral) and IP_BlatterHardavForwardProblem (vertically
  averaged ice hardness; the design Jacobian is a volume integral over the
  whole ice column).
*/
class IP_BlatterForwardProblem : public stressbalance::Blatter {
public:

  /// The function space for the design variable.
  typedef array::Scalar DesignVec;
  typedef array::Scalar1 DesignVecGhosted;

  /// The function space for the state variable (2D surface velocity).
  typedef array::Vector StateVec;
  typedef array::Vector1 StateVec1;

  IP_BlatterForwardProblem(std::shared_ptr<const Grid> grid,
                           int Mz, int coarsening_factor,
                           IPDesignVariableParameterization &tp);

  virtual ~IP_BlatterForwardProblem() = default;

  void init();

  //! Selects nodes where the design variable (more specifically \f$\zeta\f$)
  //! should not be adjusted.
  virtual void set_design_fixed_locations(array::Scalar &locations) {
    m_fixed_design_locations = &locations;
  }

  //! Returns the 2D surface velocity from the last Blatter solve.
  virtual std::shared_ptr<array::Vector> solution() {
    return m_surface_velocity;
  }

  //! Exposes the design-variable parameterization in use.
  virtual IPDesignVariableParameterization &design_param() {
    return m_design_param;
  }

  //! Sets the current value of \f$\zeta\f$ and derives the physical design
  //! variable from it. Implementations must call store_zeta().
  virtual void set_design(array::Scalar &zeta) = 0;

  //! Sets \f$\zeta\f$ and solves the Blatter system.
  virtual std::shared_ptr<TerminationReason> linearize_at(array::Scalar &zeta);

  //! Not supported for Blatter (the design Jacobian lives in the 3D state
  //! space); use apply_linearization() instead.
  virtual void apply_jacobian_design(array::Vector &u,
                                     array::Scalar &dzeta,
                                     array::Vector &du);

  //! Not supported for Blatter; use apply_linearization_transpose() instead.
  virtual void apply_jacobian_design_transpose(array::Vector &u,
                                               array::Vector &du,
                                               array::Scalar &dzeta);

  //! Applies the reduced linearization \f$DF\f$ to a design perturbation.
  virtual void apply_linearization(array::Scalar &dzeta, array::Vector &du);

  //! Applies the transpose of the reduced linearization \f$DF^T\f$ to a
  //! (surface velocity) state perturbation.
  virtual void apply_linearization_transpose(array::Vector &du,
                                             array::Scalar &dzeta);

  //! Exposes the DM for the benefit of TAO.
  petsc::DM &get_da() const;

protected:

  //! Hand the design variable (and anything else design-specific) to the
  //! Blatter solver. Called by linearize_at() after the geometry, enthalpy,
  //! `tauc` (when available in Grid::variables()) and boundary conditions
  //! have been filled in.
  virtual void fill_inputs(stressbalance::Inputs &inputs) = 0;

  //! Apply \f$J_{\rm Design}\, d\zeta\f$ in the full 3D Blatter state space.
  //! `dzeta` is ghosted, with fixed design locations already zeroed.
  virtual void apply_jacobian_design_3d(array::Scalar &dzeta, Vec result_3d) = 0;

  //! Apply \f$J_{\rm Design}^T \lambda\f$ (3D adjoint) to get a 2D `dzeta`.
  //! Implementations are responsible for multiplying by \f$g'(\zeta)\f$ and
  //! zeroing fixed design locations.
  virtual void apply_jacobian_design_transpose_3d(Vec lambda_3d, array::Scalar &dzeta) = 0;

  //! Extract the 2D surface velocity (top sigma level) from the 3D solution.
  void extract_surface_velocity();

  //! Point m_zeta at a ghosted view of `zeta` (copying if necessary), so
  //! design Jacobian kernels can read \f$\zeta\f$ at ghost-element nodes.
  void store_zeta(array::Scalar &zeta);

  //! Re-assemble the SNES (Newton) Jacobian at the current solution if needed.
  void update_state_jacobian();

  //! Zero `dzeta` at fixed design locations (owned points only).
  void apply_fixed_locations(array::Scalar &dzeta);

  /// Current value of zeta (ghosted), provided by the caller or copied.
  array::Scalar *m_zeta;

  /// Ghosted copy of zeta, used when the caller provides a ghostless one.
  array::Scalar1 m_zeta_local;

  /// Ghosted copy of dzeta with fixed locations zeroed.
  array::Scalar1 m_dzeta_local;

  /// Locations where the design variable should not be adjusted.
  array::Scalar *m_fixed_design_locations;

  /// The function taking \f$\zeta\f$ to the physical design variable.
  IPDesignVariableParameterization &m_design_param;

  /// 2D surface velocity extracted from the 3D Blatter solution.
  std::shared_ptr<array::Vector> m_surface_velocity;

  /// Picard (symmetric) Jacobian matrix used for the "incomplete" adjoint.
  petsc::Mat m_J_picard;

  /// Standalone KSP for adjoint solves (option prefix `inv_adj_`).
  petsc::KSP m_ksp;

  /// Flag indicating that the state Jacobian needs re-assembly.
  bool m_rebuild_J_state;
};

} // end of namespace inverse
} // end of namespace pism

#endif /* IP_BLATTERFORWARDPROBLEM_HH */
