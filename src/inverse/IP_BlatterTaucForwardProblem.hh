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

#ifndef IP_BLATTERTAUCFORWARDPROBLEM_HH
#define IP_BLATTERTAUCFORWARDPROBLEM_HH

#include "pism/stressbalance/blatter/Blatter.hh"
#include "pism/inverse/IPDesignVariableParameterization.hh"
#include "pism/util/TerminationReason.hh"
#include "pism/util/petscwrappers/KSP.hh"

namespace pism {
namespace inverse {

//! Implements the forward problem of the map taking \f$\tau_c\f$ to the
//! corresponding solution of the Blatter stress balance.
/*!
  Analogous to IP_SSATaucForwardProblem but using the Blatter (higher-order)
  solver instead of the SSA. The design variable is \f$\tau_c\f$ (basal yield
  stress), parameterized by \f$\zeta\f$ via an IPDesignVariableParameterization.

  The state variable for inversion is the **2D surface velocity**, extracted
  from the 3D Blatter solution. This matches the typical observation space
  (satellite-derived surface velocities).

  The Blatter solver produces 3D velocity on a sigma grid internally. The
  forward problem extracts surface velocity and provides the necessary
  linearization operators (design Jacobian, its transpose, and the reduced
  gradient) for adjoint-based inversion.

  Since \f$\tau_c\f$ enters the Blatter residual only through the basal
  boundary condition, the design Jacobian \f$J_{\rm Design} = \partial
  \mathcal{R}/\partial \zeta\f$ is nonzero only at basal (bottom-of-column)
  nodes. However, the state Jacobian couples all vertical levels, so adjoint
  solves operate on the full 3D system.
*/
class IP_BlatterTaucForwardProblem : public stressbalance::Blatter {
public:

  /// The function space for the design variable, i.e. \f$\tau_c\f$.
  typedef array::Scalar DesignVec;
  typedef array::Scalar1 DesignVecGhosted;

  /// The function space for the state variable (2D surface velocity).
  typedef array::Vector StateVec;
  typedef array::Vector1 StateVec1;

  //! Constructs from the same objects as Blatter, plus a specification of how
  //! \f$\tau_c\f$ is parameterized.
  IP_BlatterTaucForwardProblem(std::shared_ptr<const Grid> grid,
                               int Mz, int coarsening_factor,
                               IPDesignVariableParameterization &tp);

  virtual ~IP_BlatterTaucForwardProblem() = default;

  void init();

  //! Selects nodes where \f$\tau_c\f$ (more specifically \f$\zeta\f$) should
  //! not be adjusted.
  virtual void set_tauc_fixed_locations(array::Scalar &locations) {
    m_fixed_tauc_locations = &locations;
  }

  //! Returns the 2D surface velocity from the last Blatter solve.
  virtual std::shared_ptr<array::Vector> solution() {
    return m_surface_velocity;
  }

  //! Exposes the \f$\tau_c\f$ parameterization being used.
  virtual IPDesignVariableParameterization &tauc_param() {
    return m_tauc_param;
  }

  virtual void set_design(array::Scalar &zeta);

  virtual std::shared_ptr<TerminationReason> linearize_at(array::Scalar &zeta);

  virtual void apply_jacobian_design(array::Vector &u,
                                     array::Scalar &dzeta,
                                     array::Vector &du);

  virtual void apply_jacobian_design_transpose(array::Vector &u,
                                               array::Vector &du,
                                               array::Scalar &dzeta);

  virtual void apply_linearization(array::Scalar &dzeta, array::Vector &du);
  virtual void apply_linearization_transpose(array::Vector &du,
                                             array::Scalar &dzeta);

  //! Exposes the DM for the benefit of TAO.
  petsc::DM &get_da() const;

protected:

  /// Extract 2D surface velocity from the 3D Blatter solution.
  void extract_surface_velocity();

  /// Apply J_design * dzeta in the full 3D Blatter state space.
  void apply_jacobian_design_3d(array::Scalar &dzeta, Vec result_3d);

  /// Apply J_design^T * lambda_3d to get a 2D dzeta.
  void apply_jacobian_design_transpose_3d(Vec lambda_3d, array::Scalar &dzeta);

  /// Current value of zeta, provided from caller.
  array::Scalar *m_zeta;

  /// Storage for tauc (avoids modifying fields obtained via pism::Vars)
  array::Scalar2 m_tauc_copy;

  /// Locations where \f$\tau_c\f$ should not be adjusted.
  array::Scalar *m_fixed_tauc_locations;

  /// The function taking \f$\zeta\f$ to \f$\tau_c\f$.
  IPDesignVariableParameterization &m_tauc_param;

  /// 2D surface velocity extracted from the 3D Blatter solution.
  std::shared_ptr<array::Vector> m_surface_velocity;

  /// KSP for adjoint (transpose) solves, separate from the SNES's KSP
  /// because the MG+SOR smoother doesn't support transpose. Configured
  /// via -inv_ksp_type, -inv_pc_type, -inv_ksp_rtol, etc.
  petsc::KSP m_ksp;

  /// Flag indicating that the state Jacobian matrix needs rebuilding.
  bool m_rebuild_J_state;
};

} // end of namespace inverse
} // end of namespace pism

#endif /* IP_BLATTERTAUCFORWARDPROBLEM_HH */
