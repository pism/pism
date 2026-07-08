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

#ifndef IP_BLATTERHARDAVFORWARDPROBLEM_HH
#define IP_BLATTERHARDAVFORWARDPROBLEM_HH

#include "pism/stressbalance/blatter/Blatter.hh"
#include "pism/inverse/IPDesignVariableParameterization.hh"
#include "pism/util/TerminationReason.hh"
#include "pism/util/petscwrappers/KSP.hh"
#include "pism/util/petscwrappers/Mat.hh"

namespace pism {
namespace inverse {

//! Implements the forward problem of the map taking ice hardness (vertically
//! averaged Glen-law `B = A^{-1/n}`) to the corresponding solution of the
//! Blatter stress balance.
/*!
  Analogous to IP_SSAHardavForwardProblem but using the Blatter higher-order
  solver, and analogous to IP_BlatterTaucForwardProblem but with hardness as
  the design variable instead of basal yield stress.

  The hardness enters the Blatter residual only through the volume viscous
  term (effective viscosity `eta(B, gamma_dot)`), so the design Jacobian is a
  volume integral over the 3D sigma grid. By contrast, the tauc design
  Jacobian is supported only on the basal face. The state-Jacobian and
  surface-extraction machinery is identical to the tauc version.

  As in SSA, "vertically averaged" hardness means a single column-constant
  value `B(i,j)` at each grid point, replicated across the Blatter sigma
  grid at solve time. This matches the SSA-hardav semantics and the typical
  observational setting (where surface velocity is what we can fit).
*/
class IP_BlatterHardavForwardProblem : public stressbalance::Blatter {
public:

  /// The function space for the design variable (vertically-averaged hardness).
  typedef array::Scalar DesignVec;
  typedef array::Scalar1 DesignVecGhosted;

  /// The function space for the state variable (2D surface velocity).
  typedef array::Vector StateVec;
  typedef array::Vector1 StateVec1;

  IP_BlatterHardavForwardProblem(std::shared_ptr<const Grid> grid,
                                 int Mz, int coarsening_factor,
                                 IPDesignVariableParameterization &tp);

  virtual ~IP_BlatterHardavForwardProblem() = default;

  void init();

  //! Selects nodes where the design variable (more specifically `zeta`) should
  //! not be adjusted.
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

  /// Apply J_design * dzeta in the full 3D Blatter state space (volume integral).
  void apply_jacobian_design_3d(array::Scalar &dzeta, Vec result_3d);

  /// Apply J_design^T * lambda_3d to get a 2D dzeta (transpose volume integral).
  void apply_jacobian_design_transpose_3d(Vec lambda_3d, array::Scalar &dzeta);

  /// Overwrite the column-constant hardness stored on the Blatter sigma grid
  /// with the current value of m_hardav at each (i,j). Called inside
  /// linearize_at after Blatter::update() to ensure subsequent residual /
  /// Jacobian evaluations use our hardness rather than the enthalpy-derived
  /// one. (Blatter::init_ice_hardness is not virtual; this is the override
  /// hook.)
  void install_hardav_on_sigma_grid();

  /// Current value of zeta, provided from caller.
  array::Scalar *m_zeta;

  /// Storage for vertically-averaged hardness (the "design variable" in the
  /// natural physical units, as produced by m_design_param.convertToDesignVariable).
  array::Scalar2 m_hardav;

  /// Locations where the design variable should not be adjusted.
  array::Scalar *m_fixed_design_locations;

  /// The function taking `zeta` to hardness.
  IPDesignVariableParameterization &m_design_param;

  /// 2D surface velocity extracted from the 3D Blatter solution.
  std::shared_ptr<array::Vector> m_surface_velocity;

  /// Picard (symmetric) Jacobian matrix used for the adjoint solve.
  petsc::Mat m_J_picard;

  /// KSP for adjoint solves using the Picard Jacobian.
  petsc::KSP m_ksp;

  /// Flag indicating that the Jacobians need rebuilding.
  bool m_rebuild_J_state;
};

} // end of namespace inverse
} // end of namespace pism

#endif /* IP_BLATTERHARDAVFORWARDPROBLEM_HH */
