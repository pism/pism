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

#include "pism/inverse/IP_BlatterForwardProblem.hh"

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
  term (effective viscosity `eta(B, gamma)`), so the design Jacobian is a
  **volume integral over the whole 3D sigma grid**. By contrast, the tauc
  design Jacobian is supported only on the basal face. Because the effective
  viscosity is linear in `B`, `J_design * dzeta` is simply the viscous
  residual evaluated with `B` replaced by `g'(zeta) * dzeta`, and
  `J_design^T * lambda` is the column sum (over all sigma levels) of the
  viscous "stress power" of the adjoint field weighted by `eta / B`.

  "Vertically averaged" hardness means a single column-constant value
  `B(i,j)` at each grid point, replicated across the Blatter sigma grid by
  the solver itself: linearize_at passes m_hardav as
  stressbalance::Inputs::averaged_hardness, and Blatter::init_ice_hardness
  uses it (via Blatter::init_averaged_ice_hardness) in place of the
  enthalpy-derived hardness. This matches the SSA-hardav semantics and the
  typical observational setting (where surface velocity is what we can fit).

  The basal yield stress is taken from the 2D field named `tauc` in
  `Grid::variables()`; this is how alternating tauc/hardav inversions hand
  the inverted tauc to the hardav phase.
*/
class IP_BlatterHardavForwardProblem : public IP_BlatterForwardProblem {
public:

  IP_BlatterHardavForwardProblem(std::shared_ptr<const Grid> grid,
                                 int Mz, int coarsening_factor,
                                 IPDesignVariableParameterization &tp);

  virtual ~IP_BlatterHardavForwardProblem() = default;

  virtual void set_design(array::Scalar &zeta);

protected:

  virtual void fill_inputs(stressbalance::Inputs &inputs);

  /// Apply J_design * dzeta in the full 3D Blatter state space (volume integral).
  virtual void apply_jacobian_design_3d(array::Scalar &dzeta, Vec result_3d);

  /// Apply J_design^T * lambda_3d to get a 2D dzeta (column sum of the
  /// transpose volume integral).
  virtual void apply_jacobian_design_transpose_3d(Vec lambda_3d, array::Scalar &dzeta);

  /// Storage for vertically-averaged hardness (the design variable in its
  /// natural physical units, as produced by m_design_param.convertToDesignVariable).
  array::Scalar m_hardav;
};

} // end of namespace inverse
} // end of namespace pism

#endif /* IP_BLATTERHARDAVFORWARDPROBLEM_HH */
