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

#include "pism/inverse/IP_BlatterForwardProblem.hh"

namespace pism {
namespace inverse {

//! Implements the forward problem of the map taking \f$\tau_c\f$ to the
//! corresponding solution of the Blatter stress balance.
/*!
  Analogous to IP_SSATaucForwardProblem but using the Blatter (higher-order)
  solver instead of the SSA. The design variable is \f$\tau_c\f$ (basal yield
  stress), parameterized by \f$\zeta\f$ via an IPDesignVariableParameterization.

  The state variable for inversion is the **2D surface velocity**, extracted
  from the 3D Blatter solution (see IP_BlatterForwardProblem).

  Since \f$\tau_c\f$ enters the Blatter residual only through the basal
  boundary condition, the design Jacobian \f$J_{\rm Design} = \partial
  \mathcal{R}/\partial \zeta\f$ is nonzero only at basal (bottom-of-column)
  nodes. However, the state Jacobian couples all vertical levels, so adjoint
  solves operate on the full 3D system.

  If a 2D field named `hardav` is present in `Grid::variables()` it is used
  as column-constant ice hardness (instead of the enthalpy-derived hardness);
  this is how alternating tauc/hardav inversions hand the inverted hardness
  to the tauc phase.
*/
class IP_BlatterTaucForwardProblem : public IP_BlatterForwardProblem {
public:

  //! Constructs from the same objects as Blatter, plus a specification of how
  //! \f$\tau_c\f$ is parameterized.
  IP_BlatterTaucForwardProblem(std::shared_ptr<const Grid> grid,
                               int Mz, int coarsening_factor,
                               IPDesignVariableParameterization &tp);

  virtual ~IP_BlatterTaucForwardProblem() = default;

  //! Selects nodes where \f$\tau_c\f$ (more specifically \f$\zeta\f$) should
  //! not be adjusted. Alias of set_design_fixed_locations().
  virtual void set_tauc_fixed_locations(array::Scalar &locations) {
    set_design_fixed_locations(locations);
  }

  //! Exposes the \f$\tau_c\f$ parameterization being used. Alias of design_param().
  virtual IPDesignVariableParameterization &tauc_param() {
    return design_param();
  }

  virtual void set_design(array::Scalar &zeta);

protected:

  virtual void fill_inputs(stressbalance::Inputs &inputs);

  /// Apply J_design * dzeta in the full 3D Blatter state space (basal face integral).
  virtual void apply_jacobian_design_3d(array::Scalar &dzeta, Vec result_3d);

  /// Apply J_design^T * lambda_3d to get a 2D dzeta (basal face integral).
  virtual void apply_jacobian_design_transpose_3d(Vec lambda_3d, array::Scalar &dzeta);

  /// Storage for tauc (avoids modifying fields obtained via pism::Vars)
  array::Scalar2 m_tauc_copy;
};

} // end of namespace inverse
} // end of namespace pism

#endif /* IP_BLATTERTAUCFORWARDPROBLEM_HH */
