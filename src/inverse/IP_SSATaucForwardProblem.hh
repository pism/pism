// Copyright (C) 2012, 2013, 2014, 2015, 2016, 2017, 2020, 2021, 2022 David Maxwell and Constantine Khroulev
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

#ifndef IP_SSATAUCFORWARDPROBLEM_HH_4AEVR4Z
#define IP_SSATAUCFORWARDPROBLEM_HH_4AEVR4Z

#include "pism/stressbalance/ssa/SSAFEM.hh"
#include "IPDesignVariableParameterization.hh"
#include "pism/util/petscwrappers/KSP.hh"
#include "pism/util/petscwrappers/Mat.hh"

namespace pism {
namespace inverse {

//! Implements the forward problem of the map taking \f$\tau_c\f$ to the corresponding solution of the %SSA.
/*! The class SSAFEM solves the %SSA, and the solution depends on a large number of parameters.  Considering
  all of these to be fixed except for \f$\tau_c\f$, we obtain a map \f$F_{\rm SSA}\f$ taking
  \f$\tau_c\f$ to the corresponding solution \f$u_{\rm SSA}\f$ of the %SSA.  This is a forward problem which
  we would like to invert; given \f$u_{\rm SSA}\f$ determine \f$\tau_c\f$  such that
  \f$F_{\rm SSA}(\tau_c) = u_{\rm SSA}\f$.  The forward problem actually implemented by
  IP_SSATaucForwardProblem is a mild variation \f$F_{\rm SSA}\f$.

  First, given the constraint \f$\tau_c\ge 0\f$ it can be helpful to parameterize \f$\tau_c\f$ by some other
  parameter \f$\zeta\f$,
  \f[
  \tau_c = g(\zeta),
  \f]
  where the function \f$g\f$ is non-negative.  The function \f$g\f$ is specified by an instance
  of IPDesignVariableParameterization.

  Second, there may be locations where the value of \f$\tau_c\f$ (and hence \f$\zeta\f$)
  is known a-priori, and should not be adjusted.  Let \f$\pi\f$ be the map that replaces
  the values of \f$\zeta\f$ with known values at these locations.

  IP_SSATaucForwardProblem implements the forward problem
  \f[
  F(\zeta) = F_{\rm SSA}(g(\pi(\zeta))).
  \f]

  Algorithms to solve the inverse problem make use of variations on the linearization
  of \f$F\f$, which are explained below.

  The solution of the %SSA in SSAFEM is implemented using SNES to solve
  a nonlinear residual function of the form
  \f[
  \mathcal{R}_{SSA}(u;\tau_c, \text{other parameters})= \vec 0,
  \f]
  usually using Newton's method. Let
  \f[
  \mathcal{R}(u, \zeta) = \mathcal{R}_{SSA}(u; g(\pi(\zeta)), \text{other parameters}).
  \f]

  The derivative of \f$\mathcal{R}\f$ with respect to \f$ u\f$ is called the state Jacobian,
  \f$J_{\rm State}\f$. Specifically, if \f$u=[U_j]\f$ then
  \f[
  (J_{\rm State})_{ij} = \frac{d \mathcal{R}_i}{dU_j}.
  \f]
  This is exactly the same Jacobian that is computed by SSAFEM for solving the %SSA via SNES. The
  matrix for the design Jacobian can be obtained with \ref assemble_jacobian_state.

  The derivative of \f$\mathcal{R}\f$ with respect to \f$ \zeta\f$ is called the design Jacobian,
  \f$J_{\rm Design}\f$. Specifically, if \f$\zeta=[Z_k]\f$ then
  \f[
  (J_{\rm Design})_{ik} = \frac{d \mathcal{R}_i}{dZ_k}.
  \f]
  The map \f$J_{\rm Design}\f$ can be applied to a vector \f$d\zeta\f$ with
  apply_jacobian_design.  For inverse methods using adjoints, one also
  needs to be able to apply the transpose of \f$J_{\rm Design}\f$,
  which is done using \ref apply_jacobian_design_transpose.

  The forward problem map \f$F\f$ solves the implicit equation
  \f[
  \mathcal{R}(F(\zeta), \zeta) = 0.
  \f]
  Its linearisation \f$DF\f$ then satisfies
  \f[
  J_{\rm State}\; DF\; d\zeta + J_{\rm Design}\; d\zeta = 0
  \f]
  for any perturbation \f$d\zeta\f$.  Hence
  \f[
  DF = -J_{\rm State}^{-1}\circ J_{\rm Design}.
  \f]
  This derivative is sometimes called the reduced gradient in the literature.  To
  apply \f$DF\f$ to a perturbation \f$d\zeta\f$, use \ref apply_linearization.  Adjoint
  methods require the transpose of this map; to apply \f$DF^t\f$ to \f$du\f$ use
  \ref apply_linearization_transpose.
*/
class IP_SSATaucForwardProblem : public stressbalance::SSAFEM
{
public:

  /// The function space for the design variable, i.e. \f$\tau_c\f$.
  typedef array::Scalar DesignVec;
  typedef array::Scalar1 DesignVecGhosted;

  /// The function space for the state variable, \f$u_{\rm SSA}\f$.
  typedef IceModelVec2V StateVec;
  typedef Velocity1 StateVec1;

  //! Constructs from the same objects as SSAFEM, plus a specification of how \f$\tau_c\f$
  //! is parameterized.
  IP_SSATaucForwardProblem(IceGrid::ConstPtr g,
                           IPDesignVariableParameterization &tp);

  virtual ~IP_SSATaucForwardProblem() = default;

  void init();

  //! Selects nodes where \f$\tau_c\f$ (more specifically \f$\zeta\f$) should not be adjusted.
  /*! The paramter \a locations should be set to 1 at each node where \f$\tau_c\f$
    is fixed. The forward map then effectively treats the design space as the subspace
    of nodes where \a locations is 0. Tangent vectors to this subspace, as would be
    generated by, e.g., \f$J_{\rm Design}^t\f$ are represented as vectors in the full
    space with entries set to zero in the fixed locations.  These can safely be added
    to preexisting values of \f$\zeta\f$ without changing the entries of \f$\zeta\f$ at the
    fixed locations.  Inversion can be done by setting an initial value of \f$\zeta\f$
    having the desired values in the fixed locations, and using set_tauc_fixed_locations()
    to indicate the nodes that should not be changed.
  */
  virtual void set_tauc_fixed_locations(array::Scalar &locations)
  {
    m_fixed_tauc_locations = &locations;
  }

  //! Returns the last solution of the %SSA as computed by \ref linearize_at.
  virtual IceModelVec2V::Ptr solution() {
    m_velocity_shared->copy_from(m_velocity);
    return m_velocity_shared;
  }

  //! Exposes the \f$\tau_c\f$ parameterization being used.
  virtual IPDesignVariableParameterization & tauc_param() {
    return m_tauc_param;
  }

  virtual void set_design(array::Scalar &zeta);

  virtual TerminationReason::Ptr linearize_at(array::Scalar &zeta);

  virtual void assemble_residual(IceModelVec2V &u, IceModelVec2V &R);
  virtual void assemble_residual(IceModelVec2V &u, Vec R);

  virtual void assemble_jacobian_state(IceModelVec2V &u, Mat J);

  virtual void apply_jacobian_design(IceModelVec2V &u, array::Scalar &dzeta, IceModelVec2V &du);
  virtual void apply_jacobian_design(IceModelVec2V &u, array::Scalar &dzeta, Vec du);
  virtual void apply_jacobian_design(IceModelVec2V &u, array::Scalar &dzeta, Vector2 **du_a);

  virtual void apply_jacobian_design_transpose(IceModelVec2V &u, IceModelVec2V &du, array::Scalar &dzeta);
  virtual void apply_jacobian_design_transpose(IceModelVec2V &u, IceModelVec2V &du, Vec dzeta);
  virtual void apply_jacobian_design_transpose(IceModelVec2V &u, IceModelVec2V &du, double **dzeta);

  virtual void apply_linearization(array::Scalar &dzeta, IceModelVec2V &du);
  virtual void apply_linearization_transpose(IceModelVec2V &du, array::Scalar &dzeta);

  //! Exposes the DMDA of the underlying grid for the benefit of TAO.
  petsc::DM& get_da() const {
    return *m_da;
  }

protected:

  /// Current value of zeta, provided from caller.
  array::Scalar   *m_zeta;
  /// Storage for d_zeta with ghosts, if needed when an argument d_zeta is ghost-less.
  array::Scalar1   m_dzeta_local;
  /// Storage for tauc (avoids modifying fields obtained via pism::Vars)
  array::Scalar2 m_tauc_copy;

  /// Locations where \f$\tau_c\f$ should not be adjusted.
  array::Scalar *m_fixed_tauc_locations;

  /// The function taking \f$\zeta\f$ to \f$\tau_c\f$.
  IPDesignVariableParameterization &m_tauc_param;

  /// Copy of the velocity field managed using a shared pointer.
  IceModelVec2V::Ptr m_velocity_shared;

  /// Temporary storage when state vectors need to be used without ghosts.
  IceModelVec2V  m_du_global;
  /// Temporary storage when state vectors need to be used with ghosts.
  Velocity1  m_du_local;

  fem::ElementIterator m_element_index;
  fem::Q1Element2       m_element;

  /// KSP used in \ref apply_linearization and \ref apply_linearization_transpose
  petsc::KSP  m_ksp;
  /// Mat used in \ref apply_linearization and \ref apply_linearization_transpose
  petsc::Mat  m_J_state;

  SNESConvergedReason m_reason;

  /// Flag indicating that the state jacobian matrix needs rebuilding.
  bool m_rebuild_J_state;
};

} // end of namespace inverse
} // end of namespace pism

#endif /* end of include guard: IP_SSATAUCFORWARDPROBLEM_HH_4AEVR4Z */
