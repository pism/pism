// Copyright (C) 2012  David Maxwell
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

#ifndef INVSSAFORWARDPROBLEM_HH_4AEVR4Z
#define INVSSAFORWARDPROBLEM_HH_4AEVR4Z


#include "SSAFEM.hh"
#include "InvTaucParameterization.hh"

class InvSSAForwardProblem : public SSAFEM
{
public:

  InvSSAForwardProblem(IceGrid &g, IceBasalResistancePlasticLaw &b,
    EnthalpyConverter &e, InvTaucParameterization &tp,
    const NCConfigVariable &c);

  virtual ~InvSSAForwardProblem();

  virtual PetscErrorCode set_tauc_fixed_locations(IceModelVec2Int &locations)
  { 
    m_fixed_tauc_locations = &locations;
    return 0;
  }

  IceModelVec2V &solution() {
    return m_velocity;
  }

  InvTaucParameterization & tauc_param() {
    return m_tauc_param;
  }

  PetscErrorCode set_zeta( IceModelVec2S &zeta);

  PetscErrorCode linearize_at( IceModelVec2S &zeta, TerminationReason::Ptr &reason);

  PetscErrorCode assemble_residual(IceModelVec2V &u, IceModelVec2V &R);
  PetscErrorCode assemble_residual(IceModelVec2V &u, Vec R);

  PetscErrorCode assemble_jacobian_state(IceModelVec2V &u, Mat J);

  PetscErrorCode apply_jacobian_design(IceModelVec2V &u,IceModelVec2S &dzeta,IceModelVec2V &du);
  PetscErrorCode apply_jacobian_design(IceModelVec2V &u,IceModelVec2S &dzeta, Vec du);
  PetscErrorCode apply_jacobian_design(IceModelVec2V &u,IceModelVec2S &dzeta, PISMVector2 **du_a);

  PetscErrorCode apply_jacobian_design_transpose(IceModelVec2V &u,IceModelVec2V &du,IceModelVec2S &dzeta);
  PetscErrorCode apply_jacobian_design_transpose(IceModelVec2V &u,IceModelVec2V &du,Vec dzeta);
  PetscErrorCode apply_jacobian_design_transpose(IceModelVec2V &u,IceModelVec2V &du,PetscReal **dzeta);

  PetscErrorCode apply_linearization(IceModelVec2S &dzeta, IceModelVec2V &du);
  PetscErrorCode apply_linearization_transpose(IceModelVec2V &du, IceModelVec2S &dzeta);

  PetscErrorCode get_da(DM *da) {
    *da = SSADA;
    return 0;
  }
  
protected:

  PetscErrorCode construct();
  PetscErrorCode destruct();

  IceGrid &m_grid;

  IceModelVec2S   *m_zeta;
  IceModelVec2S   m_dzeta_local;

  IceModelVec2Int *m_fixed_tauc_locations;

  InvTaucParameterization &m_tauc_param;

  IceModelVec2V  m_du_global;
  IceModelVec2V  m_du_local;

  FEElementMap m_element_index;
  FEQuadrature m_quadrature;
  FEDOFMap     m_dofmap;

  KSP  m_ksp;
  Mat  m_J_state;

  SNESConvergedReason m_reason;

  bool m_rebuild_J_state;
};


#endif /* end of include guard: INVSSAFORWARDPROBLEM_HH_4AEVR4Z */
