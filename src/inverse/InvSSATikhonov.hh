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

#ifndef INVSSATIKHONOV_HH_NSO6ITYB
#define INVSSATIKHONOV_HH_NSO6ITYB

#include "SSAFEM.hh"
#include "InvTaucParameterization.hh"
#include "Functional.hh"
#include <memory> // for std::auto_ptr
#include "TikhonovProblem.hh"
#include "PythonTikhonovSVListener.hh"
#include "MeanSquareFunctional.hh"

class InvSSATikhonov : public SSAFEM
{
public:

  typedef IceModelVec2S DesignVec;
  typedef IceModelVec2V StateVec;

  InvSSATikhonov(IceGrid &g, IceBasalResistancePlasticLaw &b,
    EnthalpyConverter &e, InvTaucParameterization &tp,
    const NCConfigVariable &c);

  virtual ~InvSSATikhonov();

  virtual PetscErrorCode init(PISMVars &vars);

  virtual PetscErrorCode set_tauc_fixed_locations(IceModelVec2Int &locations)
  { 
    m_fixed_tauc_locations = &locations;
    return 0;
  }

  IceModelVec2V &solution() {
    return velocity;
  }

  std::string reasonDescription() {
    return std::string("FIXME");
  }
  
  PetscErrorCode get_da(DM *da) {
    *da = SSADA;
    return 0;
  }

  PetscErrorCode set_functionals();

  PetscErrorCode getVariableBounds(Vec lo, Vec hi);

  PetscErrorCode set_zeta( IceModelVec2S &zeta);

  PetscErrorCode linearizeAt( IceModelVec2S &c, bool &success);

  PetscErrorCode evalObjective(IceModelVec2S &dc, PetscReal *OUTPUT);
  PetscErrorCode evalGradObjective(IceModelVec2S &dc, IceModelVec2S &gradient);

  PetscErrorCode evalPenalty(IceModelVec2V &du, PetscReal *OUTPUT);
  PetscErrorCode evalGradPenalty(IceModelVec2V &du,  IceModelVec2V &gradient);
  PetscErrorCode evalGradPenaltyReduced(IceModelVec2V &dr, IceModelVec2S &gradient);
  PetscErrorCode evalGradPenaltyReducedFD(IceModelVec2V &du, IceModelVec2S &gradient);
  PetscErrorCode evalGradPenaltyReducedNoTranspose(IceModelVec2V &du, IceModelVec2S &gradient);

  PetscErrorCode assembleFunction(IceModelVec2V u, Vec RHS);
  
  PetscErrorCode assembleJacobian(IceModelVec2V u, Mat J);

  PetscErrorCode assemble_T_rhs(IceModelVec2S &zeta, IceModelVec2V &u, IceModelVec2S &dzeta, IceModelVec2V &rhs);
  PetscErrorCode assemble_T_rhs(PetscReal **zeta, PISMVector2 **u, PetscReal **dzeta, PISMVector2 **rhs);
  PetscErrorCode computeT(IceModelVec2S &dzeta,IceModelVec2V &du);
  PetscErrorCode compute_Jdesign_transpose(PetscReal **zeta_a, PISMVector2 **u_a, PISMVector2 **v_a, PetscReal **gradient_a);
  PetscErrorCode domainIP(IceModelVec2S &a, IceModelVec2S &b, PetscReal *OUTPUT);
  PetscErrorCode rangeIP(IceModelVec2V &a, IceModelVec2V &b, PetscReal *OUTPUT);

protected:


  PetscErrorCode construct();
  PetscErrorCode destruct();

  IceGrid &m_grid;
  IceModelVec2S   *m_zeta;

  IceModelVec2Int *m_fixed_tauc_locations;
  IceModelVec2S   *m_misfit_weight;

  InvTaucParameterization &m_tauc_param;
  
  IceModelVec2V  m_vGlobal;
  IceModelVec2V  m_v;
  IceModelVec2V  m_adjointRHS;
  
  FEElementMap m_element_index;
  FEQuadrature m_quadrature;
  FEDOFMap     m_dofmap;

  std::auto_ptr< Functional<IceModelVec2S> > m_designFunctional;
  std::auto_ptr< Functional<IceModelVec2V> > m_penaltyFunctional;

  std::auto_ptr< MeanSquareFunctional2S > m_domainIP;

  KSP  m_ksp;
  Mat  m_Jadjoint;
  
  SNESConvergedReason m_reason;
};

// SWIG is having problem recognising subclasses of templated superclasses
// have the correct type.  The following are workarounds to allow for
// python-based listeners to InvSSATikhonov

class InvSSAPythonListenerBridge: public TikhonovProblemListener<InvSSATikhonov> {
public:
  InvSSAPythonListenerBridge(PythonTikhonovSVListener::Ptr core) : m_core(core) { }
  PetscErrorCode iteration( TikhonovProblem<InvSSATikhonov> &,
             PetscReal eta, PetscInt iter, 
             PetscReal objectiveValue, PetscReal designValue,
             IceModelVec2S &d, IceModelVec2S &diff_d, IceModelVec2S &grad_d,
             IceModelVec2V &u,  IceModelVec2V &diff_u,  IceModelVec2S &grad_u,
             IceModelVec2S &gradient) {
    m_core->iteration(iter,eta,objectiveValue,designValue,d,diff_d,grad_d,
      u,diff_u,grad_u,gradient);
    return 0;
  }
protected:
  PythonTikhonovSVListener::Ptr m_core;
};

void InvSSATikhonovAddListener(TikhonovProblem<InvSSATikhonov> &problem, 
PythonTikhonovSVListener::Ptr listener );

#endif /* end of include guard: INVSSATIKHONOV_HH_NSO6ITYB */
