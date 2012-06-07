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

#ifndef INVSSABASICTIKHONOV_HH_5AJEULAC
#define INVSSABASICTIKHONOV_HH_5AJEULAC


#include "SNESCallbacks.hh"
#include "iceModelVec.hh"
#include "FETools.hh"
#include "TikhonovProblem.hh"
#include "TaoUtil.hh"
#include "Functional.hh"
#include "InvTaucParameterization.hh"
#include <memory> // for std::auto_ptr
#include "PythonTikhonovSVListener.hh"

class InvSSABasicTikhonov {
public:
  typedef IceModelVec2S DesignVec;
  typedef IceModelVec2V StateVec;
  
  InvSSABasicTikhonov( IceGrid  &grid, IceModelVec2V &f, PetscReal p, PetscReal q, InvTaucParameterization &tp);
  ~InvSSABasicTikhonov();

  PetscErrorCode setZeta(IceModelVec2S &zeta);

  void setEpsilon( PetscReal eps_vel, PetscReal eps_Du, PetscReal eps_nuH) {
    m_epsilon_velocity = eps_vel;
    m_epsilon_strainrate = eps_Du;
    m_epsilon_nuH = eps_nuH;
  }

  void setParamsFromConfig( NCConfigVariable &config );

  void setPseudoPlasticThreshold( PetscReal threshold) {
    m_pseudo_plastic_threshold = threshold;
  }


  void setBH(PetscReal B, PetscReal H) {
    m_B = B; m_H = H;
  }

  void setDirichletData(IceModelVec2Int &locs, IceModelVec2V &values) {
    m_dirichletLocations = &locs;
    m_dirichletValues = &values;
  }

  void setFixedDesignLocations(IceModelVec2Int &locs) {
    m_fixedDesignLocations = &locs;
  }

  void setObservationWeights(IceModelVec2S &weights) {
    m_observationWeights = &weights;
  }

  PetscErrorCode solve(bool &success);

  IceModelVec2V &solution() {
    return m_u;
  }

  SNESConvergedReason reason() {
    return m_reason;
  }
  
  std::string reasonDescription() {
    return std::string(SNESConvergedReasons[m_reason]);
  }

  PetscErrorCode assembleFunction( DMDALocalInfo *info, PISMVector2 **x, PISMVector2 **f);
  PetscErrorCode assembleJacobian( DMDALocalInfo *info, PISMVector2 **x, Mat J);

  PetscErrorCode linearizeAt( IceModelVec2S &c, bool &success);

  PetscErrorCode evalObjective(IceModelVec2S &dc, PetscReal *OUTPUT);
  PetscErrorCode evalGradObjective(IceModelVec2S &dc, IceModelVec2S &gradient);

  PetscErrorCode evalPenalty(IceModelVec2V &du, PetscReal *OUTPUT);
  PetscErrorCode evalGradPenaltyReduced(IceModelVec2V &dr, IceModelVec2S &gradient);

protected:

  PetscErrorCode construct();
  PetscErrorCode destruct();

  void computeNuH(PetscReal B, PetscReal H, const PetscReal Du[],
    PetscReal *nuH, PetscReal *dNuH );

  void computeBeta(PetscReal tauc, PISMVector2 U,
    PetscReal *beta, PetscReal *dBeta );

  IceGrid &m_grid;
  IceModelVec2S   *m_zeta;
  IceModelVec2V   *m_f;
  PetscReal m_p, m_q;
  IceModelVec2S   m_c;
  PetscReal m_B, m_H;

  PetscReal m_epsilon_velocity;
  PetscReal m_epsilon_strainrate;
  PetscReal m_epsilon_nuH;
  PetscReal m_pseudo_plastic_threshold;

  IceModelVec2Int *m_dirichletLocations;
  IceModelVec2V   *m_dirichletValues;
  PetscReal        m_dirichletWeight;

  IceModelVec2Int *m_fixedDesignLocations;
  IceModelVec2S   *m_observationWeights;

  InvTaucParameterization &m_tauc_param;
  
  IceModelVec2V  m_uGlobal;
  IceModelVec2V  m_u;
  IceModelVec2V  m_r;

  IceModelVec2V  m_vGlobal;
  IceModelVec2V  m_v;

  IceModelVec2V  m_adjointRHS;
  
  SNESDMCallbacks<InvSSABasicTikhonov,PISMVector2 **> m_callbacks;

  FEElementMap m_element_index;
  FEQuadrature m_quadrature;
  FEDOFMap     m_dofmap;

  std::auto_ptr< Functional<IceModelVec2S> > m_designFunctional;
  std::auto_ptr< Functional<IceModelVec2V> > m_penaltyFunctional;
  
  SNES m_snes;
  DM   m_da;
  KSP  m_ksp;
  Mat  m_J;
  
  SNESConvergedReason m_reason;

};

typedef TikhonovProblem<InvSSABasicTikhonov> InvSSABasicTikhonovProblem;
typedef TaoBasicSolver<InvSSABasicTikhonovProblem> InvSSABasicTikhonovSolver;


// SWIG is having problem recognising subclasses of templated superclasses
// have the correct type.  The following are workarounds to allow for
// python-based listeners to InvSSATikhonov

class InvSSABasicPythonListenerBridge: public TikhonovProblemListener<InvSSABasicTikhonov> {
public:
  InvSSABasicPythonListenerBridge(PythonTikhonovSVListener::Ptr core) : m_core(core) { }
  PetscErrorCode iteration( TikhonovProblem<InvSSABasicTikhonov> &,
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

void InvSSABasicTikhonovAddListener(TikhonovProblem<InvSSABasicTikhonov> &problem, 
                 PythonTikhonovSVListener::Ptr listener ) {
  std::tr1::shared_ptr<InvSSABasicPythonListenerBridge> bridge(new InvSSABasicPythonListenerBridge(listener) );
  problem.addListener(bridge);
}

#endif /* end of include guard: INVSSABASICTIKHONOV_HH_5AJEULAC */

