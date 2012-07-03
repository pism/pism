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

#include "InvSSATikhonovGN.hh"
#include <assert.h>

InvSSATikhonovGN::InvSSATikhonovGN( InvSSAForwardProblem &ssaforward,
DesignVec &d0, StateVec &u_obs, PetscReal eta,
IPFunctional<DesignVec> &designFunctional, IPFunctional<StateVec> &stateFunctional):
m_ssaforward(ssaforward), m_d0(d0), m_u_obs(u_obs), m_eta(eta),
m_designFunctional(designFunctional), m_stateFunctional(stateFunctional)
{
  PetscErrorCode ierr;
  ierr = this->construct();
  assert(ierr==0);
}

InvSSATikhonovGN::~InvSSATikhonovGN() {
  PetscErrorCode ierr;
  ierr = this->destruct(); CHKERRCONTINUE(ierr);
  assert(ierr==0);
}


PetscErrorCode InvSSATikhonovGN::construct() {
  PetscErrorCode ierr;
  IceGrid &grid = *m_d0.get_grid();
  m_comm = grid.com;

  PetscInt design_stencil_width = m_d0.get_stencil_width();
  PetscInt state_stencil_width = m_u_obs.get_stencil_width();

  ierr = m_x.create(grid,"x",kNoGhosts,0); CHKERRQ(ierr);
  ierr = m_y.create(grid,"y",kNoGhosts,0); CHKERRQ(ierr);

  ierr = m_tmp_D1.create(grid,"work vector",kNoGhosts,0); CHKERRQ(ierr);
  ierr = m_tmp_S1.create(grid,"work vector",kNoGhosts,0); CHKERRQ(ierr);
  ierr = m_tmp_S2.create(grid,"work vector",kNoGhosts,0); CHKERRQ(ierr);

  ierr = m_GN_rhs.create(grid,"GN_rhs",kNoGhosts,0); CHKERRQ(ierr);

  ierr = m_dGlobal.create(grid,"d (sans ghosts)",kNoGhosts,0); CHKERRQ(ierr);
  ierr = m_d.create(grid,"d",kHasGhosts,design_stencil_width); CHKERRQ(ierr);
  ierr = m_d_diff.create(grid,"d_diff",kHasGhosts,design_stencil_width); CHKERRQ(ierr);
  
  // ierr = m_uGlobal.create(grid,"u (sans ghosts)",kNoGhosts,0); CHKERRQ(ierr);
  ierr = m_u_diff.create(grid,"du",kHasGhosts,state_stencil_width); CHKERRQ(ierr);

  ierr = KSPCreate(grid.com, &m_ksp); CHKERRQ(ierr);
  ierr = KSPSetOptionsPrefix(m_ksp,"inv_gn_"); CHKERRQ(ierr);
  PetscReal ksp_rtol = 1e-2; // Soft tolerance
  ierr = KSPSetTolerances(m_ksp,ksp_rtol,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);
  ierr = KSPSetType(m_ksp,KSPCG); CHKERRQ(ierr);
  PC pc;
  ierr = KSPGetPC(m_ksp,&pc); CHKERRQ(ierr);
  ierr = PCSetType(pc,PCNONE); CHKERRQ(ierr);
  ierr = KSPSetFromOptions(m_ksp); CHKERRQ(ierr);  

  PetscInt nLocalNodes  = grid.xm*grid.ym;
  PetscInt nGlobalNodes = grid.Mx*grid.My;
  ierr = MatCreateShell(grid.com,nLocalNodes,nLocalNodes,nGlobalNodes,nGlobalNodes,this,&m_mat_GN); CHKERRQ(ierr);

  typedef MatrixMultiplyCallback<InvSSATikhonovGN,&InvSSATikhonovGN::apply_GN> multCallback;
  ierr = multCallback::connect(m_mat_GN);

  return 0;
}

PetscErrorCode InvSSATikhonovGN::destruct() {
  PetscErrorCode ierr;
  ierr = KSPDestroy(&m_ksp); CHKERRQ(ierr);
  ierr = MatDestroy(&m_mat_GN);
  return 0;
}

PetscErrorCode InvSSATikhonovGN::init() {
  PetscErrorCode ierr;
  bool success;
  ierr = m_ssaforward.linearize_at(m_d0,success); CHKERRQ(ierr);
  return 0;
}

PetscErrorCode InvSSATikhonovGN::apply_GN(IceModelVec2S &x,IceModelVec2S &y) {
  PetscErrorCode ierr;
  ierr = this->apply_GN(x.get_vec(),y.get_vec()); CHKERRQ(ierr);
  return 0; 
}

PetscErrorCode InvSSATikhonovGN::apply_GN(Vec x, Vec y) {
  PetscErrorCode ierr;
  
  // FIXME: Needless copies for now.
  ierr = m_x.copy_from(x); CHKERRQ(ierr);
  
  
  ierr = m_ssaforward.apply_linearization(m_x,m_tmp_S1); CHKERRQ(ierr);
  ierr = m_stateFunctional.interior_product(m_tmp_S1,m_tmp_S2); CHKERRQ(ierr);
  ierr = m_ssaforward.apply_linearization_transpose(m_tmp_S2,m_y); CHKERRQ(ierr);
  // ierr = verbPrintf(1,m_comm,"part 1:\n");
  // ierr = VecView(m_y.get_vec(),PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);


  ierr = m_designFunctional.interior_product(m_x,m_tmp_D1); CHKERRQ(ierr);
  // ierr = verbPrintf(1,m_comm,"part 2:\n");
  // ierr = VecView(m_tmp_D1.get_vec(),PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
  ierr = m_y.add(1./m_eta,m_tmp_D1);

  // ierr = verbPrintf(1,m_comm,"combined:\n");
  // ierr = VecView(m_y.get_vec(),PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
  
  ierr = m_y.copy_to(y); CHKERRQ(ierr);
  return 0;
}  
  
PetscErrorCode InvSSATikhonovGN::assemble_GN_rhs(DesignVec &rhs) {
  PetscErrorCode ierr;
  
  ierr = m_stateFunctional.interior_product(m_u_diff,m_tmp_S1); CHKERRQ(ierr);
  ierr = m_ssaforward.apply_linearization_transpose(m_tmp_S1,rhs); CHKERRQ(ierr);

  ierr = m_designFunctional.interior_product(m_d_diff,m_tmp_D1); CHKERRQ(ierr);
  ierr = rhs.add(1./m_eta,m_tmp_D1);

  return 0;
}

PetscErrorCode InvSSATikhonovGN::solve() {
  PetscErrorCode ierr;

  ierr = m_d_diff.copy_from(m_d); CHKERRQ(ierr);
  ierr = m_d_diff.add(-1,m_d0); CHKERRQ(ierr);
  
  bool success;
  ierr = m_ssaforward.linearize_at(m_d0,success); CHKERRQ(ierr);
  ierr = m_u_diff.copy_from(m_ssaforward.solution()); CHKERRQ(ierr);
  ierr = m_u_diff.add(-1,m_u_obs); CHKERRQ(ierr);
  
  ierr = this->assemble_GN_rhs(m_GN_rhs); CHKERRQ(ierr);
  
  ierr = KSPSetOperators(m_ksp,m_mat_GN,m_mat_GN,SAME_NONZERO_PATTERN); CHKERRQ(ierr);
  ierr = KSPSolve(m_ksp,m_GN_rhs.get_vec(),m_dGlobal.get_vec()); CHKERRQ(ierr);

  KSPConvergedReason reason;
  ierr = KSPGetConvergedReason(m_ksp,&reason); CHKERRQ(ierr);
  if(reason<0) {
    SETERRQ1(m_comm,1,"InvSSATikhonovGN::solve KSP failed %d\n",reason);    
  }

  ierr = m_d.copy_from(m_dGlobal); CHKERRQ(ierr);

  return 0;
}
