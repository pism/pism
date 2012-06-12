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

#include "InvSSATikhonov.hh"
#include "Mask.hh"
#include "basal_resistance.hh"
#include "PISMVars.hh"
#include <assert.h>
#include "H1NormFunctional.hh"
#include "MeanSquareFunctional.hh"
#include "pism_options.hh"

InvSSATikhonov::InvSSATikhonov(IceGrid &g, IceBasalResistancePlasticLaw &b,
  EnthalpyConverter &e, InvTaucParameterization &tp,
  const NCConfigVariable &c) : SSAFEM(g,b,e,c),
  m_grid(grid), m_zeta(NULL), 
  m_fixed_tauc_locations(NULL), m_misfit_weight(NULL), 
  m_tauc_param(tp), m_element_index(grid) 
{
  PetscErrorCode ierr = this->construct();
  CHKERRCONTINUE(ierr);
  assert(ierr == 0);
}

InvSSATikhonov::~InvSSATikhonov() {
  PetscErrorCode ierr = this->destruct();
  CHKERRCONTINUE(ierr);
  assert(ierr == 0);
}

// FIXME: replace this
PetscErrorCode InvSSATikhonov::set_functionals() {

  PetscReal cL2 = m_grid.config.get("inv_ssa_cL2");
  PetscReal cH1 = m_grid.config.get("inv_ssa_cH1");

  m_designFunctional.reset(new H1NormFunctional2S(m_grid,cL2,cH1,m_fixed_tauc_locations));    
  // PetscReal stress_scale = m_grid.config.get("tauc_param_tauc_scale");
  // m_designFunctional.reset(new MeanSquareFunctional2S(m_grid));
  // (reinterpret_cast<MeanSquareFunctional2S&>(*m_designFunctional)).normalize(stress_scale);

  PetscReal velocity_scale = m_grid.config.get("inv_ssa_velocity_scale")/secpera;
  m_penaltyFunctional.reset(new MeanSquareFunctional2V(m_grid,m_misfit_weight));    
  (reinterpret_cast<MeanSquareFunctional2V&>(*m_penaltyFunctional)).normalize(velocity_scale);

  m_domainIP.reset(new MeanSquareFunctional2S(m_grid));

  return 0;
}


// Initialize the solver, called once by the client before use.
PetscErrorCode InvSSATikhonov::init(PISMVars &vars) {
  PetscErrorCode ierr;

  ierr = SSAFEM::init(vars); CHKERRQ(ierr);

  m_misfit_weight = dynamic_cast<IceModelVec2S*>(vars.get("vel_misfit_weight"));
  if (m_misfit_weight == NULL){
    verbPrintf(3,grid.com,"Weight for inverse problem L2 norm not available; using standard L2 norm.\n");
  }

  return 0;
}

// m_uGlobal "=" SSAX
// m_u = velocity

PetscErrorCode InvSSATikhonov::construct() {
  PetscErrorCode ierr;
  PetscInt stencilWidth = 1;

  m_vGlobal.create(m_grid,"adjoint work vector (sans ghosts)",kNoGhosts,stencilWidth);
  m_v.create(m_grid,"adjoint work vector",kHasGhosts,stencilWidth);
  m_adjointRHS.create(m_grid,"adjoint RHS",kNoGhosts,stencilWidth);

  ierr = DMGetMatrix(SSADA, "baij", &m_Jadjoint); CHKERRQ(ierr);

  ierr = KSPCreate(m_grid.com, &m_ksp); CHKERRQ(ierr);
  PetscReal ksp_rtol = 1e-12;
  ierr = KSPSetTolerances(m_ksp,ksp_rtol,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);
  PC pc;
  ierr = KSPGetPC(m_ksp,&pc); CHKERRQ(ierr);
  ierr = PCSetType(pc,PCBJACOBI); CHKERRQ(ierr);
  ierr = KSPSetFromOptions(m_ksp); CHKERRQ(ierr);  

  m_quadrature.init(m_grid);
  return 0;
}

PetscErrorCode InvSSATikhonov::destruct() {
  PetscErrorCode ierr;
  ierr = MatDestroy(&m_Jadjoint); CHKERRQ(ierr);
  ierr = KSPDestroy(&m_ksp); CHKERRQ(ierr);
  return 0;
}

PetscErrorCode InvSSATikhonov::set_zeta(IceModelVec2S &new_zeta )
{
  PetscErrorCode ierr;
  PetscInt i,j,q;

  IceModelVec2S *m_tauc = tauc;

  m_zeta = &new_zeta;
  
  // Convert zeta to tauc.
  m_tauc_param.convertToTauc(*m_zeta,*m_tauc);

  // Cache tauc at the quadrature points in feStore.
  PetscReal **tauc_a;
  PetscReal tauc_q[FEQuadrature::Nq];
  ierr = m_tauc->get_array(tauc_a); CHKERRQ(ierr);
  PetscInt xs = m_element_index.xs, xm = m_element_index.xm,
           ys = m_element_index.ys, ym = m_element_index.ym;
  for (i=xs; i<xs+xm; i++) {
    for (j=ys;j<ys+ym; j++) {
      quadrature.computeTrialFunctionValues(i,j,dofmap,tauc_a,tauc_q);
      const PetscInt ij = m_element_index.flatten(i,j);
      FEStoreNode *feS = &feStore[ij*FEQuadrature::Nq];
      for (q=0; q<4; q++) {
        feS[q].tauc = tauc_q[q];
      }
    }
  }
  ierr = m_tauc->end_access(); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode InvSSATikhonov::linearizeAt( IceModelVec2S &zeta, bool &success) {

  PetscErrorCode ierr;
  ierr = this->set_zeta(zeta); CHKERRQ(ierr);

  ierr = this->solve(); CHKERRQ(ierr);

  success = true;
  return 0;
}

PetscErrorCode InvSSATikhonov::assembleFunction(IceModelVec2V u, Vec RHS) {
  PetscErrorCode ierr;

  PISMVector2 **u_a, **rhs_a;

  ierr = u.get_array(u_a); CHKERRQ(ierr);
  ierr = DMDAVecGetArray(SSADA,RHS,&rhs_a); CHKERRQ(ierr);    

  DMDALocalInfo *info = NULL;
  ierr = this->compute_local_function(info,const_cast<const PISMVector2 **>(u_a),rhs_a); CHKERRQ(ierr);

  ierr = DMDAVecRestoreArray(SSADA,RHS,&rhs_a); CHKERRQ(ierr);
  ierr = u.end_access(); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode InvSSATikhonov::assembleJacobian(IceModelVec2V u, Mat Jac) {
  PetscErrorCode ierr;

  PISMVector2 **u_a;
  ierr = u.get_array(u_a); CHKERRQ(ierr);

  DMDALocalInfo *info = NULL;
  ierr = this->compute_local_jacobian(info,const_cast<const PISMVector2 **>(u_a),Jac); CHKERRQ(ierr);
  
  ierr = u.end_access(); CHKERRQ(ierr);
  
  return 0;
}

PetscErrorCode InvSSATikhonov::evalObjective(IceModelVec2S &dzeta, PetscReal *OUTPUT) {
  PetscErrorCode ierr;
  ierr = m_designFunctional->valueAt(dzeta,OUTPUT);
  return 0;
}

PetscErrorCode InvSSATikhonov::evalGradObjective(IceModelVec2S &dzeta, IceModelVec2S &gradient) {
  PetscErrorCode ierr;
  ierr = m_designFunctional->gradientAt(dzeta,gradient); CHKERRQ(ierr);
  return 0;
}

PetscErrorCode InvSSATikhonov::evalPenalty(IceModelVec2V &du, PetscReal *OUTPUT) {
  PetscErrorCode ierr;
  ierr = m_penaltyFunctional->valueAt(du,OUTPUT);
  return 0;
}

PetscErrorCode InvSSATikhonov::evalGradPenalty(IceModelVec2V &du, IceModelVec2V &gradient) {
  PetscErrorCode ierr;
  ierr = m_penaltyFunctional->gradientAt(du,gradient); CHKERRQ(ierr);
  return 0;
}

PetscErrorCode InvSSATikhonov::evalGradPenaltyReducedFD(IceModelVec2V &du, IceModelVec2S &gradient) {
  PetscErrorCode ierr;
  bool success;
  PetscReal h = PETSC_SQRT_MACHINE_EPSILON;

  IceModelVec2S &zeta = *m_zeta;

  IceModelVec2V u0;
  ierr = u0.create(grid,"u",kHasGhosts,1); CHKERRQ(ierr);
  ierr = u0.copy_from(velocity); CHKERRQ(ierr);

  IceModelVec2V uprime;
  ierr = uprime.create(grid,"uprime",kHasGhosts,1); CHKERRQ(ierr);

  IceModelVec2V uprime2;
  ierr = uprime2.create(grid,"uprime2",kHasGhosts,1); CHKERRQ(ierr);

  IceModelVec2S dzeta;
  ierr = dzeta.create(grid,"dzeta",kHasGhosts,1); CHKERRQ(ierr);
  dzeta.set(0);

  IceModelVec2V gradPenalty;
  ierr = gradPenalty.create(grid,"gradPenalty",kNoGhosts,0); CHKERRQ(ierr);
  ierr = this->evalGradPenalty(du,gradPenalty); CHKERRQ(ierr);
  
  ierr = gradient.begin_access(); CHKERRQ(ierr);
  for(PetscInt i=grid.xs; i< grid.xs+grid.xm; i++) {
    for(PetscInt j=grid.ys; j< grid.ys+grid.ym; j++) {
      ierr = zeta.begin_access(); CHKERRQ(ierr);
      zeta(i,j) += h;
      ierr = zeta.end_access(); CHKERRQ(ierr);
      ierr = zeta.beginGhostComm(); CHKERRQ(ierr);
      ierr = zeta.endGhostComm(); CHKERRQ(ierr);
      
      ierr = this->linearizeAt(zeta,success); CHKERRQ(ierr);

      ierr = zeta.begin_access(); CHKERRQ(ierr);
      zeta(i,j) -= h;
      ierr = zeta.end_access(); CHKERRQ(ierr);
      ierr = zeta.beginGhostComm(); CHKERRQ(ierr);
      ierr = zeta.endGhostComm(); CHKERRQ(ierr);

      ierr = uprime.copy_from(velocity); CHKERRQ(ierr);
      ierr = uprime.add(-1,u0); CHKERRQ(ierr);
      ierr = uprime.scale(1/h); CHKERRQ(ierr);

      ierr = dzeta.begin_access(); CHKERRQ(ierr);
      dzeta(i,j) = h;
      ierr = dzeta.end_access(); CHKERRQ(ierr);
      ierr = dzeta.beginGhostComm(); CHKERRQ(ierr);
      ierr = dzeta.endGhostComm(); CHKERRQ(ierr);

      ierr = this->computeT(dzeta,uprime2);

      ierr = dzeta.begin_access(); CHKERRQ(ierr);
      dzeta(i,j) = 0;
      ierr = dzeta.end_access(); CHKERRQ(ierr);
      ierr = dzeta.beginGhostComm(); CHKERRQ(ierr);
      ierr = dzeta.endGhostComm(); CHKERRQ(ierr);

      PetscReal g=0;
      ierr = uprime.begin_access(); CHKERRQ(ierr);
      ierr = gradPenalty.begin_access(); CHKERRQ(ierr);
      for(PetscInt k=grid.xs; k< grid.xs+grid.xm; k++) {
        for(PetscInt l=grid.ys; l< grid.ys+grid.ym; l++) {
          g += gradPenalty(k,l).u*uprime(k,l).u+gradPenalty(k,l).v*uprime(k,l).v;
        }
      }
      ierr = uprime.end_access(); CHKERRQ(ierr);
      ierr = gradPenalty.end_access(); CHKERRQ(ierr);
      gradient(i,j) = g;
    }
  }
  ierr = gradient.end_access(); CHKERRQ(ierr);
  ierr = linearizeAt(zeta, success); CHKERRQ(ierr);
  return 0;
}



PetscErrorCode InvSSATikhonov::evalGradPenaltyReducedNoTranspose(IceModelVec2V &du, IceModelVec2S &gradient) {
  PetscErrorCode ierr;
  bool success;
  
  printf("NO TRANSPOSE!\n");

  IceModelVec2S &zeta = *m_zeta;

  IceModelVec2V u0;
  ierr = u0.create(grid,"u",kHasGhosts,1); CHKERRQ(ierr);
  ierr = u0.copy_from(velocity); CHKERRQ(ierr);

  IceModelVec2V uprime;
  ierr = uprime.create(grid,"uprime",kHasGhosts,1); CHKERRQ(ierr);

  IceModelVec2S dzeta;
  ierr = dzeta.create(grid,"dzeta",kHasGhosts,1); CHKERRQ(ierr);
  dzeta.set(0);

  IceModelVec2V gradPenalty;
  ierr = gradPenalty.create(grid,"gradPenalty",kNoGhosts,0); CHKERRQ(ierr);
  ierr = this->evalGradPenalty(du,gradPenalty); CHKERRQ(ierr);
  
  ierr = gradient.begin_access(); CHKERRQ(ierr);
  for(PetscInt i=grid.xs; i< grid.xs+grid.xm; i++) {
    for(PetscInt j=grid.ys; j< grid.ys+grid.ym; j++) {
      ierr = dzeta.begin_access(); CHKERRQ(ierr);
      dzeta(i,j) = 1;
      ierr = dzeta.end_access(); CHKERRQ(ierr);
      ierr = dzeta.beginGhostComm(); CHKERRQ(ierr);
      ierr = dzeta.endGhostComm(); CHKERRQ(ierr);

      ierr = this->computeT(dzeta,uprime);

      ierr = dzeta.begin_access(); CHKERRQ(ierr);
      dzeta(i,j) = 0;
      ierr = dzeta.end_access(); CHKERRQ(ierr);
      ierr = dzeta.beginGhostComm(); CHKERRQ(ierr);
      ierr = dzeta.endGhostComm(); CHKERRQ(ierr);

      PetscReal g=0;
      ierr = uprime.begin_access(); CHKERRQ(ierr);
      ierr = gradPenalty.begin_access(); CHKERRQ(ierr);
      for(PetscInt k=grid.xs; k< grid.xs+grid.xm; k++) {
        for(PetscInt l=grid.ys; l< grid.ys+grid.ym; l++) {
          g += gradPenalty(k,l).u*uprime(k,l).u+gradPenalty(k,l).v*uprime(k,l).v;
        }
      }
      ierr = uprime.end_access(); CHKERRQ(ierr);
      ierr = gradPenalty.end_access(); CHKERRQ(ierr);
      gradient(i,j) = g;
    }
  }
  ierr = gradient.end_access(); CHKERRQ(ierr);
  ierr = linearizeAt(zeta, success); CHKERRQ(ierr);
  return 0;
}

PetscErrorCode InvSSATikhonov::computeT(IceModelVec2S &dzeta,IceModelVec2V &du) {
  PetscErrorCode ierr;
  
  // Assemble the Jacobian matrix.
  DMDALocalInfo *info = NULL;
  PISMVector2 **u_a;

  ierr = velocity.get_array(u_a); CHKERRQ(ierr);
  this->compute_local_jacobian( info, const_cast<const PISMVector2**>(u_a), m_Jadjoint);
  ierr = velocity.end_access(); CHKERRQ(ierr);

  ierr = this->assemble_T_rhs(*m_zeta,velocity,dzeta,m_adjointRHS);

  KSPConvergedReason  kspreason;
  // call PETSc to solve linear system by iterative method.
  ierr = KSPSetOperators(m_ksp, m_Jadjoint, m_Jadjoint, SAME_NONZERO_PATTERN); CHKERRQ(ierr);
  ierr = KSPSolve(m_ksp, m_adjointRHS.get_vec(), m_vGlobal.get_vec()); CHKERRQ(ierr); // SOLVE

  ierr = KSPGetConvergedReason(m_ksp, &kspreason); CHKERRQ(ierr);
  
  if (kspreason < 0) {
    SETERRQ1(m_grid.com,1,"InvSSATikhonov adjoint linear solve failed (KSP reason %s)",KSPConvergedReasons[kspreason]);
  }

  ierr = du.copy_from(m_vGlobal); CHKERRQ(ierr);
  ierr = du.beginGhostComm(); CHKERRQ(ierr);
  ierr = du.endGhostComm(); CHKERRQ(ierr);
  return 0;
}

PetscErrorCode InvSSATikhonov::assemble_T_rhs(IceModelVec2S &zeta, IceModelVec2V &u, IceModelVec2S &dzeta, IceModelVec2V &rhs) {
  PetscErrorCode ierr;

  PetscReal **zeta_a;
  ierr = zeta.get_array(zeta_a); CHKERRQ(ierr);
  
  PISMVector2 **u_a;
  ierr = u.get_array(u_a); CHKERRQ(ierr);

  PetscReal **dzeta_a;
  ierr = dzeta.get_array(dzeta_a); CHKERRQ(ierr);

  PISMVector2 **rhs_a;
  ierr = rhs.get_array(rhs_a); CHKERRQ(ierr);

  ierr = this->assemble_T_rhs(zeta_a,u_a,dzeta_a,rhs_a);

  ierr = zeta.end_access(); CHKERRQ(ierr);
  ierr = u.end_access(); CHKERRQ(ierr);
  ierr = dzeta.end_access(); CHKERRQ(ierr);
  ierr = rhs.end_access(); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode InvSSATikhonov::assemble_T_rhs(PetscReal **zeta_a, PISMVector2 **u_a, PetscReal **dzeta_a, PISMVector2 **rhs_a) {
  PetscInt         i,j;
  PetscErrorCode ierr;

  // Zero out the portion of the function we are responsible for computing.
  for (i=grid.xs; i<grid.xs+grid.xm; i++) {
    for (j=grid.ys; j<grid.ys+grid.ym; j++) {
      rhs_a[i][j].u = rhs_a[i][j].v = 0;
    }
  }

  // Aliases to help with notation consistency below.
  IceModelVec2Int *m_dirichletLocations = bc_locations;
  IceModelVec2V   *m_dirichletValues    = vel_bc;
  PetscReal        m_dirichletWeight    = dirichletScale;

  PISMVector2 u_e[FEQuadrature::Nk];
  PISMVector2 u_q[FEQuadrature::Nq];

  PISMVector2 rhs_e[FEQuadrature::Nk];

  PetscReal dzeta_e[FEQuadrature::Nk];

  PetscReal zeta_e[FEQuadrature::Nk];

  PetscReal dtauc_e[FEQuadrature::Nk];
  PetscReal dtauc_q[FEQuadrature::Nq];

  // An Nq by Nk array of test function values.
  const FEFunctionGerm (*test)[FEQuadrature::Nk] = m_quadrature.testFunctionValues();

  DirichletData dirichletBC;
  ierr = dirichletBC.init(m_dirichletLocations,m_dirichletValues,m_dirichletWeight); CHKERRQ(ierr);
  DirichletData fixedZeta;
  ierr = fixedZeta.init(m_fixed_tauc_locations);

  // Jacobian times weights for quadrature.
  PetscScalar JxW[FEQuadrature::Nq];
  m_quadrature.getWeightedJacobian(JxW);

  // Mask (query?) for determining where ice is grounded.
  Mask M;

  // Loop through all elements.
  PetscInt xs = m_element_index.xs, xm = m_element_index.xm,
           ys = m_element_index.ys, ym = m_element_index.ym;
  for (i=xs; i<xs+xm; i++) {
    for (j=ys; j<ys+ym; j++) {

      // Zero out the element-local residual in prep for updating it.
      for(PetscInt k=0;k<FEQuadrature::Nk;k++){
        rhs_e[k].u = 0;
        rhs_e[k].v = 0;
      }

      // Index into coefficient storage in feStore
      const PetscInt ij = element_index.flatten(i,j);

      // Initialize the map from global to local degrees of freedom for this element.
      m_dofmap.reset(i,j,m_grid);
      
      // Obtain the value of the solution at the nodes adjacent to the element,
      // fix dirichlet values, and compute values at quad pts.
      m_dofmap.extractLocalDOFs(i,j,u_a,u_e);
      if(dirichletBC) dirichletBC.update(m_dofmap,u_e);
      m_quadrature.computeTrialFunctionValues(u_e,u_q);

      // Compute dzeta at the nodes
      m_dofmap.extractLocalDOFs(i,j,dzeta_a,dzeta_e);

      // Compute the change in tau_c with respect to zeta at the quad points.
      m_dofmap.extractLocalDOFs(i,j,zeta_a,zeta_e);
      for(PetscInt k=0;k<FEQuadrature::Nk;k++){
        m_tauc_param.toTauc(zeta_e[k],NULL,dtauc_e + k);
        dtauc_e[k]*=dzeta_e[k];
      }
      if(fixedZeta) fixedZeta.updateHomogeneous(m_dofmap,dtauc_e);
      m_quadrature.computeTrialFunctionValues(dtauc_e,dtauc_q);

      for (PetscInt q=0; q<FEQuadrature::Nq; q++) {
        PISMVector2 u_qq = u_q[q];

        const FEStoreNode *feS = &feStore[ij*FEQuadrature::Nq+q];

        // Determine "dbeta/dtauc" at the quadrature point
        PetscReal dbeta = 0;
        if( M.grounded_ice(feS->mask) ) {
          dbeta = basal.drag(dtauc_q[q],u_qq.u,u_qq.v);
        }

        for (PetscInt k=0; k<FEQuadrature::Nk; k++) {
          rhs_e[k].u += -JxW[q]*dbeta*u_qq.u*test[q][k].val;
          rhs_e[k].v += -JxW[q]*dbeta*u_qq.v*test[q][k].val;
        }
      } // q
      m_dofmap.addLocalResidualBlock(rhs_e,rhs_a);
    } // j
  } // i

  if(dirichletBC) dirichletBC.fixResidualHomogeneous(rhs_a);
  ierr = dirichletBC.finish(); CHKERRQ(ierr);

  return 0;
}

/*
PetscErrorCode InvSSATikhonov::computeJacobianDesign(IceModelVec2S &zeta, IceModelVec2V &u, Mat Jac) {

  PetscInt         i,j;
  PetscErrorCode   ierr;

  IceModelVec2Int *m_dirichletLocations = bc_locations;
  IceModelVec2V   *m_dirichletValues    = vel_bc;
  PetscReal        m_dirichletWeight    = dirichletScale;

  // Zero out the Jacobian in preparation for updating it.
  ierr = MatZeroEntries(Jac);CHKERRQ(ierr);

  // Jacobian times weights for quadrature.
  PetscScalar JxW[FEQuadrature::Nq];
  quadrature.getWeightedJacobian(JxW);

  PISMVector2 **u_a;
  PISMVector2  u_e[FEQuadrature::Nk];
  PISMVector2  u_q[FEQuadrature::Nq];
  ierr = u.get_array(u_a); CHKERRQ(ierr);

  PetscScalar **zeta_a;
  PetscScalar zeta_e[FEQuadrature::Nk];
  zeta.get_array(zeta_a);

  PetscScalar dtauc_dzeta_e[FEQuadrature::Nk];

  DirichletData dirichletBC;
  ierr = dirichletBC.init(m_dirichletLocations,m_dirichletValues,m_dirichletWeight); CHKERRQ(ierr);
  DirichletData fixedZeta;
  ierr = fixedZeta.init(m_fixed_tauc_locations);

  // Values of the finite element test functions at the quadrature points.
  // This is an Nq by Nk array of function germs (Nq=#of quad pts, Nk=#of test functions).
  const FEFunctionGerm (*test)[FEQuadrature::Nk] = quadrature.testFunctionValues();

  // Loop through all the elements.
  PetscInt xs = element_index.xs, xm = element_index.xm,
           ys = element_index.ys, ym = element_index.ym;
  for (i=xs; i<xs+xm; i++) {
    for (j=ys; j<ys+ym; j++) {

      // Element-local Jacobian matrix.
      PetscReal      K[(2*FEQuadrature::Nk)*(FEQuadrature::Nk)];

      // Index into the coefficient storage array.
      const PetscInt ij = element_index.flatten(i,j);

      // Initialize the map from global to local degrees of freedom for this element.
      dofmap.reset(i,j,grid);

      // Obtain the value of the solution at the adjacent nodes to the element.
      dofmap.extractLocalDOFs(i,j,u_a,u_e);
      if(dirichletData) dirichletData.update(u_e);
      quadrature.computeTrialFunctionValues(u_e,u_q);

      // Compute the change in tau_c with respect to zeta at the element nodes.
      m_dofmap.extractLocalDOFs(i,j,zeta_a,zeta_e);
      for(PetscInt k=0;k<FEQuadrature::Nk;k++){
        m_tauc_param.toTauc(zeta_e[k],NULL,dtauc_dzeta_e + k);
      }
      if(fixedZeta) fixedZeta.updateHomogeneous(m_dofmap,dtauc_dzeta_e);

      // Build the element-local Jacobian.
      ierr = PetscMemzero(K,sizeof(K));CHKERRQ(ierr);
      for (PetscInt q=0; q<FEQuadrature::Nq; q++) {

        // Shorthand for values and derivatives of the solution at the single quadrature point.
        PISMVector2 &u_qq = u_q[q];

        // Coefficients evaluated at the single quadrature point.
        const FEStoreNode *feS = &feStore[ij*Quadrature::Nq + q];
        const PetscReal    jw  = JxW[q];

        // Determine "dbeta/dtauc" at the quadrature point
        PetscReal dbeta_dtauc = 0;
        if( M.grounded_ice(feS->mask) ) {
          dbeta_dtauc = basal.drag(1.,u_qq.u,u_qq.v);
        }

        for (PetscInt k=0; k<FEQuadrature::Nk; k++) {
          for (PetscInt l=0; l<FEQuadrature::Nk; l++) {
            K[2*k*Quadrature::Nk +2*l  ] += JxW[q]*dbeta_dtauc*u_qq.u*test[q][k].val*test[q][l].val*dtauc_dzeta_e[l];
            K[2*k*Quadrature::Nk +2*l+1] += JxW[q]*dbeta_dtauc*u_qq.v*test[q][k].val*test[q][l].val*dtauc_dzeta_e[l];
          }
        }
      } // q

      

????      ierr = dofmap.addLocalJacobianBlock(K,Jac);

    } // j
  } // i

  ierr = MatAssemblyBegin(Jac,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(Jac,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

  ierr = MatSetOption(Jac,MAT_NEW_NONZERO_LOCATION_ERR,PETSC_TRUE);CHKERRQ(ierr);

  return 0;
}
*/

PetscErrorCode InvSSATikhonov::evalGradPenaltyReduced(IceModelVec2V &du, IceModelVec2S &gradient) {
  PetscErrorCode ierr;

  bool do_fd;
  ierr = PISMOptionsIsSet("-inv_ssa_tik_fd","compute tikhonov forward linearization via finite differences",do_fd); CHKERRQ(ierr);
  if(do_fd) {
    ierr = evalGradPenaltyReducedFD(du,gradient); CHKERRQ(ierr);
    return 0;
  }

  // Some aliases to help with notation consistency below.
  IceModelVec2V   &m_u                  = velocity;
  IceModelVec2Int *m_dirichletLocations = bc_locations;
  IceModelVec2V   *m_dirichletValues    = vel_bc;
  PetscReal        m_dirichletWeight    = dirichletScale;

  // Clear the gradient in prep for updating it.
  gradient.set(0);

  PISMVector2 **u_a;
  ierr = m_u.get_array(u_a); CHKERRQ(ierr);  

  // Assemble the Jacobian matrix.
  DMDALocalInfo *info = NULL;
  this->compute_local_jacobian( info, const_cast<const PISMVector2**>(u_a), m_Jadjoint);

  ierr = m_penaltyFunctional->gradientAt(du,m_adjointRHS); CHKERRQ(ierr);
  // ierr = gradientFD(*m_penaltyFunctional,du,m_adjointRHS); CHKERRQ(ierr);

  if(m_dirichletLocations) {
    DirichletData dirichletBC;
    ierr = dirichletBC.init(m_dirichletLocations,m_dirichletValues,m_dirichletWeight); CHKERRQ(ierr);
    PISMVector2 **rhs_a;
    ierr = m_adjointRHS.get_array(rhs_a); CHKERRQ(ierr);
    dirichletBC.fixResidualHomogeneous(rhs_a);
    ierr = dirichletBC.finish(); CHKERRQ(ierr); 
    ierr = m_adjointRHS.end_access(); CHKERRQ(ierr);
  }

  KSPConvergedReason  kspreason;
  // call PETSc to solve linear system by iterative method.
  ierr = KSPSetOperators(m_ksp, m_Jadjoint, m_Jadjoint, SAME_NONZERO_PATTERN); CHKERRQ(ierr);
  ierr = KSPSolve(m_ksp, m_adjointRHS.get_vec(), m_vGlobal.get_vec()); CHKERRQ(ierr); // SOLVE

  ierr = KSPGetConvergedReason(m_ksp, &kspreason); CHKERRQ(ierr);
  
  if (kspreason < 0) {
    SETERRQ1(m_grid.com,1,"InvSSATikhonov adjoint linear solve failed (KSP reason %s)",KSPConvergedReasons[kspreason]);
  }

  ierr = m_v.copy_from(m_vGlobal); CHKERRQ(ierr);

  PISMVector2 **v_a;
  ierr = m_v.get_array(v_a); CHKERRQ(ierr);

  PetscReal **gradient_a;
  ierr = gradient.get_array(gradient_a); CHKERRQ(ierr);
  
  PetscReal **zeta_a;
  ierr = m_zeta->get_array(zeta_a); CHKERRQ(ierr);
  
  ierr = this->compute_Jdesign_transpose(zeta_a,u_a,v_a,gradient_a); CHKERRQ(ierr);


  ierr = m_v.end_access(); CHKERRQ(ierr);
  ierr = m_u.end_access(); CHKERRQ(ierr);
  ierr = gradient.end_access(); CHKERRQ(ierr);
  ierr = m_zeta->end_access(); CHKERRQ(ierr);
  
  return 0;
}

PetscErrorCode InvSSATikhonov::compute_Jdesign_transpose(PetscReal **zeta_a, PISMVector2 **u_a, PISMVector2 **v_a, PetscReal **gradient_a) {
  PetscInt         i,j;
  PetscErrorCode ierr;

  PISMVector2 v_e[FEQuadrature::Nk];
  PISMVector2 v_q[FEQuadrature::Nq];

  PISMVector2 u_e[FEQuadrature::Nk];
  PISMVector2 u_q[FEQuadrature::Nq];

  PetscReal gradient_e[FEQuadrature::Nk];

  // An Nq by Nk array of test function values.
  const FEFunctionGerm (*test)[FEQuadrature::Nk] = m_quadrature.testFunctionValues();

  // Jacobian times weights for quadrature.
  PetscScalar JxW[FEQuadrature::Nq];
  m_quadrature.getWeightedJacobian(JxW);

  // Mask (query?) for determining where ice is grounded.
  Mask M;

  PetscInt xs = m_element_index.xs, xm = m_element_index.xm,
           ys = m_element_index.ys, ym = m_element_index.ym;
  for (i=xs; i<xs+xm; i++) {
    for (j=ys; j<ys+ym; j++) {

      // Index into coefficient storage in feStore
      const PetscInt ij = element_index.flatten(i,j);

      // Initialize the map from global to local degrees of freedom for this element.
      m_dofmap.reset(i,j,m_grid);
      
      // Obtain the value of the solution at the nodes adjacent to the element.
      // Compute the solution values and symmetric gradient at the quadrature points.
      m_dofmap.extractLocalDOFs(i,j,v_a,v_e);
      m_quadrature.computeTrialFunctionValues(v_e,v_q);
      m_dofmap.extractLocalDOFs(i,j,u_a,u_e);
      m_quadrature.computeTrialFunctionValues(u_e,u_q);

      // Zero out the element-local residual in prep for updating it.
      for(PetscInt k=0;k<FEQuadrature::Nk;k++){
        gradient_e[k] = 0;
      }

      for (PetscInt q=0; q<FEQuadrature::Nq; q++) {
        PISMVector2 v_qq = v_q[q];
        PISMVector2 u_qq = u_q[q];

        const FEStoreNode *feS = &feStore[ij*FEQuadrature::Nq+q];

        // Determine "dbeta/dtauc" at the quadrature point
        PetscReal dbeta_dtauc = 0;
        if( M.grounded_ice(feS->mask) ) {
          dbeta_dtauc = basal.drag(1.,u_qq.u,u_qq.v);
        }

        for (PetscInt k=0; k<FEQuadrature::Nk; k++) {
          gradient_e[k] += -JxW[q]*dbeta_dtauc*(v_qq.u*u_qq.u+v_qq.v*u_qq.v)*test[q][k].val;
        }
      } // q
      m_dofmap.addLocalResidualBlock(gradient_e,gradient_a);
    } // j
  } // i

  for( i=m_grid.xs;i<m_grid.xs+m_grid.xm;i++){
    for( j=m_grid.ys;j<m_grid.ys+m_grid.ym;j++){
      PetscReal dtauc_dzeta;
      m_tauc_param.toTauc(zeta_a[i][j],NULL,&dtauc_dzeta);
      gradient_a[i][j] *= dtauc_dzeta;
    }
  }
  
  if(m_fixed_tauc_locations) {
    DirichletData dirichletBC;
    ierr = dirichletBC.init(m_fixed_tauc_locations); CHKERRQ(ierr);
    dirichletBC.fixResidualHomogeneous(gradient_a);
    ierr = dirichletBC.finish(); CHKERRQ(ierr);
  }

  return 0;
}
PetscErrorCode InvSSATikhonov::domainIP(IceModelVec2S &a, IceModelVec2S &b, PetscReal *OUTPUT) {
  PetscErrorCode ierr;
  ierr = m_domainIP->dot(a,b,OUTPUT); CHKERRQ(ierr);
  return 0;
}

PetscErrorCode InvSSATikhonov::rangeIP(IceModelVec2V &a, IceModelVec2V &b, PetscReal *OUTPUT) {
  PetscErrorCode ierr;
  ierr = dynamic_cast<IPFunctional<IceModelVec2V> *>(m_penaltyFunctional.get())->dot(a,b,OUTPUT); CHKERRQ(ierr);
  return 0;
}

void InvSSATikhonovAddListener(TikhonovProblem<InvSSATikhonov> &problem, 
                 PythonTikhonovSVListener::Ptr listener ) {
  std::tr1::shared_ptr<InvSSAPythonListenerBridge> bridge(new InvSSAPythonListenerBridge(listener) );
  problem.addListener(bridge);
}
