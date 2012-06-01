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
#include "MeanSquareObservationFunctional.hh"

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
  
  PetscReal cL2 = config.get("inv_ssa_domain_l2_coeff");
  PetscReal cH1 = config.get("inv_ssa_domain_h1_coeff");
  
  m_designFunctional.reset(new H1NormFunctional2S(m_grid,cL2,cH1,m_fixed_tauc_locations));    

  m_penaltyFunctional.reset(new MeanSquareObservationFunctional2V(m_grid,m_misfit_weight));    
  (reinterpret_cast<MeanSquareObservationFunctional2V&>(*m_penaltyFunctional)).normalize();  

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
  PetscReal **zeta_a;  
  ierr = m_zeta->get_array(zeta_a);CHKERRQ(ierr);
  ierr = m_tauc->begin_access(); CHKERRQ(ierr);
  for (i=grid.xs; i<m_grid.xs+m_grid.xm; i++) {
    for (j=grid.ys;j<m_grid.ys+m_grid.ym; j++) {
      PetscReal tmp;
      m_tauc_param.toTauc((*m_zeta)(i,j),&tmp,NULL);
      (*tauc)(i,j) = tmp;
    }
  }  
  ierr = m_tauc->end_access(); CHKERRQ(ierr);
  ierr = m_zeta->end_access(); CHKERRQ(ierr);

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


PetscErrorCode InvSSATikhonov::evalObjective(IceModelVec2S &dzeta, PetscReal *OUTPUT) {
  PetscErrorCode ierr;
  ierr = m_designFunctional->valueAt(dzeta,OUTPUT);
  return 0;
}

PetscErrorCode InvSSATikhonov::evalGradObjective(IceModelVec2S &dzeta, IceModelVec2S &gradient) {
  PetscErrorCode ierr;
  ierr = m_designFunctional->gradientAt(dzeta,gradient);
  return 0;
}

PetscErrorCode InvSSATikhonov::evalPenalty(IceModelVec2V &du, PetscReal *OUTPUT) {
  PetscErrorCode ierr;
  ierr = m_penaltyFunctional->valueAt(du,OUTPUT);
  return 0;
}

PetscErrorCode InvSSATikhonov::evalGradPenaltyReduced(IceModelVec2V &du, IceModelVec2S &gradient) {
  PetscErrorCode ierr;

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

  m_penaltyFunctional->gradientAt(du,m_adjointRHS);
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



  PetscInt         i,j;

  PISMVector2 **v_a;
  PISMVector2 v_e[FEQuadrature::Nk];
  PISMVector2 v_q[FEQuadrature::Nq];

  PISMVector2 u_e[FEQuadrature::Nk];
  PISMVector2 u_q[FEQuadrature::Nq];

  PetscReal **gradient_a;
  PetscReal gradient_e[FEQuadrature::Nk];

  ierr = m_v.get_array(v_a); CHKERRQ(ierr);
  ierr = gradient.get_array(gradient_a); CHKERRQ(ierr);

  // An Nq by Nk array of test function values.
  const FEFunctionGerm (*test)[FEQuadrature::Nk] = m_quadrature.testFunctionValues();

  // Jacobian times weights for quadrature.
  PetscScalar JxW[FEQuadrature::Nq];
  m_quadrature.getWeightedJacobian(JxW);

  // Mask (query?) for determining where ice is grounded.
  Mask M;

  // Loop through LOCAL elements.
  PetscInt xs = m_element_index.lxs, xm = m_element_index.lxm,
           ys = m_element_index.lys, ym = m_element_index.lym;
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

        // Determine "beta" at the quadrature point
        PetscReal notquitebeta = 0;
        if( M.grounded_ice(feS->mask) ) {
          notquitebeta = basal.drag(1.,u_qq.u,u_qq.v);
        }

        for (PetscInt k=0; k<FEQuadrature::Nk; k++) {
          gradient_e[k] += -JxW[q]*notquitebeta*(v_qq.u*u_qq.u+v_qq.v*u_qq.v)*test[q][k].val;
        }
      } // q
      m_dofmap.addLocalResidualBlock(gradient_e,gradient_a);
    } // j
  } // i

  PetscReal **zeta_a;
  ierr = m_zeta->get_array(zeta_a); CHKERRQ(ierr);
  for( i=m_grid.xs;i<m_grid.xs+m_grid.xm;i++){
    for( j=m_grid.ys;j<m_grid.ys+m_grid.ym;j++){
      PetscReal dtauc_dzeta;
      m_tauc_param.toTauc(zeta_a[i][j],NULL,&dtauc_dzeta);
      gradient_a[i][j] *= dtauc_dzeta;
    }
  }
  ierr = m_zeta->end_access(); CHKERRQ(ierr);

  if(m_fixed_tauc_locations) {
    DirichletData dirichletBC;
    ierr = dirichletBC.init(m_fixed_tauc_locations); CHKERRQ(ierr);
    dirichletBC.fixResidualHomogeneous(gradient_a);
    ierr = dirichletBC.finish(); CHKERRQ(ierr);
  }

  ierr = m_v.end_access(); CHKERRQ(ierr);
  ierr = m_u.end_access(); CHKERRQ(ierr);
  ierr = gradient.end_access(); CHKERRQ(ierr);
  return 0;
}
