// Copyright (C) 2009--2011 Jed Brown and Ed Bueler and Constantine Khroulev and David Maxwell
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

#include "SSAFEM_Forward.hh"

PetscErrorCode SSAFEM_Forward::allocate_ksp()
{
  PetscErrorCode ierr;

  // Storage for vector unknowns.
  ierr = DAGetMatrix(SSADA, "baij", &m_MatA); CHKERRQ(ierr);
  // ierr = DACreateLocalVector(SSADA, &m_VecU); CHKERRQ(ierr);
  ierr = DACreateGlobalVector(SSADA, &m_VecZ); CHKERRQ(ierr);
  ierr = DACreateGlobalVector(SSADA, &m_VecRHS2); CHKERRQ(ierr);

  // Storage for scalar unknowns.
  ierr = DAGetMatrix(grid.da2, "baij", &m_MatB); CHKERRQ(ierr);
  printf("After making B: %ld\n",(long int)m_MatB);
  ierr = DACreateGlobalVector(grid.da2, &m_VecRHS); CHKERRQ(ierr);
  ierr = DACreateGlobalVector(grid.da2, &m_VecV); CHKERRQ(ierr);
  
  ierr = KSPCreate(grid.com, &m_KSP); CHKERRQ(ierr);
  PC pc;
  ierr = KSPGetPC(m_KSP,&pc); CHKERRQ(ierr);
  ierr = PCSetType(pc,PCBJACOBI); CHKERRQ(ierr);
  ierr = KSPSetFromOptions(m_KSP); CHKERRQ(ierr);

  ierr = KSPCreate(grid.com, &m_KSP_B); CHKERRQ(ierr);
  PC pc_B;
  ierr = KSPGetPC(m_KSP_B,&pc_B); CHKERRQ(ierr);
  ierr = PCSetType(pc_B,PCBJACOBI); CHKERRQ(ierr);
  ierr = KSPSetFromOptions(m_KSP_B); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode SSAFEM_Forward::deallocate_ksp()
{
  PetscErrorCode ierr;

  if (m_KSP != PETSC_NULL) {
    ierr = KSPDestroy(m_KSP); CHKERRQ(ierr);
  }

  if (m_KSP_B != PETSC_NULL) {
    ierr = KSPDestroy(m_KSP_B); CHKERRQ(ierr);
  }

  if (m_MatA != PETSC_NULL) {
    ierr = MatDestroy(m_MatA); CHKERRQ(ierr);
  }

  if (m_MatB != PETSC_NULL) {
    ierr = MatDestroy(m_MatB); CHKERRQ(ierr);
  }

  if (m_VecRHS2 != PETSC_NULL) {
    ierr = VecDestroy(m_VecRHS2); CHKERRQ(ierr);
  }

  if (m_VecRHS != PETSC_NULL) {
    ierr = VecDestroy(m_VecRHS); CHKERRQ(ierr);
  }

  if (m_VecZ != PETSC_NULL) {
    ierr = VecDestroy(m_VecZ); CHKERRQ(ierr);
  }

  if (m_VecV != PETSC_NULL) {
    ierr = VecDestroy(m_VecV); CHKERRQ(ierr);
  }

  return 0;
}

PetscErrorCode SSAFEM_Forward::set_initial_velocity_guess(  IceModelVec2V &v )
{
  PetscErrorCode ierr;
  ierr = v.copy_to(SSAX); CHKERRQ(ierr);
  return 0;
}

PetscErrorCode SSAFEM_Forward::setup_vars()
{
  PetscErrorCode ierr;
  ierr = setup(); CHKERRQ(ierr);
  ierr = assemble_DomainNorm_matrix(); CHKERRQ(ierr);
  
  return 0;
}

PetscErrorCode SSAFEM_Forward::set_tauc(IceModelVec2S &new_tauc )
{
  PetscErrorCode ierr;
  PetscInt i,j,q;

  PetscReal **tauc_array;
  ierr = new_tauc.get_array(tauc_array);CHKERRQ(ierr);
  PetscInt xs = element_index.xs, xm = element_index.xm,
           ys = element_index.ys, ym = element_index.ym;  
  for (i=xs; i<xs+xm; i++) {
    for (j=ys;j<ys+ym; j++) {
      PetscReal taucq[FEQuadrature::Nq];
      quadrature.computeTrialFunctionValues(i,j,dofmap,tauc_array,taucq);
      const PetscInt ij = element_index.flatten(i,j);
      FEStoreNode *feS = &feStore[4*ij];
      for (q=0; q<4; q++) {
        feS[q].tauc = taucq[q]; // Or exp(taucq[q])!
      }
    }
  }
  ierr = new_tauc.end_access();CHKERRQ(ierr);

  m_reassemble_T_matrix_needed = true;
  m_forward_F_needed = true;

  return 0;
}

PetscErrorCode SSAFEM_Forward::solveF_core()
{
  PetscErrorCode ierr;
  if(m_forward_F_needed)
  {
    // Set the SNES callbacks to call into our compute_local_function and compute_local_jacobian
    // methods via SSAFEFunction and SSAFEJ
    ierr = DASetLocalFunction(SSADA,(DALocalFunction1)SSAFEFunction);CHKERRQ(ierr);
    ierr = DASetLocalJacobian(SSADA,(DALocalFunction1)SSAFEJacobian);CHKERRQ(ierr);
    callback_data.da = SSADA;  callback_data.ssa = this;
    ierr = SNESSetFunction(snes, r,    SNESDAFormFunction,   &callback_data);CHKERRQ(ierr);
    ierr = SNESSetJacobian(snes, J, J, SNESDAComputeJacobian,&callback_data);CHKERRQ(ierr);

    // Solve:
    ierr = SNESSolve(snes,NULL,SSAX);CHKERRQ(ierr);

    // See if it worked.
    SNESConvergedReason reason;
    ierr = SNESGetConvergedReason( snes, &reason); CHKERRQ(ierr);
    if(reason < 0) {
      SETERRQ1(1, 
        "SSAFEM solve failed to converge (SNES reason %s)\n\n", SNESConvergedReasons[reason]);
    }
    verbPrintf(2,grid.com,"SSAFEM solve converged (SNES reason %s)\n\n", SNESConvergedReasons[reason]);
    m_forward_F_needed = false;
  }  
  return 0;
}

PetscErrorCode SSAFEM_Forward::solveF(IceModelVec2V &result)
{
  PetscErrorCode ierr;
  ierr = solveF_core(); CHKERRQ(ierr);

  ierr = result.copy_from(SSAX); CHKERRQ(ierr);
  ierr = result.beginGhostComm(); CHKERRQ(ierr);
  ierr = result.endGhostComm(); CHKERRQ(ierr);
  return 0;
}


PetscErrorCode SSAFEM_Forward::solveT( IceModelVec2S &dtauc, IceModelVec2V &result)
{
  KSPConvergedReason  reason;
  PetscErrorCode ierr;
  PetscScalar **dtauc_a;
  PISMVector2 **vel;
  PISMVector2 **rhs;

  ierr = solveF_core(); CHKERRQ(ierr);

  if(m_reassemble_T_matrix_needed)
  {
    assemble_T_matrix();
  }

  dtauc.get_array(dtauc_a);  
  ierr = DAVecGetArray(SSADA,SSAX,&vel); CHKERRQ(ierr);
  ierr = DAVecGetArray(SSADA,m_VecRHS2,&rhs); CHKERRQ(ierr);

  ierr = assemble_T_rhs(vel,dtauc_a,rhs); CHKERRQ(ierr);

  dtauc.end_access();
  ierr = DAVecRestoreArray(SSADA, SSAX, &vel); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(SSADA, m_VecRHS2, &rhs); CHKERRQ(ierr);

  PetscTruth monitorFunction;
  ierr = PetscOptionsHasName(NULL,"-ssaf_monitor_Trhs",&monitorFunction);CHKERRQ(ierr);
  if (monitorFunction) {
    PetscViewer    viewer;
    ierr = PetscViewerCreate(grid.com, &viewer);CHKERRQ(ierr);
    ierr = PetscViewerSetType(viewer, PETSC_VIEWER_ASCII);CHKERRQ(ierr);
    dtauc.copy_to(m_VecV);
    ierr = VecView(m_VecV,viewer);
    
    // ierr = MatView(m_MatA,viewer); CHKERRQ(ierr);
    ierr = PetscViewerDestroy(viewer);CHKERRQ(ierr);
    
  }


  // call PETSc to solve linear system by iterative method.
  ierr = KSPSetOperators(m_KSP, m_MatA, m_MatA, SAME_NONZERO_PATTERN); CHKERRQ(ierr);
  ierr = KSPSolve(m_KSP, m_VecRHS2, m_VecZ); CHKERRQ(ierr); // SOLVE

  ierr = KSPGetConvergedReason(m_KSP, &reason); CHKERRQ(ierr);
  if (reason < 0) {
    SETERRQ1(1, 
      "SSAFEM_Forward::solveT solve failed to converge (KSP reason %s)\n\n", KSPConvergedReasons[reason]);
  }
  else  
  {
    verbPrintf(2,grid.com,"SSAFEM_Forward::solveT converged (KSP reason %s)\n", KSPConvergedReasons[reason] );
  }

  // Extract the solution and communicate.
  ierr = result.copy_from(m_VecZ); CHKERRQ(ierr);
  ierr = result.beginGhostComm(); CHKERRQ(ierr);
  ierr = result.endGhostComm(); CHKERRQ(ierr);  
  
  return 0;
}

PetscErrorCode SSAFEM_Forward::solveTStar( IceModelVec2V &residual, IceModelVec2S &result)
{
  KSPConvergedReason  reason;
  PetscErrorCode ierr;
  PISMVector2 **R;
  PISMVector2 **U;
  PISMVector2 **RHS2;
  PISMVector2 **Z;
  PetscScalar **RHS;

  // Solve the nonlinear forward problem if this has not yet been done.
  ierr = solveF_core(); CHKERRQ(ierr);

  // If tauc has been updated, we'll need to update the linearized forward map.
  if(m_reassemble_T_matrix_needed)
  {
    assemble_T_matrix();
  }

  // Assemble the right-hand side for the first step.
  residual.get_array(R);  
  ierr = DAVecGetArray(SSADA,m_VecRHS2,&RHS2); CHKERRQ(ierr);  
  ierr = assemble_TStarA_rhs(R,RHS2); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(SSADA, m_VecRHS2, &RHS2); CHKERRQ(ierr);
  residual.end_access();
  

  // call PETSc to solve linear system by iterative method.
  ierr = KSPSetOperators(m_KSP, m_MatA, m_MatA, SAME_NONZERO_PATTERN); CHKERRQ(ierr);
  ierr = KSPSolve(m_KSP, m_VecRHS2, m_VecZ); CHKERRQ(ierr); // SOLVE
  
  ierr = KSPGetConvergedReason(m_KSP, &reason); CHKERRQ(ierr);
  if (reason < 0) {
    SETERRQ1(1, 
      "SSAFEM_Forward::solveTStarA solve failed to converge (KSP reason %s)\n\n", KSPConvergedReasons[reason]);
  }
  else  
  {
    verbPrintf(2,grid.com,"SSAFEM_Forward::solveTStarA converged (KSP reason %s)\n", KSPConvergedReasons[reason] );
  }
  
  // Assemble the right-hand side for the second step.
  ierr = DAVecGetArray(SSADA,m_VecZ,&Z); CHKERRQ(ierr);  
  ierr = DAVecGetArray(SSADA,SSAX,&U); CHKERRQ(ierr);  
  ierr = DAVecGetArray(grid.da2,m_VecRHS,&RHS); CHKERRQ(ierr);  
  ierr = assemble_TStarB_rhs(Z,U,RHS); CHKERRQ(ierr);
  
  ierr = DAVecRestoreArray(grid.da2,m_VecRHS,&RHS); CHKERRQ(ierr);  
  ierr = DAVecRestoreArray(SSADA,SSAX,&U); CHKERRQ(ierr);  
  ierr = DAVecRestoreArray(SSADA,m_VecZ,&Z); CHKERRQ(ierr);  
  
  // call PETSc to solve linear system by iterative method.
  ierr = KSPSetOperators(m_KSP_B, m_MatB, m_MatB, SAME_NONZERO_PATTERN); CHKERRQ(ierr);
  ierr = KSPSolve(m_KSP_B, m_VecRHS, m_VecV); CHKERRQ(ierr); // SOLVE
  ierr = KSPGetConvergedReason(m_KSP, &reason); CHKERRQ(ierr);
  if (reason < 0) {
    SETERRQ1(1, 
      "SSAFEM_Forward::solveTStarB solve failed to converge (KSP reason %s)\n\n", KSPConvergedReasons[reason]);
  }
  else  {
    verbPrintf(2,grid.com,"SSAFEM_Forward::solveTStarB converged (KSP reason %s)\n", KSPConvergedReasons[reason] );
  }
  
  
  // Extract the solution and communicate.
  ierr = result.copy_from(m_VecV); CHKERRQ(ierr);
  ierr = result.beginGhostComm(); CHKERRQ(ierr);
  ierr = result.endGhostComm(); CHKERRQ(ierr);  
  
  return 0;
}

PetscErrorCode SSAFEM_Forward::assemble_T_matrix()
{
  PISMVector2 **vel;
  PetscErrorCode ierr;

  ierr = DAVecGetArray(SSADA,SSAX,&vel); CHKERRQ(ierr);
  
  DALocalInfo *info = NULL;
  ierr = compute_local_jacobian(info,const_cast<const PISMVector2**>(vel),m_MatA);
  
  ierr = DAVecRestoreArray(SSADA, SSAX, &vel); CHKERRQ(ierr);
  return 0;
}

PetscErrorCode SSAFEM_Forward::assemble_DomainNorm_matrix()
{
  // Just L2 matrix for now.
  PetscInt         i,j;
  PetscErrorCode   ierr;

  // Zero out the Jacobian in preparation for updating it.
  ierr = MatZeroEntries(m_MatB);CHKERRQ(ierr);

  // Jacobian times weights for quadrature.
  PetscScalar JxW[FEQuadrature::Nq];
  quadrature.getWeightedJacobian(JxW);

  // Values of the finite element test functions at the quadrature points.
  // This is an Nq by Nk array of function germs (Nq=#of quad pts, Nk=#of test functions).
  const FEFunctionGerm (*test)[FEQuadrature::Nk] = quadrature.testFunctionValues();

  // Loop through all the elements.
  PetscInt xs = element_index.xs, xm = element_index.xm,
           ys = element_index.ys, ym = element_index.ym;
  for (i=xs; i<xs+xm; i++) {
    for (j=ys; j<ys+ym; j++) {
      // Element-local Jacobian matrix (there are FEQuadrature::Nk vector valued degrees 
      // of freedom per elment, for a total of (2*FEQuadrature::Nk)*(2*FEQuadrature::Nk) = 16
      // entries in the local Jacobian.
      PetscReal      K[FEQuadrature::Nk][FEQuadrature::Nk];

      // Initialize the map from global to local degrees of freedom for this element.
      dofmap.reset(i,j,grid);

      // Build the element-local Jacobian.
      ierr = PetscMemzero(K,sizeof(K));CHKERRQ(ierr);
      for (PetscInt q=0; q<FEQuadrature::Nq; q++) {
        for (PetscInt k=0; k<4; k++) {   // Test functions
          for (PetscInt l=0; l<4; l++) { // Trial functions
            const FEFunctionGerm &test_qk=test[q][k];
            const FEFunctionGerm &test_ql=test[q][l];
            K[k][l]     += JxW[q]*test_qk.val*test_ql.val;
          } // l
        } // k
      } // q
      ierr = dofmap.addLocalJacobianBlock(&K[0][0],m_MatB);
    } // j
  } // i
  
  ierr = MatAssemblyBegin(m_MatB,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(m_MatB,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    
  return 0;
}

PetscErrorCode SSAFEM_Forward::domainIP(IceModelVec2S &avec, IceModelVec2S &bvec, PetscScalar *OUTPUT)
{
  // Just L2 matrix for now.
  PetscInt         i,j;
  PetscErrorCode   ierr;

  // The value of the inner product.
  PetscReal IP = 0;

  PetscReal **A, **B;
  PetscReal a[FEQuadrature::Nq], b[FEQuadrature::Nq];

  avec.get_array(A);
  bvec.get_array(B);

  // Jacobian times weights for quadrature.
  PetscScalar JxW[FEQuadrature::Nq];
  quadrature.getWeightedJacobian(JxW);

  // Loop through all LOCAL elements.
  PetscInt xs = element_index.lxs, xm = element_index.lxm,
           ys = element_index.lys, ym = element_index.lym;
  for (i=xs; i<xs+xm; i++) {
    for (j=ys; j<ys+ym; j++) {
      quadrature.computeTrialFunctionValues(i,j,dofmap,A,a);
      quadrature.computeTrialFunctionValues(i,j,dofmap,B,b);
      for (PetscInt q=0; q<FEQuadrature::Nq; q++) {
        IP += JxW[q]*a[q]*b[q];
      } // q
    } // j
  } // i
  
  ierr = PetscGlobalSum(&IP, OUTPUT, grid.com); CHKERRQ(ierr);

  avec.end_access();
  bvec.end_access();
  return 0;
}

PetscErrorCode SSAFEM_Forward::rangeIP(IceModelVec2V &avec, IceModelVec2V &bvec, PetscScalar *OUTPUT)
{
  PetscInt         i,j;
  PetscErrorCode   ierr;

  // The value of the inner product.
  PetscReal IP = 0;

  PISMVector2 **A, **B;
  PISMVector2 a[FEQuadrature::Nq], b[FEQuadrature::Nq];

  avec.get_array(A);
  bvec.get_array(B);

  // Jacobian times weights for quadrature.
  PetscScalar JxW[FEQuadrature::Nq];
  quadrature.getWeightedJacobian(JxW);

  // Loop through all LOCAL elements.
  PetscInt xs = element_index.lxs, xm = element_index.lxm,
           ys = element_index.lys, ym = element_index.lym;
  for (i=xs; i<xs+xm; i++) {
    for (j=ys; j<ys+ym; j++) {
      quadrature.computeTrialFunctionValues(i,j,dofmap,A,a);
      quadrature.computeTrialFunctionValues(i,j,dofmap,B,b);
      for (PetscInt q=0; q<FEQuadrature::Nq; q++) {
        IP += JxW[q]*(a[q].u*b[q].u+a[q].v+b[q].v);
      } // q
    } // j
  } // i
  
  ierr = PetscGlobalSum(&IP, OUTPUT, grid.com); CHKERRQ(ierr);

  avec.end_access();
  bvec.end_access();
  return 0;
}

PetscErrorCode SSAFEM_Forward::assemble_T_rhs( PISMVector2 **gvel, PetscReal **gdtauc, PISMVector2 **grhs)
{  
  PetscInt         i,j,k,q;
  PetscErrorCode   ierr;
  PetscReal        **bc_mask;

  // Zero out the portion of the function we are responsible for computing.
  for (i=grid.xs; i<grid.xs+grid.xm; i++) {
    for (j=grid.ys; j<grid.ys+grid.ym; j++) {
      grhs[i][j].u = 0;
      grhs[i][j].v = 0;
    }
  }

  // Start access of Dirichlet data, if present.
  if (bc_locations && vel_bc) {
    ierr = bc_locations->get_array(bc_mask);CHKERRQ(ierr);
  }

  // Jacobian times weights for quadrature.
  PetscScalar JxW[FEQuadrature::Nq];
  quadrature.getWeightedJacobian(JxW);

  // Storage for velocity and dtauc at quadrature points.
  PISMVector2 u[FEQuadrature::Nq];
  PetscReal dtauc[FEQuadrature::Nq];

  // An Nq by Nk array of test function values.
  const FEFunctionGerm (*test)[FEQuadrature::Nk] = quadrature.testFunctionValues();

  // Flags for each vertex in an element that determine if explicit Dirichlet data has
  // been set.
  PetscReal local_bc_mask[FEQuadrature::Nk];

  // Iterate over the elements.
  PetscInt xs = element_index.xs, xm = element_index.xm,
           ys = element_index.ys, ym = element_index.ym;
  for (i=xs; i<xs+xm; i++) {
    for (j=ys; j<ys+ym; j++) {
      // Storage for element-local data
      PISMVector2     x[4],y[4]; 
      PetscReal       w[4];

      // Index into coefficient storage in feStore
      const PetscInt ij = element_index.flatten(i,j);

      // Initialize the map from global to local degrees of freedom for this element.
      dofmap.reset(i,j,grid);

      // Obtain the value of the solution at the adjacent nodes to the element.
      dofmap.extractLocalDOFs(i,j,gvel,x);

      // These values now need to be adjusted if some nodes in the element have
      // Dirichlet data.
      if(bc_locations && vel_bc) {
        dofmap.extractLocalDOFs(i,j,bc_mask,local_bc_mask);
        for (k=0; k<4; k++) {
          if (PismIntMask(local_bc_mask[k]) == 1) { // Dirichlet node
            PetscInt ii, jj;
            dofmap.localToGlobal(k,&ii,&jj);
            x[k].u = 0;
            x[k].v = 0;
            // Mark any kind of Dirichlet node as not to be touched
            dofmap.markRowInvalid(k);
            dofmap.markColInvalid(k);
          }
        }
      }

      quadrature.computeTrialFunctionValues(x,u);

      dofmap.extractLocalDOFs(i,j,gdtauc,w);
      quadrature.computeTrialFunctionValues(w,dtauc);

      // Zero out the element-local residual in prep for updating it.
      for(k=0;k<FEQuadrature::Nk;k++){ 
        y[k].u = 0; y[k].v = 0;
      }

      for (q=0; q<FEQuadrature::Nq; q++) {     // loop over quadrature points on this element.

        // Coefficients and weights for this quadrature point.
        const FEStoreNode *feS = &feStore[ij*FEQuadrature::Nq+q];
        const PetscReal    jw  = JxW[q];

        // Determine dbeta/dtauc at the quadrature point
        PetscReal dbeta_dtauc = 0;
        if(feS->mask == MASK_GROUNDED ) {
          dbeta_dtauc = basal.drag(dtauc[q],u[q].u,u[q].v);
        }
        // dbeta_dtauc *= exp(dtauc[q]) 

        for(k=0; k<FEQuadrature::Nk;k++) {  // loop over the test functions.
          const FEFunctionGerm &testqk = test[q][k];
          //Note the -= (not +=) in the following lines.
          y[k].u -= jw*dbeta_dtauc*testqk.val*u[q].u;
          y[k].v -= jw*dbeta_dtauc*testqk.val*u[q].v;
        } // k
      } // q

      dofmap.addLocalResidualBlock(y,grhs);
    } // j-loop
  } // i-loop

  // Until now we have not touched rows in the residual corresponding to Dirichlet data.
  // We fix this now.
  if (bc_locations && vel_bc) {
    // Enforce Dirichlet conditions strongly
    for (i=grid.xs; i<grid.xs+grid.xm; i++) {
      for (j=grid.ys; j<grid.ys+grid.ym; j++) {
        if (bc_locations->value(i,j) == 1) {
          // Enforce explicit homogeneous dirichlet data.
          grhs[i][j].u = 0;
          grhs[i][j].v = 0;
        }
      }
    }
    ierr = bc_locations->end_access();CHKERRQ(ierr);
  }

  return 0;
}

PetscErrorCode SSAFEM_Forward::assemble_TStarA_rhs( PISMVector2 **R, PISMVector2 **RHS)
{  
  PetscInt         i,j,k,q;
  PetscErrorCode   ierr;
  PetscReal        **bc_mask;

  // Zero out the portion of the function we are responsible for computing.
  for (i=grid.xs; i<grid.xs+grid.xm; i++) {
    for (j=grid.ys; j<grid.ys+grid.ym; j++) {
      RHS[i][j].u = 0;
      RHS[i][j].v = 0;
    }
  }

  // Start access of Dirichlet data, if present.
  if (bc_locations && vel_bc) {
    ierr = bc_locations->get_array(bc_mask);CHKERRQ(ierr);
  }

  // Jacobian times weights for quadrature.
  PetscScalar JxW[FEQuadrature::Nq];
  quadrature.getWeightedJacobian(JxW);

  // Storage for R at quadrature points.
  PISMVector2 res[FEQuadrature::Nq];

  // An Nq by Nk array of test function values.
  const FEFunctionGerm (*test)[FEQuadrature::Nk] = quadrature.testFunctionValues();

  // Iterate over the elements.
  PetscInt xs = element_index.xs, xm = element_index.xm,
           ys = element_index.ys, ym = element_index.ym;
  for (i=xs; i<xs+xm; i++) {
    for (j=ys; j<ys+ym; j++) {
      // Storage for element-local data
      PISMVector2 x[FEQuadrature::Nk], y[FEQuadrature::Nk];

      // Initialize the map from global to local degrees of freedom for this element.
      dofmap.reset(i,j,grid);

      // Obtain the value of the solution at the adjacent nodes to the element.
      dofmap.extractLocalDOFs(i,j,R,x);
      quadrature.computeTrialFunctionValues(x,res);

      // Zero out the element-local residual in prep for updating it.
      for(k=0;k<FEQuadrature::Nk;k++){ 
        y[k].u = 0; y[k].v = 0;
      }

      for (q=0; q<FEQuadrature::Nq; q++) {     // loop over quadrature points on this element.
        const FEFunctionGerm &testqk = test[q][k];
        // Coefficients and weights for this quadrature point.
        const PetscReal    jw  = JxW[q];
        for(k=0; k<FEQuadrature::Nk;k++) {  // loop over the test functions.
          y[k].u += jw*testqk.val*res[q].u;
          y[k].v += jw*testqk.val*res[q].v;
        }
      } // q

      dofmap.addLocalResidualBlock(y,RHS);
    } // j-loop
  } // i-loop

  // Until now we have not touched rows in the residual corresponding to Dirichlet data.
  // We fix this now.
  if (bc_locations && vel_bc) {
    // Enforce Dirichlet conditions strongly
    for (i=grid.xs; i<grid.xs+grid.xm; i++) {
      for (j=grid.ys; j<grid.ys+grid.ym; j++) {
        if (bc_locations->value(i,j) == 1) {
          // Enforce explicit homogeneous dirichlet data.
          RHS[i][j].u = 0;
          RHS[i][j].v = 0;
        }
      }
    }
    ierr = bc_locations->end_access();CHKERRQ(ierr);
  }

  return 0;
}


PetscErrorCode SSAFEM_Forward::assemble_TStarB_rhs( PISMVector2 **Z, 
                                                    PISMVector2 **U,
                                                    PetscScalar **RHS)
{
  PetscInt         i,j,k,q;
  PetscErrorCode   ierr;

  // Zero out the portion of the function we are responsible for computing.
  for (i=grid.xs; i<grid.xs+grid.xm; i++) {
    for (j=grid.ys; j<grid.ys+grid.ym; j++) {
      RHS[i][j] = 0;
    }
  }

  // Jacobian times weights for quadrature.
  PetscScalar JxW[FEQuadrature::Nq];
  quadrature.getWeightedJacobian(JxW);

  // Storage for z, u at quadrature points.
  PISMVector2 z[FEQuadrature::Nq];
  PISMVector2 u[FEQuadrature::Nq];

  // An Nq by Nk array of test function values.
  const FEFunctionGerm (*test)[FEQuadrature::Nk] = quadrature.testFunctionValues();

  // Iterate over the elements.
  PetscInt xs = element_index.xs, xm = element_index.xm,
           ys = element_index.ys, ym = element_index.ym;
  for (i=xs; i<xs+xm; i++) {
    for (j=ys; j<ys+ym; j++) {
      // Storage for element-local data
      PISMVector2 x[FEQuadrature::Nk];
      PetscScalar y[FEQuadrature::Nk];

      // Index into coefficient storage in feStore
      const PetscInt ij = element_index.flatten(i,j);

      // Initialize the map from global to local degrees of freedom for this element.
      dofmap.reset(i,j,grid);

      // Obtain the values of Z and VEL at the element quad points
      dofmap.extractLocalDOFs(i,j,Z,x);
      quadrature.computeTrialFunctionValues(x,z);
      dofmap.extractLocalDOFs(i,j,U,x);
      quadrature.computeTrialFunctionValues(x,u);

      // Zero out the element-local residual in prep for updating it.
      for(k=0;k<FEQuadrature::Nk;k++){ 
        y[k] = 0;
      }

      for (q=0; q<FEQuadrature::Nq; q++) {     // loop over quadrature points on this element.
        const FEFunctionGerm &testqk = test[q][k];
        // Coefficients and weights for this quadrature point.
        const FEStoreNode *feS = &feStore[ij*FEQuadrature::Nq+q];
        const PetscReal    jw  = JxW[q];

        // Determine "beta" at the quadrature point
        PetscReal notquitebeta = 0;
        if(feS->mask == MASK_GROUNDED ) {
          notquitebeta = basal.drag(1.,u[q].u,u[q].v);
        }
        // dbeta_dtauc *= exp(dtauc[q]) 

        for(k=0; k<FEQuadrature::Nk;k++) {  // loop over the test functions.
          y[k] -= jw*notquitebeta*testqk.val*(u[q].u*z[q].u+u[q].v*z[q].v);
        }
      } // q

      dofmap.addLocalResidualBlock(y,RHS);
    } // j-loop
  } // i-loop
  return 0;  

  ierr = MatAssemblyBegin(m_MatB,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(m_MatB,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

  return 0;
}

