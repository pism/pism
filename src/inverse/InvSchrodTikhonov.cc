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

#include "InvSchrodTikhonov.hh"
#include "H1NormFunctional.hh"

InvSchrodTikhonov::InvSchrodTikhonov( IceGrid  &grid, IceModelVec2V &f) :
m_grid(grid), m_c(NULL), m_f(&f), 
m_dirichletLocations(NULL), m_dirichletValues(NULL), m_dirichletWeight(1.), 
m_fixedDesignLocations(NULL), m_element_index(grid) {
  PetscErrorCode ierr = this->construct();
  CHKERRCONTINUE(ierr);
  if(ierr) {
    PetscPrintf(m_grid.com, "FATAL ERROR: InvSchrodTikhonov allocation failed.\n");
    PISMEnd();    
  }
}

InvSchrodTikhonov::~InvSchrodTikhonov( ) {
  PetscErrorCode ierr = this->destruct();
  CHKERRCONTINUE(ierr);
  if(ierr) {
    PetscPrintf(m_grid.com, "FATAL ERROR: InvSchrodTikhonov deallocation failed.\n");
    PISMEnd();    
  }
}

PetscErrorCode InvSchrodTikhonov::construct() {
  PetscErrorCode ierr;
  PetscInt stencilWidth = 1;

  // m_designFunctional.reset(new L2NormFunctional2S(m_grid));

  m_uGlobal.create(m_grid,"Schrodinger solution (sans ghosts)",kNoGhosts,stencilWidth);
  m_u.create(m_grid,"Schrodinger solution",kHasGhosts,stencilWidth);
  m_r.create(m_grid,"Schrodinger residual",kNoGhosts,stencilWidth);

  m_vGlobal.create(m_grid,"adjoint work vector (sans ghosts)",kNoGhosts,stencilWidth);
  m_v.create(m_grid,"adjoint work vector",kHasGhosts,stencilWidth);

  m_adjointRHS.create(m_grid,"adjoint RHS",kNoGhosts,stencilWidth);
  

  // Create a DA with 2 degres of freedom.
  PetscInt dof=2, stencil_width=1;
  ierr = DMDACreate2d( m_grid.com,
                      DMDA_BOUNDARY_PERIODIC, DMDA_BOUNDARY_PERIODIC,
                      DMDA_STENCIL_BOX,
                      m_grid.My, m_grid.Mx,
                      m_grid.Ny, m_grid.Nx,
                      dof, stencil_width,
                      &m_grid.procs_y[0], &m_grid.procs_x[0],
                      &m_da); CHKERRQ(ierr);

  ierr = DMGetMatrix(m_da, "baij", &m_J); CHKERRQ(ierr);

  ierr = KSPCreate(m_grid.com, &m_ksp); CHKERRQ(ierr);
  PetscReal ksp_rtol = 1e-12;
  ierr = KSPSetTolerances(m_ksp,ksp_rtol,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);
  PC pc;
  ierr = KSPGetPC(m_ksp,&pc); CHKERRQ(ierr);
  ierr = PCSetType(pc,PCBJACOBI); CHKERRQ(ierr);
  ierr = KSPSetFromOptions(m_ksp); CHKERRQ(ierr);  

  ierr = SNESCreate(m_grid.com, &m_snes);CHKERRQ(ierr);

  // Default of maximum 200 iterations; possibly overridded by commandline
  PetscInt snes_max_it = 200;
  ierr = SNESSetTolerances(m_snes,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT,
    snes_max_it,PETSC_DEFAULT); CHKERRQ(ierr);
  ierr = SNESSetFromOptions(m_snes);CHKERRQ(ierr);

  ierr = m_callbacks.connect(m_snes,*this,m_da,m_r.get_vec(),m_J); CHKERRQ(ierr);

  m_quadrature.init(m_grid);
  return 0;
}

PetscErrorCode InvSchrodTikhonov::solve(bool &success) {
  PetscErrorCode ierr;

  PetscReal cL2 = 1;
  ierr = PetscOptionsReal("-inv_schrod_cL2",
                          "design variable L2 weight",
                          "",
                          cL2,
                          &cL2,NULL);CHKERRQ(ierr);
  PetscReal cH1 = 1;
  ierr = PetscOptionsReal("-inv_schrod_cH1",
                          "design variable H1 weight",
                          "",
                          cH1,
                          &cH1,NULL);CHKERRQ(ierr);

  m_designFunctional.reset(new H1NormFunctional2S(m_grid,cL2,cH1,m_fixedDesignLocations));    

  m_penaltyFunctional.reset(new L2NormFunctional2V(m_grid));

  success=false;
  ierr = SNESSolve(m_snes,NULL,m_uGlobal.get_vec()); CHKERRQ(ierr);

  ierr = SNESGetConvergedReason( m_snes, &m_reason); CHKERRQ(ierr);
  if(m_reason>0){
    success = true;
  } else {

    KSP ksp;
    ierr = SNESGetKSP(m_snes,&ksp); CHKERRQ(ierr);

    PetscInt n = 40, neig;
    
    PetscReal *r, *c;
    ierr = PetscMalloc(2*n*sizeof(PetscReal),&r);CHKERRQ(ierr);
    c = r + n;
    ierr = KSPComputeEigenvalues(ksp,n,r,c,&neig);CHKERRQ(ierr);
    if (!m_grid.rank) {
      PetscDraw   draw;
      PetscDrawSP drawsp;
      PetscViewer viewer;

      ierr = PetscViewerDrawOpen(PETSC_COMM_SELF,0,"Iteratively Computed Eigenvalues",PETSC_DECIDE,PETSC_DECIDE,300,300,&viewer);CHKERRQ(ierr);
      ierr = PetscViewerDrawGetDraw(viewer,0,&draw);CHKERRQ(ierr);
      ierr = PetscDrawSPCreate(draw,1,&drawsp);CHKERRQ(ierr);
      for (PetscInt i=0; i<neig; i++) {
        ierr = PetscDrawSPAddPoint(drawsp,r+i,c+i);CHKERRQ(ierr);
      }
      ierr = PetscDrawSPDraw(drawsp);CHKERRQ(ierr);
      ierr = PetscDrawSPDestroy(&drawsp);CHKERRQ(ierr);
      ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);
    }
    ierr = PetscFree(r);CHKERRQ(ierr);

    PetscViewer viewer;
    ierr = PetscViewerASCIIOpen(m_grid.com,"InvSchrodTikhonovFailed.m",&viewer); CHKERRQ(ierr);
    ierr = PetscViewerSetFormat(viewer, PETSC_VIEWER_ASCII_MATLAB); CHKERRQ(ierr);
    ierr = MatView(m_J,viewer); CHKERRQ(ierr);
    ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);    
  }
  m_u.copy_from(m_uGlobal);
  return 0;
}

PetscErrorCode InvSchrodTikhonov::destruct() {
  PetscErrorCode ierr;
  ierr = DMDestroy(&m_da); CHKERRQ(ierr);
  ierr = MatDestroy(&m_J); CHKERRQ(ierr);
  ierr = SNESDestroy(&m_snes); CHKERRQ(ierr);
  ierr = KSPDestroy(&m_ksp); CHKERRQ(ierr);
  return 0;
}

PetscErrorCode InvSchrodTikhonov::linearizeAt( IceModelVec2S &c, bool &success) {
  this->set_c(c);
  PetscErrorCode ierr;
  ierr = this->solve(success); CHKERRQ(ierr);
  return 0;
}

PetscErrorCode InvSchrodTikhonov::assembleFunction( DMDALocalInfo *info, PISMVector2 **xg, PISMVector2 **yg) {
  PetscInt         i,j,k,q;
  PetscErrorCode   ierr;

  (void) info; // Avoid compiler warning.

  // Zero out the portion of the function we are responsible for computing.
  for (i=m_grid.xs; i<m_grid.xs+m_grid.xm; i++) {
    for (j=m_grid.ys; j<m_grid.ys+m_grid.ym; j++) {
      yg[i][j].u = yg[i][j].v = 0;
    }
  }

  // Start access of Dirichlet data, if present.
  DirichletData dirichletBC;
  ierr = dirichletBC.init(m_dirichletLocations,m_dirichletValues,m_dirichletWeight); CHKERRQ(ierr);

  // Jacobian times weights for quadrature.
  PetscScalar JxW[FEQuadrature::Nq];
  m_quadrature.getWeightedJacobian(JxW);

  // Storage for the current solution at quadrature points.
  PISMVector2 u_q[FEQuadrature::Nq];
  PISMVector2 Dux_q[FEQuadrature::Nq];
  PISMVector2 Duy_q[FEQuadrature::Nq];

  PISMVector2 **f;
  ierr = m_f->get_array(f); CHKERRQ(ierr);
  PISMVector2 f_q[FEQuadrature::Nq];

  PetscReal **c;
  ierr = m_c->get_array(c); CHKERRQ(ierr);
  PetscReal   c_q[FEQuadrature::Nq];

  // An Nq by Nk array of test function values.
  const FEFunctionGerm (*test)[FEQuadrature::Nk] = m_quadrature.testFunctionValues();

  // Iterate over the elements.
  PetscInt xs = m_element_index.xs, xm = m_element_index.xm,
           ys = m_element_index.ys, ym = m_element_index.ym;
  for (i=xs; i<xs+xm; i++) {
    for (j=ys; j<ys+ym; j++) {
      // Storage for element-local solution and residuals.
      PISMVector2     x_e[FEQuadrature::Nk], y_e[FEQuadrature::Nk],
                      f_e[FEQuadrature::Nk];
      PetscReal       c_e[FEQuadrature::Nk];

      // Initialize the map from global to local degrees of freedom for this element.
      m_dofmap.reset(i,j,m_grid);

      // Obtain the value of the solution at the nodes adjacent to the element.
      m_dofmap.extractLocalDOFs(i,j,xg,x_e);
      m_dofmap.extractLocalDOFs(i,j,f,f_e);
      m_dofmap.extractLocalDOFs(i,j,c,c_e);

      if(dirichletBC) dirichletBC.update(m_dofmap,x_e);

      // Zero out the element-local residual in prep for updating it.
      for(k=0;k<FEQuadrature::Nk;k++){
        y_e[k].u = 0; y_e[k].v = 0;
      }

      // Compute the solution values and symmetric gradient at the quadrature points.
      m_quadrature.computeTrialFunctionValues(x_e,u_q,Dux_q,Duy_q);
      m_quadrature.computeTrialFunctionValues(f_e,f_q);
      m_quadrature.computeTrialFunctionValues(c_e,c_q);

      for (q=0; q<FEQuadrature::Nq; q++) {     // loop over quadrature points on this element.

        // Symmetric gradient at the quadrature point.
        PISMVector2 Dux_qq = Dux_q[q];
        PISMVector2 Duy_qq = Duy_q[q];

        const PetscReal    jw  = JxW[q];

        PetscReal Fu_q =c_q[q]*u_q[q].u-f_q[q].u;
        PetscReal Fv_q =c_q[q]*u_q[q].v-f_q[q].v;
        for(k=0; k<4;k++) {  // loop over the test functions.
          const FEFunctionGerm &testqk = test[q][k];
          y_e[k].u += jw*(testqk.dx*Dux_qq.u+testqk.dy*Duy_qq.u + testqk.val*Fu_q);
          y_e[k].v += jw*(testqk.dx*Dux_qq.v+testqk.dy*Duy_qq.v + testqk.val*Fv_q);
        }
      } // q
      m_dofmap.addLocalResidualBlock(y_e,yg);
    } // j-loop
  } // i-loop

  // Until now we have not touched rows in the residual corresponding to Dirichlet data.
  // We fix this now.
  if(dirichletBC) dirichletBC.fixResidual(xg,yg);
  ierr = dirichletBC.finish(); CHKERRQ(ierr);

  ierr = m_f->end_access(); CHKERRQ(ierr);
  ierr = m_c->end_access(); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode InvSchrodTikhonov::assembleJacobian( DMDALocalInfo *info, PISMVector2 **x, Mat J) {
  PetscInt         i,j;
  PetscErrorCode   ierr;

  // Avoid compiler warning.
  (void) info;
  (void) x;

  PetscReal **c;
  ierr = m_c->get_array(c); CHKERRQ(ierr);
  PetscReal   c_q[FEQuadrature::Nq];
  
  // Zero out the Jacobian in preparation for updating it.
  ierr = MatZeroEntries(J);CHKERRQ(ierr);

  DirichletData dirichletBC;
  dirichletBC.init(m_dirichletLocations, m_dirichletValues, m_dirichletWeight);

  // Jacobian times weights for quadrature.
  PetscScalar JxW[FEQuadrature::Nq];
  m_quadrature.getWeightedJacobian(JxW);

  // Values of the finite element test functions at the quadrature points.
  // This is an Nq by Nk array of function germs (Nq=#of quad pts, Nk=#of test functions).
  const FEFunctionGerm (*test)[FEQuadrature::Nk] = m_quadrature.testFunctionValues();

  // Loop through all the elements.
  PetscInt xs = m_element_index.xs, xm = m_element_index.xm,
           ys = m_element_index.ys, ym = m_element_index.ym;
  for (i=xs; i<xs+xm; i++) {
    for (j=ys; j<ys+ym; j++) {
      // Values of the solution at the nodes of the current element.
      PISMVector2    x_e[FEQuadrature::Nk];
      PetscReal      c_e[FEQuadrature::Nk];

      // Element-local Jacobian matrix (there are FEQuadrature::Nk vector valued degrees
      // of freedom per elment, for a total of (2*FEQuadrature::Nk)*(2*FEQuadrature::Nk) = 16
      // entries in the local Jacobian.
      PetscReal      K[(2*FEQuadrature::Nk)*(2*FEQuadrature::Nk)];

      // Initialize the map from global to local degrees of freedom for this element.
      m_dofmap.reset(i,j,m_grid);

      // Obtain the value of the solution at the adjacent nodes to the element.
      m_dofmap.extractLocalDOFs(i,j, x,x_e);
      m_dofmap.extractLocalDOFs(i,j, c,c_e);
      
      if(dirichletBC) dirichletBC.update(m_dofmap,x_e);

      // Compute the values of the solution at the quadrature points.
      m_quadrature.computeTrialFunctionValues(c_e,c_q);

      // Build the element-local Jacobian.
      ierr = PetscMemzero(K,sizeof(K));CHKERRQ(ierr);
      for (PetscInt q=0; q<FEQuadrature::Nq; q++) {

        // Shorthand for values and derivatives of the solution at the single quadrature point.

        PetscReal c_qq = c_q[q];

        const PetscReal    jw  = JxW[q];
        
        for (PetscInt k=0; k<4; k++) {   // Test functions
          for (PetscInt l=0; l<4; l++) { // Trial functions

            const FEFunctionGerm &test_qk=test[q][k];
            const FEFunctionGerm &test_ql=test[q][l];

            // // u-u coupling
            K[k*16+l*2]     += jw*((test_qk.dx*test_ql.dx + test_qk.dy*test_ql.dy) + c_qq*test_qk.val*test_ql.val);
            // // v-v coupling
            K[k*16+8+l*2+1] += jw*((test_qk.dx*test_ql.dx + test_qk.dy*test_ql.dy) + c_qq*test_qk.val*test_ql.val);
          } // l
        } // k
      } // q
      ierr = m_dofmap.addLocalJacobianBlock(K,J);
    } // j
  } // i

  if(dirichletBC) {
    ierr = dirichletBC.fixJacobian2V(J); CHKERRQ(ierr);
  }
  ierr = dirichletBC.finish(); CHKERRQ(ierr);

  ierr = MatAssemblyBegin(J,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(J,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

  ierr = MatSetOption(J,MAT_NEW_NONZERO_LOCATION_ERR,PETSC_TRUE);CHKERRQ(ierr);

  ierr = m_c->end_access(); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode InvSchrodTikhonov::evalObjective(IceModelVec2S &dc, PetscReal *OUTPUT) {

  PetscErrorCode ierr;
  ierr = m_designFunctional->valueAt(dc,OUTPUT);
  return 0;
}

PetscErrorCode InvSchrodTikhonov::evalGradObjective(IceModelVec2S &dc, IceModelVec2S &gradient) {
  PetscErrorCode ierr;
  ierr = m_designFunctional->gradientAt(dc,gradient);
  return 0;
}

PetscErrorCode InvSchrodTikhonov::evalPenalty(IceModelVec2V &du, PetscReal *OUTPUT) {
  PetscErrorCode ierr;
  ierr = m_penaltyFunctional->valueAt(du,OUTPUT);
  return 0;
}

PetscErrorCode InvSchrodTikhonov::evalGradPenaltyReduced(IceModelVec2V &du, IceModelVec2S &gradient) {
  PetscErrorCode ierr;

  // Clear the gradient in prep for updating it.
  gradient.set(0);

  // Assemble the Jacobian matrix.
  PISMVector2 **u_a;
  ierr = m_u.get_array(u_a); CHKERRQ(ierr);  

  DMDALocalInfo *info = NULL;
  this->assembleJacobian( info, u_a, m_J);

  m_penaltyFunctional->gradientAt(du,m_adjointRHS);
  if(m_dirichletLocations) {
    // Start access of Dirichlet data, if present.
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
  ierr = KSPSetOperators(m_ksp, m_J, m_J, SAME_NONZERO_PATTERN); CHKERRQ(ierr);
  ierr = KSPSolve(m_ksp, m_adjointRHS.get_vec(), m_vGlobal.get_vec()); CHKERRQ(ierr); // SOLVE

  ierr = KSPGetConvergedReason(m_ksp, &kspreason); CHKERRQ(ierr);
  
  if (kspreason < 0) {
    SETERRQ1(m_grid.com,1,"InvSchrodTikhonov adjoint linear solve failed (KSP reason %s)",KSPConvergedReasons[kspreason]);
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

  // Loop through all LOCAL elements.
  PetscInt xs = m_element_index.lxs, xm = m_element_index.lxm,
           ys = m_element_index.lys, ym = m_element_index.lym;
  for (i=xs; i<xs+xm; i++) {
    for (j=ys; j<ys+ym; j++) {

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
        for (PetscInt k=0; k<FEQuadrature::Nk; k++) {
          gradient_e[k] += -JxW[q]*(v_qq.u*u_qq.u+v_qq.v*u_qq.v)*test[q][k].val;
        }
      } // q
      m_dofmap.addLocalResidualBlock(gradient_e,gradient_a);
    } // j
  } // i

  if(m_fixedDesignLocations) {
    DirichletData dirichletBC;
    ierr = dirichletBC.init(m_fixedDesignLocations); CHKERRQ(ierr);
    dirichletBC.fixResidualHomogeneous(gradient_a);
    ierr = dirichletBC.finish(); CHKERRQ(ierr);
  }

  ierr = m_v.end_access(); CHKERRQ(ierr);
  ierr = m_u.end_access(); CHKERRQ(ierr);
  ierr = gradient.end_access(); CHKERRQ(ierr);
  return 0;

}
