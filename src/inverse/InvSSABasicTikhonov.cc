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

#include "InvSSABasicTikhonov.hh"
#include "H1NormFunctional.hh"
#include "MeanSquareFunctional.hh"

InvSSABasicTikhonov::InvSSABasicTikhonov( IceGrid  &grid, IceModelVec2V &f, PetscReal p, PetscReal q,  InvTaucParameterization &tp) :
m_grid(grid), m_zeta(NULL), m_f(&f), m_p(p), m_q(q), m_B(1), m_H(1),
m_epsilon_velocity(1e-3), m_epsilon_strainrate(1e-3), m_epsilon_nuH(1e-3), m_pseudo_plastic_threshold(1e-1),
m_dirichletLocations(NULL), m_dirichletValues(NULL), m_dirichletWeight(1.), 
m_fixedDesignLocations(NULL), m_observationWeights(NULL), 
m_tauc_param(tp), m_element_index(grid) {
  PetscErrorCode ierr = this->construct();
  CHKERRCONTINUE(ierr);
  if(ierr) {
    PetscPrintf(m_grid.com, "FATAL ERROR: InvSSABasicTikhonov allocation failed.\n");
    PISMEnd();    
  }
}

InvSSABasicTikhonov::~InvSSABasicTikhonov( ) {
  PetscErrorCode ierr = this->destruct();
  CHKERRCONTINUE(ierr);
  if(ierr) {
    PetscPrintf(m_grid.com, "FATAL ERROR: InvSSABasicTikhonov deallocation failed.\n");
    PISMEnd();    
  }
}

PetscErrorCode InvSSABasicTikhonov::construct() {
  PetscErrorCode ierr;
  PetscInt stencilWidth = 1;

  m_c.create(m_grid,"pseudo yield stress",kHasGhosts,stencilWidth);

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

PetscErrorCode InvSSABasicTikhonov::solve(bool &success) {
  PetscErrorCode ierr;

  PetscReal cL2 = m_grid.config.get("inv_ssa_cL2");
  PetscReal cH1 = m_grid.config.get("inv_ssa_cH1");

  m_designFunctional.reset(new H1NormFunctional2S(m_grid,cL2,cH1,m_fixedDesignLocations));    

  m_penaltyFunctional.reset(new MeanSquareFunctional2V(m_grid,m_observationWeights));    
  PetscReal velocity_scale = m_grid.config.get("inv_ssa_velocity_scale");  
  (reinterpret_cast<MeanSquareFunctional2V&>(*m_penaltyFunctional)).normalize(velocity_scale);

  success=false;
  ierr = SNESSolve(m_snes,NULL,m_uGlobal.get_vec()); CHKERRQ(ierr);

  ierr = SNESGetConvergedReason( m_snes, &m_reason); CHKERRQ(ierr);
  if(m_reason>0){
    success = true;
  } else {

    KSP ksp;
    ierr = SNESGetKSP(m_snes,&ksp); CHKERRQ(ierr);

    KSPConvergedReason kspreason;
    ierr = KSPGetConvergedReason(ksp,&kspreason); CHKERRQ(ierr);
    PetscPrintf(PETSC_COMM_WORLD,"snes converged reason %d\n",m_reason);
    if(m_reason == SNES_DIVERGED_LINEAR_SOLVE) {
      PetscPrintf(PETSC_COMM_WORLD,"ksp converged reason %d\n",kspreason);

      /*
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
      ierr = PetscViewerASCIIOpen(m_grid.com,"InvSSABasicTikhonovFailed.m",&viewer); CHKERRQ(ierr);
      ierr = PetscViewerSetFormat(viewer, PETSC_VIEWER_ASCII_MATLAB); CHKERRQ(ierr);
      Mat J, Jpc;
      MatStructure ms;
      ierr = KSPGetOperators(ksp,&J,&Jpc,&ms);
      ierr = MatView(J,viewer); CHKERRQ(ierr);
      ierr = MatView(Jpc,viewer); CHKERRQ(ierr);
      ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);
      */
    }
  }
  m_u.copy_from(m_uGlobal);
  return 0;
}

PetscErrorCode InvSSABasicTikhonov::destruct() {
  PetscErrorCode ierr;
  ierr = DMDestroy(&m_da); CHKERRQ(ierr);
  ierr = MatDestroy(&m_J); CHKERRQ(ierr);
  ierr = SNESDestroy(&m_snes); CHKERRQ(ierr);
  ierr = KSPDestroy(&m_ksp); CHKERRQ(ierr);
  return 0;
}

void InvSSABasicTikhonov::setParamsFromConfig( NCConfigVariable &config ) {
  m_epsilon_velocity         = config.get("plastic_regularization") / secpera;

  PetscReal schoofLen = config.get("Schoof_regularizing_length", "km", "m"); // convert to meters
  PetscReal schoofVel = config.get("Schoof_regularizing_velocity", "m/year", "m/s"); // convert to m/s
  m_epsilon_strainrate = schoofVel/schoofLen;

  m_epsilon_nuH = config.get("epsilon_ssa");

  m_pseudo_plastic_threshold = config.get("pseudo_plastic_uthreshold") / secpera;
}

PetscErrorCode InvSSABasicTikhonov::setZeta(IceModelVec2S &zeta) {

  PetscErrorCode ierr;
  
  m_zeta = &zeta;

  PetscReal **zeta_a;
  PetscReal **c_a;
  ierr = m_zeta->get_array(zeta_a); CHKERRQ(ierr);
  ierr = m_c.get_array(c_a); CHKERRQ(ierr);
  for(PetscInt i=m_grid.xs;i<m_grid.xs+m_grid.xm;i++){
    for(PetscInt j=m_grid.ys;j<m_grid.ys+m_grid.ym;j++){
      PetscReal tauc;
      m_tauc_param.toTauc(zeta_a[i][j],&tauc,NULL);
      c_a[i][j] = tauc;
    }
  }
  ierr = m_c.end_access(); CHKERRQ(ierr);
  ierr = m_zeta->end_access(); CHKERRQ(ierr);
  ierr = m_c.beginGhostComm(); CHKERRQ(ierr);
  ierr = m_c.endGhostComm(); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode InvSSABasicTikhonov::linearizeAt( IceModelVec2S &zeta, bool &success) {

  PetscErrorCode ierr;
  ierr = this->setZeta(zeta); CHKERRQ(ierr);
  ierr = this->solve(success); CHKERRQ(ierr);
  return 0;

}

void InvSSABasicTikhonov::computeNuH(PetscReal B, PetscReal H, const PetscReal Du[],
  PetscReal *nuH, PetscReal *dNuH) {
  const PetscReal alpha =  PetscSqr(m_epsilon_strainrate) + 0.5 * (PetscSqr(Du[0]) + PetscSqr(Du[1]) + PetscSqr(Du[0]+Du[1]) + 2*PetscSqr(Du[2]));
  const PetscReal power = (m_p-2)*0.5;
  
  // Strictly speaking, nuH has a factor of .5, and the SSA has coefficents
  // of the form 2*nuH.  So we just cancel these dualing factors here
  *nuH = B * H * pow(alpha, power);
  *nuH += 2*m_epsilon_nuH;
  
  if(dNuH) *dNuH = power * (*nuH) / alpha;

}


void InvSSABasicTikhonov::computeBeta(PetscReal tauc, PISMVector2 U,PetscReal *beta, PetscReal *dbeta) {
  const PetscScalar Usq = PetscSqr(m_epsilon_velocity) + PetscSqr(U.u) + PetscSqr(U.v);
  *beta = tauc * pow(Usq, 0.5*(m_q - 1)) * pow(m_pseudo_plastic_threshold, -m_q);
  if (dbeta) *dbeta = (m_q - 1) * *beta / Usq;
}


PetscErrorCode InvSSABasicTikhonov::assembleFunction( DMDALocalInfo *info, PISMVector2 **u_a, PISMVector2 **y_a) {
  PetscInt         i,j,k,q;
  PetscErrorCode   ierr;

  (void) info; // Avoid compiler warning.

  // Zero out the portion of the function we are responsible for computing.
  for (i=m_grid.xs; i<m_grid.xs+m_grid.xm; i++) {
    for (j=m_grid.ys; j<m_grid.ys+m_grid.ym; j++) {
      y_a[i][j].u = y_a[i][j].v = 0;
    }
  }

  // Start access of Dirichlet data, if present.
  DirichletData dirichletBC;
  ierr = dirichletBC.init(m_dirichletLocations,m_dirichletValues,m_dirichletWeight); CHKERRQ(ierr);
  // Jacobian times weights for quadrature.
  PetscScalar JxW[FEQuadrature::Nq];
  m_quadrature.getWeightedJacobian(JxW);

  // Storage for the current solution at quadrature points.
  PISMVector2 u_e[FEQuadrature::Nk];
  PISMVector2 u_q[FEQuadrature::Nq];
  PetscReal  Du_q[FEQuadrature::Nq][3];

  PISMVector2 **f_a;
  ierr = m_f->get_array(f_a); CHKERRQ(ierr);
  PISMVector2 f_q[FEQuadrature::Nq];

  PetscReal **c_a;
  ierr =      m_c.get_array(c_a); CHKERRQ(ierr);
  PetscReal   c_q[FEQuadrature::Nq];

  PISMVector2   y_e[FEQuadrature::Nk];

  // An Nq by Nk array of test function values.
  const FEFunctionGerm (*test)[FEQuadrature::Nk] = m_quadrature.testFunctionValues();

  PetscReal testx_DD[FEQuadrature::Nq][FEQuadrature::Nk][4];
  PetscReal testy_DD[FEQuadrature::Nq][FEQuadrature::Nk][4];
  for(q=0; q<FEQuadrature::Nq;q++) {
    for(k=0;k<FEQuadrature::Nk;k++) {
      testx_DD[q][k][0] = test[q][k].dx;     testy_DD[q][k][0] = 0;
      testx_DD[q][k][1] = 0;                 testy_DD[q][k][1] = test[q][k].dy;
      testx_DD[q][k][2] = 0.5*test[q][k].dy; testy_DD[q][k][2] = 0.5*test[q][k].dx;
      testx_DD[q][k][3] = test[q][k].dx;     testy_DD[q][k][3] = test[q][k].dy;
    }
  }
  
  
  // Iterate over the elements.
  PetscInt xs = m_element_index.xs, xm = m_element_index.xm,
           ys = m_element_index.ys, ym = m_element_index.ym;
  for (i=xs; i<xs+xm; i++) {
    for (j=ys; j<ys+ym; j++) {

      // Zero out the element-local residual in prep for updating it.
      for(k=0;k<FEQuadrature::Nk;k++){
        y_e[k].u = 0; y_e[k].v = 0;
      }

      // Initialize the map from global to local degrees of freedom for this element.
      m_dofmap.reset(i,j,m_grid);

      // Obtain the value of the solution at the nodes adjacent to the element.
      m_dofmap.extractLocalDOFs(i,j,u_a,u_e);
      if(dirichletBC) dirichletBC.update(m_dofmap,u_e);
      m_quadrature.computeTrialFunctionValues(u_e,u_q,Du_q);

      m_quadrature.computeTrialFunctionValues(i,j,m_dofmap,f_a,f_q);
      m_quadrature.computeTrialFunctionValues(i,j,m_dofmap,c_a,c_q);

      for (q=0; q<FEQuadrature::Nq; q++) {     // loop over quadrature points on this element.

        const PetscReal    jw  = JxW[q];
        PetscReal beta;
        this->computeBeta(c_q[q],u_q[q],&beta, NULL);

        PetscReal *Du_qq = Du_q[q];
        PetscReal nuH;
        computeNuH(m_B,m_H,Du_qq,&nuH,NULL);

        PetscReal Fu_q =beta*u_q[q].u-f_q[q].u;
        PetscReal Fv_q =beta*u_q[q].v-f_q[q].v;
        for(k=0; k<4;k++) {  // loop over the test functions.
          const FEFunctionGerm &testqk = test[q][k];

          const PetscReal *tDDx = testx_DD[q][k];
          const PetscReal *tDDy = testy_DD[q][k];

          y_e[k].u += jw*(nuH*(tDDx[0]*Du_qq[0]+tDDx[1]*Du_qq[1]+2*tDDx[2]*Du_qq[2]+ tDDx[3]*(Du_qq[0]+Du_qq[1])) +
            testqk.val*Fu_q);
          y_e[k].v += jw*(nuH*(tDDy[0]*Du_qq[0]+tDDy[1]*Du_qq[1]+2*tDDy[2]*Du_qq[2]+ tDDy[3]*(Du_qq[0]+Du_qq[1])) +
            testqk.val*Fv_q);
        }
      } // q
      m_dofmap.addLocalResidualBlock(y_e,y_a);
    } // j-loop
  } // i-loop

  // Until now we have not touched rows in the residual corresponding to Dirichlet data.
  // We fix this now.
  if(dirichletBC) dirichletBC.fixResidual(u_a,y_a);
  ierr = dirichletBC.finish(); CHKERRQ(ierr);

  ierr = m_f->end_access(); CHKERRQ(ierr);
  ierr = m_c.end_access(); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode InvSSABasicTikhonov::assembleJacobian( DMDALocalInfo *info, PISMVector2 **u_a, Mat J) {
  PetscInt         i,j;
  PetscErrorCode   ierr;

  // Avoid compiler warning.
  (void) info;

  PetscReal      **c_a;
  PetscReal      c_q[FEQuadrature::Nq];
  ierr = m_c.get_array(c_a); CHKERRQ(ierr);

  // Storage for the current solution at quadrature points.
  PISMVector2    u_q[FEQuadrature::Nq];
  PISMVector2    u_e[FEQuadrature::Nk];
  PetscReal     Du_q[FEQuadrature::Nq][3];
  
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

  PetscReal testx_DD[FEQuadrature::Nq][FEQuadrature::Nk][4];
  PetscReal testy_DD[FEQuadrature::Nq][FEQuadrature::Nk][4];
  for(PetscInt q=0; q<FEQuadrature::Nq;q++) {
    for(PetscInt k=0;k<FEQuadrature::Nk;k++) {
      testx_DD[q][k][0] = test[q][k].dx;     testy_DD[q][k][0] = 0;
      testx_DD[q][k][1] = 0;                 testy_DD[q][k][1] = test[q][k].dy;
      testx_DD[q][k][2] = 0.5*test[q][k].dy; testy_DD[q][k][2] = 0.5*test[q][k].dx;
      testx_DD[q][k][3] = test[q][k].dx;     testy_DD[q][k][3] = test[q][k].dy;
    }
  }

  // Loop through all the elements.
  PetscInt xs = m_element_index.xs, xm = m_element_index.xm,
           ys = m_element_index.ys, ym = m_element_index.ym;
  for (i=xs; i<xs+xm; i++) {
    for (j=ys; j<ys+ym; j++) {

      // Element-local Jacobian matrix (there are FEQuadrature::Nk vector valued degrees
      // of freedom per elment, for a total of (2*FEQuadrature::Nk)*(2*FEQuadrature::Nk) = 16
      // entries in the local Jacobian.
      PetscReal      K[(2*FEQuadrature::Nk)*(2*FEQuadrature::Nk)];

      // Initialize the map from global to local degrees of freedom for this element.
      m_dofmap.reset(i,j,m_grid);

      // Obtain the value of the solution at the adjacent nodes to the element.
      m_dofmap.extractLocalDOFs(i,j, u_a,u_e);
      if(dirichletBC) dirichletBC.update(m_dofmap,u_e);
      m_quadrature.computeTrialFunctionValues(u_e,u_q,Du_q);

      // Compute the values of the solution at the quadrature points.
      m_quadrature.computeTrialFunctionValues(i,j,m_dofmap,c_a,c_q);

      // Build the element-local Jacobian.
      ierr = PetscMemzero(K,sizeof(K));CHKERRQ(ierr);
      for (PetscInt q=0; q<FEQuadrature::Nq; q++) {

        // Shorthand for values and derivatives of the solution at the single quadrature point.

        const PISMVector2 &u_qq = u_q[q];

        PetscReal beta, dbeta;
        this->computeBeta(c_q[q],u_qq,&beta, &dbeta);

        PetscReal nuH, dNuH;
        PetscReal *Du_qq = Du_q[q];
        this->computeNuH(m_B,m_H,Du_qq,&nuH,&dNuH);

        const PetscReal    jw  = JxW[q];

        PetscReal Dtestx_dot_Du[FEQuadrature::Nk];
        PetscReal Dtesty_dot_Du[FEQuadrature::Nk];
        for(PetscInt k=0;k<FEQuadrature::Nk;k++) {
          PetscReal *Dtestxk=testx_DD[q][k];
          PetscReal *Dtestyk=testy_DD[q][k];
          Dtestx_dot_Du[k] = (Dtestxk[0]*Du_qq[0]+Dtestxk[1]*Du_qq[1]+2*Dtestxk[2]*Du_qq[2]+ Dtestxk[3]*(Du_qq[0]+Du_qq[1]));
          Dtesty_dot_Du[k] = (Dtestyk[0]*Du_qq[0]+Dtestyk[1]*Du_qq[1]+2*Dtestyk[2]*Du_qq[2]+ Dtestyk[3]*(Du_qq[0]+Du_qq[1]));
        }
        
        for (PetscInt k=0; k<4; k++) {   // Test functions
          for (PetscInt l=0; l<4; l++) { // Trial functions

            const FEFunctionGerm &test_qk=test[q][k];
            const FEFunctionGerm &test_ql=test[q][l];

            PetscReal *DDkx = testx_DD[q][k];
            PetscReal *DDlx = testx_DD[q][l];
            PetscReal *DDky = testy_DD[q][k];
            PetscReal *DDly = testy_DD[q][l];
            
            PetscReal Dtestkx_dot_Dtestlx = (DDkx[0]*DDlx[0]+DDkx[1]*DDlx[1]+2*DDkx[2]*DDlx[2]+DDkx[3]*DDlx[3]);
            PetscReal Dtestky_dot_Dtestlx = (DDky[0]*DDlx[0]+DDky[1]*DDlx[1]+2*DDky[2]*DDlx[2]+DDky[3]*DDlx[3]);
            PetscReal Dtestkx_dot_Dtestly = (DDkx[0]*DDly[0]+DDkx[1]*DDly[1]+2*DDkx[2]*DDly[2]+DDkx[3]*DDly[3]);
            PetscReal Dtestky_dot_Dtestly = (DDky[0]*DDly[0]+DDky[1]*DDly[1]+2*DDky[2]*DDly[2]+DDky[3]*DDly[3]);
            PetscReal Dtestkx_dot_Du = Dtestx_dot_Du[k];
            PetscReal Dtestlx_dot_Du = Dtestx_dot_Du[l];
            PetscReal Dtestky_dot_Du = Dtesty_dot_Du[k];
            PetscReal Dtestly_dot_Du = Dtesty_dot_Du[l];

            K[k*16+l*2]     += jw*( nuH*Dtestkx_dot_Dtestlx + dNuH* Dtestlx_dot_Du * Dtestkx_dot_Du +
              beta*test_qk.val*test_ql.val + dbeta*(test_qk.val*u_qq.u*test_ql.val*u_qq.u) );
            // u-v coupling
            K[k*16+l*2+1]   += jw*( nuH*Dtestkx_dot_Dtestly + dNuH* Dtestly_dot_Du * Dtestkx_dot_Du +
              dbeta*(test_qk.val*u_qq.u*test_ql.val*u_qq.v) );
            // v-u coupling
            K[k*16+8+l*2]   += jw*( nuH*Dtestky_dot_Dtestlx + dNuH* Dtestlx_dot_Du * Dtestky_dot_Du +
              dbeta*(test_qk.val*u_qq.v*test_ql.val*u_qq.u) );
            // v-v coupling
            K[k*16+8+l*2+1] += jw*( nuH*Dtestky_dot_Dtestly + dNuH* Dtestly_dot_Du * Dtestky_dot_Du +
              beta*test_qk.val*test_ql.val + dbeta*(test_qk.val*u_qq.v*test_ql.val*u_qq.v));
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

  ierr = m_c.end_access(); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode InvSSABasicTikhonov::evalObjective(IceModelVec2S &dzeta, PetscReal *OUTPUT) {

  PetscErrorCode ierr;
  ierr = m_designFunctional->valueAt(dzeta,OUTPUT);
  return 0;
}

PetscErrorCode InvSSABasicTikhonov::evalGradObjective(IceModelVec2S &dzeta, IceModelVec2S &gradient) {
  PetscErrorCode ierr;
  ierr = m_designFunctional->gradientAt(dzeta,gradient);
  return 0;
}

PetscErrorCode InvSSABasicTikhonov::evalPenalty(IceModelVec2V &du, PetscReal *OUTPUT) {

  PetscErrorCode ierr;
  ierr = m_penaltyFunctional->valueAt(du,OUTPUT);
  return 0;
}

PetscErrorCode InvSSABasicTikhonov::evalGradPenalty(IceModelVec2V &du, IceModelVec2V &gradient) {
  PetscErrorCode ierr;
  ierr = m_penaltyFunctional->gradientAt(du,gradient); CHKERRQ(ierr);
  return 0;
}

PetscErrorCode InvSSABasicTikhonov::evalGradPenaltyReduced(IceModelVec2V &du, IceModelVec2S &gradient) {
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
    SETERRQ1(m_grid.com,1,"InvSSABasicTikhonov adjoint linear solve failed (KSP reason %s)",KSPConvergedReasons[kspreason]);
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
  PetscInt xs = m_element_index.xs, xm = m_element_index.xm,
           ys = m_element_index.ys, ym = m_element_index.ym;
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

        PetscReal dbeta_dc;
        this->computeBeta(1,u_qq,&dbeta_dc,NULL);

        for (PetscInt k=0; k<FEQuadrature::Nk; k++) {
          gradient_e[k] += -JxW[q]*dbeta_dc*(v_qq.u*u_qq.u+v_qq.v*u_qq.v)*test[q][k].val;
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
