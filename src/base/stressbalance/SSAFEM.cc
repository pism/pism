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

#include "SSAFEM.hh"
#include "FETools.hh"


//! \brief Allocating SSAFEM-specific objects; called by the constructor.
PetscErrorCode SSAFEM::allocate_fem() {
  PetscErrorCode ierr;

  dirichletScale = 1.0;
  ocean_rho = config.get("sea_water_density");
  earth_grav = config.get("standard_gravity");

  ierr = DACreateGlobalVector(SSADA, &r);CHKERRQ(ierr);
  ierr = DAGetMatrix(SSADA, "baij", &J); CHKERRQ(ierr);

  ierr = SNESCreate(grid.com,&snes);CHKERRQ(ierr);
  // ierr = SNESSetOptionsPrefix(snes,((PetscObject)this)->prefix);CHKERRQ(ierr);

  // Set the SNES callbacks to call into our compute_local_function and compute_local_jacobian
  // methods via SSAFEFunction and SSAFEJ
  ierr = DASetLocalFunction(SSADA,(DALocalFunction1)SSAFEFunction);CHKERRQ(ierr);
  ierr = DASetLocalJacobian(SSADA,(DALocalFunction1)SSAFEJacobian);CHKERRQ(ierr);
  ctx.da = SSADA;  ctx.ssa = this;
  ierr = SNESSetFunction(snes, r,    SNESDAFormFunction,   &ctx);CHKERRQ(ierr);
  ierr = SNESSetJacobian(snes, J, J, SNESDAComputeJacobian,&ctx);CHKERRQ(ierr);


  ierr = SNESSetFromOptions(snes);CHKERRQ(ierr);
  
  PetscInt nElements = element_index.element_count();
  feStore = new FEStoreNode[FEQuadrature::Nq*nElements];

  // hardav IceModelVec2S is not used (so far).
  const PetscScalar power = 1.0 / ice.exponent();
  char unitstr[TEMPORARY_STRING_LENGTH];
  snprintf(unitstr, sizeof(unitstr), "Pa s%f", power);
  ierr = hardav.create(grid, "hardav", true); CHKERRQ(ierr);
  ierr = hardav.set_attrs("internal", "vertically-averaged ice hardness", unitstr, ""); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode SSAFEM::deallocate_fem() {
  PetscErrorCode ierr;

  ierr = SNESDestroy(snes);CHKERRQ(ierr);
  delete feStore;
  ierr = VecDestroy(r); CHKERRQ(ierr);
  ierr = MatDestroy(J); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode SSAFEM::init(PISMVars &vars) {
  PetscErrorCode ierr;

  ierr = SSA::init(vars); CHKERRQ(ierr);
  ierr = verbPrintf(2,grid.com,
           "  [using the SNES-based finite element method implementation]\n");
           CHKERRQ(ierr);

  ierr = setFromOptions(); CHKERRQ(ierr);

  // On restart, SSA::init() reads the SSA velocity from a PISM output file
  // into IceModelVec2V "velocity". We use that field as an initial guess.
  // If we are not restarting from a PISM file, "velocity" is identically zero,
  // and the call below clears SSAX.

  ierr = velocity.copy_to(SSAX); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode SSAFEM::setFromOptions()
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = PetscOptionsHead("SSA FEM options");CHKERRQ(ierr);
  ierr = PetscOptionsReal("-ssa_fe_dirichlet_scale",
                          "Enforce Dirichlet conditions with this additional scaling",
                          "",
                          dirichletScale,
                          &dirichletScale,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsTail();CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode SSAFEM::solve()
{
  PetscErrorCode ierr;
  PetscViewer    viewer;
  char           filename[PETSC_MAX_PATH_LEN];
  PetscTruth     flg;

  PetscFunctionBegin;
  ierr = PetscOptionsGetString(NULL, "-ssa_view", filename,
                               PETSC_MAX_PATH_LEN, &flg); CHKERRQ(ierr);
  if (flg) {
    ierr = PetscViewerASCIIOpen(grid.com,filename,&viewer);
             CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(viewer,"SNES before SSASolve_FE\n");
             CHKERRQ(ierr);
    ierr = SNESView(snes,viewer);CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(viewer,"solution vector before SSASolve_FE\n");
             CHKERRQ(ierr);
    ierr = VecView(SSAX,viewer);CHKERRQ(ierr);
    ierr = PetscViewerDestroy(viewer);CHKERRQ(ierr);
  }

  stdout_ssa = "";
  if (getVerbosityLevel() >= 2) 
    stdout_ssa = "  SSA: ";
  
  // Set up the system to solve:
  ierr = setup(); CHKERRQ(ierr);

  // Solve:
  ierr = SNESSolve(snes,NULL,SSAX);CHKERRQ(ierr);

  // See if it worked.
  SNESConvergedReason reason;
  ierr = SNESGetConvergedReason( snes, &reason); CHKERRQ(ierr);
  if(reason < 0)
  {
    SETERRQ1(1, 
      "SSAFEM solve failed to converge (SNES reason %s)\n\n", SNESConvergedReasons[reason]);
  }
  else if(getVerbosityLevel() > 2)
  {
    stdout_ssa += "SSAFEM converged (SNES reason ";
    stdout_ssa += SNESConvergedReasons[reason];
    stdout_ssa += ")\n";
  }

  // Extract the solution back from SSAX to velocity and communicate.
  ierr = velocity.copy_from(SSAX); CHKERRQ(ierr);
  ierr = velocity.beginGhostComm(); CHKERRQ(ierr);
  ierr = velocity.endGhostComm(); CHKERRQ(ierr);

  ierr = PetscOptionsHasName(NULL,"-ssa_view_solution",&flg);CHKERRQ(ierr);
  if (flg) {
    ierr = PetscViewerASCIIOpen(grid.com,filename,&viewer);CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(viewer,"solution vector after SSASolve\n");
             CHKERRQ(ierr);
    ierr = VecView(SSAX,viewer);CHKERRQ(ierr);
    ierr = PetscViewerDestroy(viewer);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

// This is called when surface and/or temperature have changed and sets up the
// system solved using SNES.
//
// Since this interfaces with the rest of PISM and doesn't touch velocity or
// stresses, everything is fully dimensional
PetscErrorCode SSAFEM::setup()
{
  // DALocalInfo      info_struct,
  //                 *info = &info_struct; // DAGetLocalInfo writes into our memory
  PetscReal      **h,
                 **H,
                 **topg,
                 **tauc_array,
                  *Enth_e[4],
                  *Enth_q[4];
  PetscInt         i,j,k,q,p,
                   Mz = grid.Mz;
  PetscErrorCode   ierr;

  for(q=0;q<FEQuadrature::Nq;q++)
  {
    Enth_q[q] = new PetscReal[Mz];
  }

  ierr = enthalpy->begin_access();CHKERRQ(ierr);
  ierr = surface->get_array(h);CHKERRQ(ierr);
  ierr = thickness->get_array(H);CHKERRQ(ierr);
  ierr = bed->get_array(topg);CHKERRQ(ierr);
  ierr = tauc->get_array(tauc_array);CHKERRQ(ierr);

  PetscInt xs = element_index.xs, xm = element_index.xm,
           ys = element_index.ys, ym = element_index.ym;  
  for (i=xs; i<xs+xm; i++) {
    for (j=ys;j<ys+ym; j++) {
      
      // Extract coefficient values at the quadrature points.
      PetscReal hq[FEQuadrature::Nq],hxq[FEQuadrature::Nq],hyq[FEQuadrature::Nq];
      quadrature.computeTrialFunctionValues(i,j,dofmap,h,hq,hxq,hyq);

      PetscReal Hq[FEQuadrature::Nq], bq[FEQuadrature::Nq], taucq[FEQuadrature::Nq];
      quadrature.computeTrialFunctionValues(i,j,dofmap,H,Hq);
      quadrature.computeTrialFunctionValues(i,j,dofmap,topg,bq);
      quadrature.computeTrialFunctionValues(i,j,dofmap,tauc_array,taucq);

      const PetscInt ij = element_index.flatten(i,j);
      FEStoreNode *feS = &feStore[4*ij];
      for (q=0; q<4; q++) {
        feS[q].h  = hq[q];
        feS[q].H  = Hq[q];
        feS[q].b  = bq[q];
        feS[q].tauc = taucq[q];
        feS[q].hx = hxq[q];
        feS[q].hy = hyq[q];
      }

      // Surface and thickness information is stored, now do the thermal stuff
      ierr = enthalpy->getInternalColumn(i,j,&Enth_e[0]);
      ierr = enthalpy->getInternalColumn(i+1,j,&Enth_e[1]);CHKERRQ(ierr);
      ierr = enthalpy->getInternalColumn(i+1,j+1,&Enth_e[2]);CHKERRQ(ierr);
      ierr = enthalpy->getInternalColumn(i,j+1,&Enth_e[3]);CHKERRQ(ierr);
      // Interpolate to quadrature points at every vertical level
      const FEFunctionGerm (*test)[FEQuadrature::Nk] = quadrature.testFunctionValues();
      for (k=0; k<Mz; k++) { // This loop could be cut short at the surface.
        Enth_q[0][k] = Enth_q[1][k] = Enth_q[2][k] = Enth_q[3][k] = 0;
        for (q=0; q<FEQuadrature::Nq; q++) {
          for (p=0; p<FEQuadrature::Nk; p++) {
            Enth_q[q][k] += test[q][k].val * Enth_e[p][k];
          }
        }
      }
      for (q=0; q<4; q++) {
        // Evaluate column integrals in flow law at every quadrature point's column
        feS[q].B = ice.averagedHardness_from_enth(feS[q].H, grid.kBelowHeight(feS[q].H),
                                                  grid.zlevels, Enth_q[q]);
      }
    }
  }
  ierr = surface->end_access();CHKERRQ(ierr);
  ierr = thickness->end_access();CHKERRQ(ierr);
  ierr = bed->end_access();CHKERRQ(ierr);
  ierr = tauc->end_access();CHKERRQ(ierr);
  ierr = enthalpy->end_access();CHKERRQ(ierr);

  for(q=0;q<4;q++)
  {
    delete Enth_q[q];
  }
  
  return 0;
}

inline PetscErrorCode SSAFEM::PointwiseNuHAndBeta(const FEStoreNode *feS,
                                                  const PISMVector2 *u,const PetscReal Du[],
                                                  PetscReal *nuH, PetscReal *dNuH,
                                                  PetscReal *beta, PetscReal *dbeta)
{
  if (feS->H < strength_extension->get_min_thickness()) {
    *nuH = strength_extension->get_notional_strength();
    if (dNuH) *dNuH = 0;
  } else {
    // PetscReal dimDu[3];
    // for (int i=0; i<3; i++) dimDu[i] = ref.StrainRate() * Du[i];
    ice.effectiveViscosity_with_derivative(feS->B, Du, nuH, dNuH);
    *nuH  *= feS->H;
    if (dNuH) *dNuH *= feS->H;
  }
  *nuH  *=  2; // / ref.IntegratedViscosity();
  if (dNuH) *dNuH *= 2; // * ref.StrainRate2() / ref.IntegratedViscosity(); // The derivative with respect to second invariant
  
  if (Floating(ice,ocean_rho,feS->H,feS->b)) {
    // The ice is floating here so there is no friction. Note that the purpose
    // of checking flotation this way is to get subgrid resolution of stress in
    // the vicinity of the grounding line. According to Goldberg et. al. 2009
    // (probably will be published in 2009...) this is important to loosen the
    // resolution requirements near the grounding line.
    *beta = 0;
    if (dbeta) *dbeta = 0;
    // SETERRQ(1,"Not tested yet");
  } else {
    basal.dragWithDerivative(feS->tauc,u->u,u->v,beta,dbeta);
  }
  // *beta /= ref.Drag();
  // if (dbeta) *dbeta *= ref.Velocity2() / ref.Drag();
 
  return 0;
}

//! \brief Sets Dirichlet boundary conditions. Called from SSAFEFunction and
//! SSAFEJacobian.
/*! The arrays lmask, row, col, and x are all of length N, the number of
nodes in an element.  If for some node lmask indicates that it is a 
Dirichlet node, the values of x from the node is set from the Dirichlet
data BC_vel, and the row and column stencils are made invalid.  Calls
to functions such as MatSetStencilBlock will then ignore Dirichlet rows
and columns.
*/
void SSAFEM::FixDirichletValues(PetscReal lmask[],PISMVector2 **BC_vel,
                                PISMVector2 x[], FEDOFMap &dofmap)
{
  for (PetscInt k=0; k<4; k++) {
    if (PismIntMask(lmask[k]) == MASK_SHEET) {
      PetscInt ii, jj;
      dofmap.localToGlobal(k,&ii,&jj);
      x[k].u = BC_vel[ii][jj].u;
      x[k].v = BC_vel[ii][jj].v;
      dofmap.markRowInvalid(k);
      dofmap.markColInvalid(k);
    }
  }
}


PetscErrorCode SSAFEM::compute_local_function(DALocalInfo *info, const PISMVector2 **xg, PISMVector2 **yg)
{
  PetscInt         i,j,k,q;
  PetscReal        **mask;
  PISMVector2        **BC_vel;
  PetscErrorCode   ierr;

  // Zero out the portion of the function we are responsible for computing.
  for (i=grid.xs; i<grid.xs+grid.xm; i++) {
    for (j=grid.ys; j<grid.ys+grid.ym; j++) {
      yg[i][j].u = yg[i][j].v = 0;
    }
  }
  
  if (bc_locations && vel_bc) {
    ierr = bc_locations->get_array(mask);CHKERRQ(ierr);
    ierr = vel_bc->get_array(BC_vel); CHKERRQ(ierr);
  }

  PetscScalar JxW[FEQuadrature::Nq];
  quadrature.getWeightedJacobian(JxW);
  PISMVector2 u[FEQuadrature::Nq];
  PetscScalar Du[FEQuadrature::Nq][3];
  const FEFunctionGerm (*test)[FEQuadrature::Nk] = quadrature.testFunctionValues();

  PetscInt xs = element_index.xs, xm = element_index.xm,
           ys = element_index.ys, ym = element_index.ym;
  for (i=xs; i<xs+xm; i++) {
    for (j=ys; j<ys+ym; j++) {
      const PetscInt ij = element_index.flatten(i,j);
      PISMVector2     x[4],y[4];

      dofmap.extractLocalDOFs(i,j,xg,x);
      dofmap.reset(i,j,grid);

      for(k=0;k<FEQuadrature::Nk;k++){ 
        y[k].u = 0; y[k].v=0;
      }

      if (bc_locations && vel_bc) {
        PetscReal lmask[4];
        dofmap.extractLocalDOFs(i,j,mask,lmask);
        FixDirichletValues(lmask,BC_vel,x,dofmap);
      }

      quadrature.computeTrialFunctionValues(x,u,Du);
      for (q=0; q<FEQuadrature::Nq; q++) {     // loop over quadrature points on this element.
        const FEStoreNode *feS = &feStore[ij*4+q];
        const PetscReal    jw  = JxW[q];

        PetscReal nuH, beta;
        ierr = PointwiseNuHAndBeta(feS,u+q,Du[q],&nuH,NULL,&beta,NULL);CHKERRQ(ierr);

        PetscScalar *Duq = Du[q];

        PISMVector2 f;
        f.u = beta*u[q].u + ice.rho*earth_grav*feS->H*feS->hx;
        f.v = beta*u[q].v + ice.rho*earth_grav*feS->H*feS->hy;
        
        for(k=0; k<4;k++) {  // loop over the test functions.
          const FEFunctionGerm &testqk = test[q][k];
          y[k].u += jw*(nuH*(testqk.dx*(2*Duq[0]+Duq[1]) + testqk.dy*Duq[2]) + testqk.val*f.u);
          y[k].v += jw*(nuH*(testqk.dy*(2*Duq[1]+Duq[0]) + testqk.dx*Duq[2]) + testqk.val*f.v);
        }
      }
      dofmap.addLocalResidualBlock(y,yg);
    } // j-loop
  } // i-loop

  if (bc_locations && vel_bc) {
    // Enforce Dirichlet conditions strongly
    for (i=grid.xs; i<grid.xs+grid.xm; i++) {
      for (j=grid.ys; j<grid.ys+grid.ym; j++) {
        if (bc_locations->value(i,j) == MASK_SHEET) {
          // Enforce zero sliding strongly
          yg[i][j].u = dirichletScale * (xg[i][j].u - BC_vel[i][j].u);
          yg[i][j].v = dirichletScale * (xg[i][j].v - BC_vel[i][j].v);
        }
      }
    }
    ierr = bc_locations->end_access();CHKERRQ(ierr);
    ierr = vel_bc->end_access(); CHKERRQ(ierr);
  }

  PetscTruth monitorFunction;
  ierr = PetscOptionsHasName(NULL,"-ssa_monitor_function",&monitorFunction);CHKERRQ(ierr);
  if (monitorFunction) {
    ierr = PetscPrintf(grid.com,"SSA Solution and Function values (pointwise residuals)\n");CHKERRQ(ierr);
    for (i=grid.xs; i<grid.xs+grid.xm; i++) {
      for (j=grid.ys; j<grid.ys+grid.ym; j++) {
        ierr = PetscSynchronizedPrintf(grid.com,
                                       "[%2d,%2d] u=(%12.4e,%12.4e)  f=(%12.4e,%12.4e)\n",
                                       i,j,xg[i][j].u,xg[i][j].v,yg[i][j].u,yg[i][j].v);CHKERRQ(ierr);
      }
    }
    ierr = PetscSynchronizedFlush(grid.com);CHKERRQ(ierr);
  }
  return 0;
}

PetscErrorCode SSAFEM::compute_local_jacobian(DALocalInfo *info, const PISMVector2 **xg, Mat J )
{
  PetscReal      **mask;
  PISMVector2    **BC_vel;
  PetscInt         i,j;
  PetscErrorCode   ierr;

  ierr = MatZeroEntries(J);CHKERRQ(ierr);

  if (bc_locations && vel_bc) {
    ierr = bc_locations->get_array(mask);CHKERRQ(ierr);
    ierr = vel_bc->get_array(BC_vel); CHKERRQ(ierr); 
  }

  // Jacobian times weights for quadrature.
  PetscScalar JxW[FEQuadrature::Nq];
  quadrature.getWeightedJacobian(JxW);

  // Storage for the solution at quadrature points.
  PISMVector2 w[FEQuadrature::Nq];
  PetscScalar Dw[FEQuadrature::Nq][3];

  // Values of the finite element test functions at the quadrature points.
  // This is an Nq by Nk array of function germs (Nq=#of quad pts, Nk=#of test functions).
  const FEFunctionGerm (*test)[FEQuadrature::Nk] = quadrature.testFunctionValues();

  // Loop through all the elements.
  PetscInt xs = element_index.xs, xm = element_index.xm,
           ys = element_index.ys, ym = element_index.ym;
  for (i=xs; i<xs+xm; i++) {
    for (j=ys; j<ys+ym; j++) {

      const PetscInt ij = element_index.flatten(i,j);
      PISMVector2    x[4];
      PetscReal      K[8*8];
      
      // Initialize the map from global to local degrees of freedom for this element.
      dofmap.reset(i,j,grid);

      // Obtain the value of the solution at the adjacent nodes to the element.
      dofmap.extractLocalDOFs(i,j,xg,x);

      // Adjust these values if needed based on Dirichlet data.
      if (bc_locations && vel_bc) {
        PetscReal lmask[4];
        dofmap.extractLocalDOFs(i,j,mask,lmask);
        FixDirichletValues(lmask,BC_vel,x,dofmap);
      }

      // Compute the values of the solution at the quadrature points.
      quadrature.computeTrialFunctionValues(x,w,Dw);

      // Build the local interaction matrix (K).
      ierr = PetscMemzero(K,sizeof(K));CHKERRQ(ierr);
      for (PetscInt q=0; q<FEQuadrature::Nq; q++) {

        // Shorthand for values and derivatives of the solution at the single quadrature point.
        PISMVector2 &wq = w[q];
        PetscReal *Dwq = Dw[q];

        // Coefficients evaluated at the single quadrature point.
        const FEStoreNode *feS = &feStore[ij*4+q];
        const PetscReal    jw  = JxW[q];
        PetscReal nuH,dNuH,beta,dbeta;
        ierr = PointwiseNuHAndBeta(feS,&wq,Dwq,&nuH,&dNuH,&beta,&dbeta);CHKERRQ(ierr);


        for (PetscInt k=0; k<4; k++) {   // Test functions
          for (PetscInt l=0; l<4; l++) { // Trial functions

            const FEFunctionGerm &test_qk=test[q][k];
            const FEFunctionGerm &test_ql=test[q][l];

            const PetscReal ht = test_qk.val,h = test_ql.val,
                  dxt = test_qk.dx, dyt = test_qk.dy,
                  dx = test_ql.dx, dy = test_ql.dy,

            // Cross terms appearing with beta'
            bvx = ht*wq.u,bvy = ht*wq.v,bux = wq.u*h,buy = wq.v*h,
            // Cross terms appearing with nuH'
            cvx = dxt*(2*Dwq[0]+Dwq[1]) + dyt*Dwq[2],
            cvy = dyt*(2*Dwq[1]+Dwq[0]) + dxt*Dwq[2],
            cux = (2*Dwq[0]+Dwq[1])*dx + Dwq[2]*dy,
            cuy = (2*Dwq[1]+Dwq[0])*dy + Dwq[2]*dx;
            
            // u-u coupling
            K[k*16+l*2]     += jw*(beta*ht*h + dbeta*bvx*bux + nuH*(2*dxt*dx + dyt*0.5*dy) + dNuH*cvx*cux);
            // u-v coupling
            K[k*16+l*2+1]   += jw*(dbeta*bvx*buy + nuH*(0.5*dyt*dx + dxt*dy) + dNuH*cvx*cuy);
            // v-u coupling
            K[k*16+8+l*2]   += jw*(dbeta*bvy*bux + nuH*(0.5*dxt*dy + dyt*dx) + dNuH*cvy*cux);
            // v-v coupling
            K[k*16+8+l*2+1] += jw*(beta*ht*h + dbeta*bvy*buy + nuH*(2*dyt*dy + dxt*0.5*dx) + dNuH*cvy*cuy);
          }
        }
      }
      ierr = dofmap.addLocalJacobianBlock(K,J);
      
      // ierr = MatSetValuesBlockedStencil(J,4,row,4,col,K,ADD_VALUES);CHKERRQ(ierr);
    }
  }
  if (bc_locations && vel_bc) {
    for (i=grid.xs; i<grid.xs+grid.xm; i++) {
      for (j=grid.ys; j<grid.ys+grid.ym; j++) {
        if (bc_locations->value(i,j) == MASK_SHEET) {
          const PetscReal ident[4] = {dirichletScale,0,0,dirichletScale};
          MatStencil row;
          // FIXME: Transpose shows up here!
          row.j = i; row.i = j;
          ierr = MatSetValuesBlockedStencil(J,1,&row,1,&row,ident,ADD_VALUES);CHKERRQ(ierr);
        }
      }
    }
  }

  if(bc_locations) {
    ierr = bc_locations->end_access();CHKERRQ(ierr);
  }
  if(vel_bc) {
    ierr = vel_bc->end_access(); CHKERRQ(ierr);
  }
  
  ierr = MatAssemblyBegin(J,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(J,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

  PetscTruth monitor_jacobian;
  ierr = PetscOptionsHasName(NULL,"-ssa_monitor_jacobian",&monitor_jacobian);CHKERRQ(ierr);
  if (monitor_jacobian) {
    ierr = PetscPrintf(grid.com,
                       "SSA Jacobian\n");
    CHKERRQ(ierr);
    ierr = MatView(J,PETSC_VIEWER_STDOUT_WORLD);
  }

  ierr = MatSetOption(J,MAT_NEW_NONZERO_LOCATION_ERR,PETSC_TRUE);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "SSAFEM::view"
PetscErrorCode SSAFEM::view(PetscViewer viewer)
{
  PetscErrorCode ierr;
  PetscTruth iascii;

  PetscFunctionBegin;
  ierr = PetscTypeCompare((PetscObject)viewer,PETSC_VIEWER_ASCII,&iascii);
           CHKERRQ(ierr);
  if (!iascii) {
    SETERRQ1(PETSC_ERR_SUP,"Viewer type %s not supported for SSA FE",
             ((PetscObject)viewer)->type_name);
  }
  ierr = SNESView(snes,viewer);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}


PetscErrorCode SSAFEFunction(DALocalInfo *info,
                             const PISMVector2 **xg, PISMVector2 **yg,
                             FECTX *fe)
{
  return fe->ssa->compute_local_function(info,xg,yg);
}


PetscErrorCode SSAFEJacobian(DALocalInfo *info,
                             const PISMVector2 **xg, Mat J,
                             FECTX *fe)
{
  return fe->ssa->compute_local_jacobian(info,xg,J);
}

