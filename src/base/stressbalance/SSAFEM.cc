// Copyright (C) 2009--2011 Jed Brown and Ed Bueler and Constantine Khroulev
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
#include "SSAFEM_util.hh"


IceGridElementIndexer::IceGridElementIndexer(const IceGrid &g)
{
  // Start by assuming ghost elements exist in all directions.
  // If there are 'mx' by 'my' points in the grid and no 
  // ghosts, there are (mx-1) by (my-1) elements.  In the presence
  // of ghosts there are then (mx+1) by (my+1) elements.
  xs= g.xs-1; xm = g.xm+1;
  ys= g.ys-1; ym = g.ym+1;

  // Now correct if needed. The only way there will not be ghosts is if the 
  // grid is not periodic and we are up against the grid boundary.
  
  if( !(g.periodicity & X_PERIODIC) )
  {
    // First element has x-index 0.
    if(xs < 0){
      xs = 0;
    }
    // Total number of elements is Mx-1, so xs+xm should be no larger.
    if(xs+xm > g.Mx-1) {
      xm = g.Mx-1-xs;
    }
  }

  if( !(g.periodicity & Y_PERIODIC) )
  {
    if(ys < 0){
      ys = 0;
    }
    if(ys+ym > g.My-1) {
      ym = g.My-1-ys;
    }
  }
}

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
  {
    DALocalInfo info;
    PetscInt    nElements;
    ierr = DAGetLocalInfo(SSADA,&info);CHKERRQ(ierr);

    nElements = element_index.element_count();

    // We have a struct for the feStore at each quadrature point
    ierr = PetscMalloc(4*nElements*sizeof(feStore[0]),&feStore);CHKERRQ(ierr);

    // // sbs probably refers to "store block size". In the current code it is
    // // equal to 1, i.e. one value of vertically-averaged ice hardness per
    // // corner (or quadrature point?).
    sbs = 1;
    ierr = PetscMalloc(4*nElements*sbs*sizeof(PetscReal),&integratedStore);CHKERRQ(ierr);
  }

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
  ierr = PetscFree(integratedStore);CHKERRQ(ierr);
  ierr = PetscFree(feStore);CHKERRQ(ierr);
  ierr = VecDestroy(r); CHKERRQ(ierr);
  ierr = MatDestroy(J); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode SSAFEM::init(PISMVars &vars) {
  PetscErrorCode ierr;

  ierr = SSA::init(vars); CHKERRQ(ierr);
  ierr = verbPrintf(2,grid.com,
                    "  [using the finite element method implementation by Jed Brown]\n"); CHKERRQ(ierr);

  ierr = setFromOptions(); CHKERRQ(ierr);

  // On restart, SSA::init() reads the SSA velocity from a PISM output file
  // into IceModelVec2V "velocity". We use that field as an initial guess:

  ierr = velocity.copy_to(SSAX); CHKERRQ(ierr);

  // If we are not restarting from a PISM file, "velocity" is identically zero,
  // and the call above clears SSAX.

  return 0;
}

#undef __FUNCT__
#define __FUNCT__ "setFromOptions"
PetscErrorCode SSAFEM::setFromOptions()
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = PetscOptionsHead("SSA FE options");CHKERRQ(ierr);
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
  
  // Set up the system to solve:
  ierr = setup(); CHKERRQ(ierr);

  // Solve:
  ierr = SNESSolve(snes,NULL,SSAX);CHKERRQ(ierr);
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

  if (getVerbosityLevel() >= 2) 
    stdout_ssa = "  SSA: " + stdout_ssa;

  ierr = velocity.copy_from(SSAX); CHKERRQ(ierr);

  // communicate so that the ghost values are updated (for the geometry update,
  // etc)
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
#undef __FUNCT__
#define __FUNCT__ "setup"
PetscErrorCode SSAFEM::setup()
{
  DALocalInfo      info_struct,
                  *info = &info_struct; // DAGetLocalInfo writes into our memory
  PetscReal      **h,
                 **H,
                 **topg,
                 **tauc_array,
                  *Enth_e[4],
                  *Enth_q[4];
  PetscReal        jinvDiag[2];
  PetscInt         i,j,k,q,p,
                   Mz = grid.Mz;
  PetscErrorCode   ierr;

  PetscFunctionBegin;
  // We only differentiate dimensional quantities (not velocity) here
  jinvDiag[0] = 2/grid.dx;
  jinvDiag[1] = 2/grid.dy;

  ierr = PetscMalloc4(Mz,PetscReal,&Enth_q[0],
                      Mz,PetscReal,&Enth_q[1],
                      Mz,PetscReal,&Enth_q[2],
                      Mz,PetscReal,&Enth_q[3]);CHKERRQ(ierr);
  ierr = enthalpy->begin_access();CHKERRQ(ierr);
  ierr = surface->get_array(h);CHKERRQ(ierr);
  ierr = thickness->get_array(H);CHKERRQ(ierr);
  ierr = bed->get_array(topg);CHKERRQ(ierr);
  ierr = tauc->get_array(tauc_array);CHKERRQ(ierr);
  ierr = DAGetLocalInfo(SSADA,info);CHKERRQ(ierr);
  // See SSAFEFunction for discussion of communication

  PetscInt xs = element_index.xs, xm = element_index.xm,
           ys = element_index.ys, ym = element_index.ym;  
  for (i=xs; i<xs+xm; i++) {
    for (j=ys;j<ys+ym; j++) {
      const PetscInt ij = element_index.flatten(i,j);
      // feStore is interpreted as an array of FE elements using C-style
      // ordering, with 4 quadrature points per element.
      // 
      // Therefore 4*ij is the index of the first (0-th) point corresponding to
      // the element i,j, 4*ij + 1 is the index of the second point, etc
      PetscReal he[4],He[4],be[4],hq[4],Hq[4],bq[4],taue[4],tauq[4],hxq[4],hyq[4];
      FEStoreNode *feS;

      QuadExtractScalar(i,j,h,he); // surface height
      QuadExtractScalar(i,j,H,He); // thickness
      QuadExtractScalar(i,j,topg,be); // bed elevation
      QuadExtractScalar(i,j,tauc_array,taue); // basal friction
      QuadMatMultScalar(interp,he,hq);
      QuadMatMultScalar(interp,He,Hq);
      QuadMatMultScalar(interp,be,bq);
      QuadMatMultScalar(interp,taue,tauq);
      QuadMatMultScalar(derivx,he,hxq);
      QuadMatMultScalar(derivy,he,hyq);
      feS = &feStore[4*ij];
      for (q=0; q<4; q++) {
        feS[q].h  = hq[q];
        feS[q].H  = Hq[q];
        feS[q].b  = bq[q];
        feS[q].tauc = tauq[q];
        feS[q].hx = jinvDiag[0] * hxq[q];
        feS[q].hy = jinvDiag[1] * hyq[q];
      }

      // Surface and thickness information is stored, now do the thermal stuff
      ierr = enthalpy->getInternalColumn(i,j,&Enth_e[0]);CHKERRQ(ierr);
      ierr = enthalpy->getInternalColumn(i+1,j,&Enth_e[1]);CHKERRQ(ierr);
      ierr = enthalpy->getInternalColumn(i+1,j+1,&Enth_e[2]);CHKERRQ(ierr);
      ierr = enthalpy->getInternalColumn(i,j+1,&Enth_e[3]);CHKERRQ(ierr);
      // Interpolate to quadrature points at every vertical level
      for (k=0; k<Mz; k++) { // This loop could be cut short at the surface.
        Enth_q[0][k] = Enth_q[1][k] = Enth_q[2][k] = Enth_q[3][k] = 0;
        for (q=0; q<4; q++) {
          for (p=0; p<4; p++) {
            Enth_q[q][k] += interp[q*4+p] * Enth_e[p][k];
          }
        }
      }
      for (q=0; q<4; q++) {
        // Evaluate column integrals in flow law at every quadrature point's column
        PetscReal *iS = &integratedStore[sbs*(ij*4+q)]; // Location to put the stored data
        *iS = ice.averagedHardness_from_enth(feS[q].H, grid.kBelowHeight(feS[q].H),
          grid.zlevels, Enth_q[q]);
      }
    }
  }
  ierr = surface->end_access();CHKERRQ(ierr);
  ierr = thickness->end_access();CHKERRQ(ierr);
  ierr = bed->end_access();CHKERRQ(ierr);
  ierr = tauc->end_access();CHKERRQ(ierr);
  ierr = enthalpy->end_access();CHKERRQ(ierr);
  ierr = PetscFree4(Enth_q[0],Enth_q[1],Enth_q[2],Enth_q[3]);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

// The stores are dimensional, everything after \a u (inclusive) is nondimensional
#undef __FUNCT__
#define __FUNCT__ "PointwiseNuHAndBeta"
inline PetscErrorCode SSAFEM::PointwiseNuHAndBeta(const FEStoreNode *feS,const PetscReal *iS,
                                                  const PISMVector2 *u,const PetscReal Du[],
                                                  PetscReal *nuH, PetscReal *dNuH,
                                                  PetscReal *beta, PetscReal *dbeta)
{
  if (feS->H < strength_extension->get_min_thickness()) {
    *nuH = strength_extension->get_notional_strength();
    if (dNuH) *dNuH = 0;
  } else {
    PetscReal dimDu[3];
    for (int i=0; i<3; i++) dimDu[i] = ref.StrainRate() * Du[i];
    ice.effectiveViscosity_with_derivative(*iS, dimDu, nuH, dNuH);
    *nuH  *= feS->H;
    if (dNuH) *dNuH *= feS->H;
  }
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
    basal.dragWithDerivative(feS->tauc,
                             u->u*ref.Velocity(),u->v*ref.Velocity(),
                             beta,dbeta);
  }
  // Return nondimensional values, the factor of 2 comes from a different
  //   definition of integrated viscosity
  *nuH  *=  2 / ref.IntegratedViscosity();
  if (dNuH) *dNuH *= 2 * ref.StrainRate2() / ref.IntegratedViscosity(); // The derivative with respect to second invariant
  *beta /= ref.Drag();
  if (dbeta) *dbeta *= ref.Velocity2() / ref.Drag();
  PismValidFriction(*beta);
  PetscFunctionReturn(0);
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
                               MatStencil row[],MatStencil col[],PISMVector2 x[])
{
  for (PetscInt k=0; k<4; k++) {
    if (PismIntMask(lmask[k]) == MASK_SHEET) {
      x[k].u = BC_vel[col[k].j][col[k].i].u / ref.Velocity();
      x[k].v = BC_vel[col[k].j][col[k].i].v / ref.Velocity();
      // FIXME (DAM 2/11/11): Find the right negative number to use to indicate an invalid row/column.  
      row[k].j = row[k].i = PETSC_MIN_INT/10;
      col[k].j = col[k].i = PETSC_MIN_INT/10;
    }
  }
}


PetscErrorCode SSAFEM::compute_local_function(DALocalInfo *info, const PISMVector2 **xg, PISMVector2 **yg)
{
  PetscInt         i,j,k,q;
  PetscReal        jacDiag[2],jinvDiag[2],jdet;
  PetscReal      **mask, **H, **topg;
  PISMVector2        **BC_vel;
  PetscTruth       flg;
  PetscErrorCode   ierr;

  PetscInt rank;
  ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank);CHKERRQ(ierr);

  PetscFunctionBegin;
  if (!finite(xg[info->ys][info->xs].u)) SETERRQ(1,__FUNCT__ " called with non-finite value");
  // ierr = verbPrintf(5,grid.com,"In %s\n",__FUNCT__);CHKERRQ(ierr);

  for (i=grid.xs; i<grid.xs+grid.xm; i++) {
    for (j=grid.ys; j<grid.ys+grid.ym; j++) {
      yg[i][j].u = yg[i][j].v = 0;
    }
  }
  
  // Since we use uniform cartesian coordinates, the Jacobian is constant and diagonal on every element.
  // Note that the reference element is \f$ [-1,1]^2 \f$ hence the extra factor of 2.
  jacDiag[0] = 0.5*grid.dx/ref.Length();
  jacDiag[1] = 0.5*grid.dy/ref.Length();
  //if (fe->debug.rescale) jacDiag[0] = jacDiag[1] = 1;
  jinvDiag[0] = 1/jacDiag[0];
  jinvDiag[1] = 1/jacDiag[1];
  jdet = jacDiag[0]*jacDiag[1];

  ierr = PetscOptionsHasName(NULL,"-ssa_monitor_solution",&flg);CHKERRQ(ierr);
  if (flg) {
    ierr = PetscPrintf(grid.com,
                       "SSA Pointwise solution values (evaluating function at this value)\n");
    CHKERRQ(ierr);
    for (i=grid.xs; i<grid.xs+grid.xm; i++) {
      for (j=grid.ys; j<grid.ys+grid.ym; j++) {
        ierr = PetscPrintf(grid.com,"[%2d,%2d] (%g,%g)\n",i,j,
                           xg[i][j].u,xg[i][j].v);CHKERRQ(ierr);
      }
    }
  }
  // Consider the node layout in 1 dimension for ranks A,B,C
  // Owning rank:     A   |   B   |   C
  // Node index :   0 1 2 | 3 4 5 | 6 7 8
  //
  // Rank B has up-to-date values of the current iterate for [2,3,4,5,6] but we
  // only set function values for the owned portion [3,4,5]. In the finite
  // element context, rank B will evaluate weak forms on elements
  // [(2,3),(3,4),(4,5),(5,6)] but will only insert values for the owned nodes.
  // This is consistent with the communication pattern natural to finite
  // difference methods for which the DA object is designed.
  //
  // In the finite element context, the usual communication pattern is to have
  // the numbers (0..8) above refer to elements and decide on ownership of
  // interface degrees of freedom. Then each rank evaluates weak forms over the
  // locally owned elements and performs communication to sum into the global
  // vectors and matrices. With FD communication, a rank only inserts into
  // locally owned portions of the residual vector and rows of the global
  // matrix while with FE communication each rank adds values into locally
  // owned plus remotely owned interface parts of the residual vector and rows
  // of the global matrix.
  //
  // Note that PETSc uses the opposite meaning of x,y directions in the DA
  ierr = thickness->get_array(H);CHKERRQ(ierr);
  ierr = bed->get_array(topg);CHKERRQ(ierr);

  if (bc_locations && vel_bc) {
    ierr = bc_locations->get_array(mask);CHKERRQ(ierr);
    ierr = vel_bc->get_array(BC_vel); CHKERRQ(ierr);
  }

  PetscInt xs = element_index.xs, xm = element_index.xm,
           ys = element_index.ys, ym = element_index.ym;
  for (i=xs; i<xs+xm; i++) {
    for (j=ys; j<ys+ym; j++) {
      const PetscInt ij = element_index.flatten(i,j);
      PISMVector2     x[4],y[4];
      MatStencil row[4],col[4];

      QuadExtractVel(i,j,xg,x);
      ierr = QuadGetStencils(info,i,j,row,col);CHKERRQ(ierr);
      QuadZeroVel(y);

      if (bc_locations && vel_bc) {
        PetscReal lmask[4];
        QuadExtractScalar(i,j,mask,lmask);
        FixDirichletValues(lmask,BC_vel,row,col,x);
      }
      // \a x now contains correct velocities at Dirichlet nodes and \a row and
      // \a col indices for those nodes are set to -1


      for (q=0; q<numQuadPoints; q++) {     // loop over quadrature points
        const FEStoreNode *feS = &feStore[ij*4+q];
        const PetscReal *iS  = &integratedStore[(ij*4+q)*sbs];
        const PetscReal    jw  = jdet * quadWeights[q];
        PISMVector2 u,v;
        PetscReal nuH,beta,Du[3],Dv[3];
        ierr = QuadEvaluateVel(x,q,jinvDiag,&u,Du);CHKERRQ(ierr);
        ierr = PointwiseNuHAndBeta(feS,iS,&u,Du,&nuH,NULL,&beta,NULL);CHKERRQ(ierr);
        
          // Coefficients of test function (v) in weak form.
        v.u = beta*u.u + ice.rho * earth_grav * feS->H * feS->hx
          / ref.DrivingStress();
        v.v = beta*u.v + ice.rho * earth_grav * feS->H * feS->hy
          / ref.DrivingStress();
        // printf("draging %g driving %g\n",beta*u.v,ice.rho * earth_grav * feS->H * feS->hy
        //   / ref.DrivingStress());

        Dv[0] = nuH * Du[0];
        Dv[1] = nuH * Du[1];
        Dv[2] = nuH * Du[2];

        // Sum residuals over test functions
        for (k=0; k<4; k++) {   // Loop over test functions
          const PetscInt qk = q*4+k; // this indexing transposes matrices
          y[k].u += interp[qk] * jw * v.u                     // basal drag and driving force
            + derivx[qk] * jinvDiag[0] * jw * (2*Dv[0]+Dv[1]) // diagonal term in colon plus the trace term
            + derivy[qk] * jinvDiag[1] * jw * Dv[2];          // off-diagonal term "half is summed twice"
          y[k].v += interp[qk] * jw * v.v                     // basal drag and driving force
            + derivy[qk] * jinvDiag[1] * jw * (2*Dv[1]+Dv[0]) // diagonal term in colon plus the trace term
            + derivx[qk] * jinvDiag[0] * jw * Dv[2];          // off-diagonal term "half is summed twice"
        }
      }
      QuadInsertVel(row,y,yg);

    } // j-loop
  } // i-loop
  ierr = thickness->end_access();CHKERRQ(ierr);
  ierr = bed->end_access();CHKERRQ(ierr);

  if (bc_locations && vel_bc) {
    // Enforce Dirichlet conditions strongly
    for (i=grid.xs; i<grid.xs+grid.xm; i++) {
      for (j=grid.ys; j<grid.ys+grid.ym; j++) {
        if (bc_locations->value(i,j) == MASK_SHEET) {
          // Enforce zero sliding strongly
          yg[i][j].u = dirichletScale * (xg[i][j].u - BC_vel[i][j].u)
            / ref.Velocity();
          yg[i][j].v = dirichletScale * (xg[i][j].v - BC_vel[i][j].v)
            / ref.Velocity();
        }
      }
    }

    ierr = bc_locations->end_access();CHKERRQ(ierr);
    ierr = vel_bc->end_access(); CHKERRQ(ierr);
  }

  ierr = PetscOptionsHasName(NULL,"-ssa_monitor_function",&flg);CHKERRQ(ierr);
  if (flg) {
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
  PetscFunctionReturn(0);
}

PetscErrorCode SSAFEM::compute_local_jacobian(DALocalInfo *info, const PISMVector2 **xg, Mat J )
{
  PetscReal        jacDiag[2],jinvDiag[2],jdet;
  PetscReal      **mask;
  PISMVector2    **BC_vel;
  PetscInt         i,j;
  PetscErrorCode   ierr;
  PetscTruth     flg;

  PetscInt rank;
  ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank);CHKERRQ(ierr);

  // Since we use uniform cartesian coordinates, the Jacobian is constant and diagonal on every element.
  // Note that the reference element is \f$ [-1,1]^2 \f$ hence the extra factor of 2.

  jacDiag[0] = 0.5*grid.dx / ref.Length();
  jacDiag[1] = 0.5*grid.dy / ref.Length();
  //if (fe->debug.rescale) jacDiag[0] = jacDiag[1] = 1;
  jinvDiag[0] = 1/jacDiag[0];
  jinvDiag[1] = 1/jacDiag[1];
  jdet = jacDiag[0]*jacDiag[1];

  ierr = MatZeroEntries(J);CHKERRQ(ierr);

  if (bc_locations && vel_bc) {
    ierr = bc_locations->get_array(mask);CHKERRQ(ierr);
    ierr = vel_bc->get_array(BC_vel); CHKERRQ(ierr); 
  }
  
  

  PetscInt xs = element_index.xs, xm = element_index.xm,
           ys = element_index.ys, ym = element_index.ym;

         // for(i=xs;i<xs+xm;i++)
         // {
         //   for(j=ys;j<ys+ym;j++)
         //   {
         //     PISMVector2 pv = xg[i][j];
         //     printf("xg[%d,%d] %g %g\n",i,j,pv.u,pv.v);
         //   }
         // }


  for (i=xs; i<xs+xm; i++) { // Include ghost cells on either end
    for (j=ys; j<ys+ym; j++) {
      const PetscInt ij = element_index.flatten(i,j);
      MatStencil     row[4],col[4];
      PISMVector2        x[4],xh[4];
      PetscReal      K[8*8],lmask[4];
      const PetscScalar dh = 1.e-10;
      

      QuadExtractVel(i,j,xg,x);

      ierr = QuadGetStencils(info,i,j,row,col);CHKERRQ(ierr);
      if (bc_locations && vel_bc) {
        QuadExtractScalar(i,j,mask,lmask);
        FixDirichletValues(lmask,BC_vel,row,col,x);
      }
      ierr = PetscMemzero(K,sizeof(K));CHKERRQ(ierr);

      for (PetscInt q=0; q<numQuadPoints; q++) {
        const FEStoreNode *feS = &feStore[ij*4+q];
        const PetscReal   *iS  = &integratedStore[(ij*4+q)*sbs];
        const PetscReal    jw  = jdet*quadWeights[q];
        PISMVector2 w;
        PetscReal nuH,dNuH,beta,dbeta,Dw[3];
        ierr = QuadEvaluateVel(x,q,jinvDiag,&w,Dw);CHKERRQ(ierr);
        ierr = PointwiseNuHAndBeta(feS,iS,&w,Dw,&nuH,&dNuH,&beta,&dbeta);CHKERRQ(ierr);

        // J[k][l] = \partial r_k / \partial c_l

        for (PetscInt k=0; k<4; k++) {   // Test functions
          for (PetscInt l=0; l<4; l++) { // Trial functions

            const PetscInt qk = q*4+k,ql = q*4+l; // transpose and regular

            for(PetscInt ell=0;ell<4;ell++)
            {
              xh[ell]=x[ell];
            }
            PISMVector2 wh;
            PetscReal nuHh,dNuHh,betah,dbetah,Dwh[3];

            PetscReal Fu =   nuH*jw*(derivx[qk] * jinvDiag[0]  * (2*Dw[0]+Dw[1]) + derivy[qk] * jinvDiag[1] *  Dw[2] );
            PetscReal Fv =   nuH*jw*(derivy[qk] * jinvDiag[1] * (2*Dw[1]+Dw[0]) + derivx[qk] * jinvDiag[0] * Dw[2]);

            // Derivatives w.r.t u dof
            xh[l].u = x[l].u+dh; xh[l].v = x[l].v;
            ierr = QuadEvaluateVel(xh,q,jinvDiag,&wh,Dwh);CHKERRQ(ierr);
            ierr = PointwiseNuHAndBeta(feS,iS,&wh,Dwh,&nuHh,&dNuHh,&betah,&dbetah);CHKERRQ(ierr);
            PetscReal Fuhu = nuHh*jw*(derivx[qk] * jinvDiag[0]  * (2*Dwh[0]+Dwh[1]) + derivy[qk] * jinvDiag[1] *  Dwh[2] );
            PetscReal Fvhu = nuHh*jw*(derivy[qk] * jinvDiag[1] * (2*Dwh[1]+Dwh[0]) + derivx[qk] * jinvDiag[0] * Dwh[2]);

            // Derivatives w.r.t v dof
            xh[l].v = x[l].v+dh; xh[l].u = x[l].u;
            ierr = QuadEvaluateVel(xh,q,jinvDiag,&wh,Dwh);CHKERRQ(ierr);
            ierr = PointwiseNuHAndBeta(feS,iS,&wh,Dwh,&nuHh,&dNuHh,&betah,&dbetah);CHKERRQ(ierr);
            PetscReal Fuhv = nuHh*jw*(derivx[qk] * jinvDiag[0]  * (2*Dwh[0]+Dwh[1]) + derivy[qk] * jinvDiag[1] *  Dwh[2] );
            PetscReal Fvhv = nuHh*jw*(derivy[qk] * jinvDiag[1] * (2*Dwh[1]+Dwh[0]) + derivx[qk] * jinvDiag[0] * Dwh[2]);

            // Compute all finite differences.
            PetscReal dJFD[4];
            dJFD[0] = (Fuhu-Fu)/dh;
            dJFD[1] = (Fuhv-Fu)/dh;
            dJFD[2] = (Fvhu-Fv)/dh;
            dJFD[3] = (Fvhv-Fv)/dh;


            const PetscReal ht = interp[qk],h = interp[ql],
            dxt = derivx[qk]*jinvDiag[0],dyt = derivy[qk]*jinvDiag[1],
            dx = jinvDiag[0]*derivx[ql], dy = jinvDiag[1]*derivy[ql],
            // Cross terms appearing with beta'
            bvx = ht*w.u,bvy = ht*w.v,bux = w.u*h,buy = w.v*h,
            // Cross terms appearing with nuH'
            cvx = dxt*(2*Dw[0]+Dw[1]) + dyt*Dw[2],
            cvy = dyt*(2*Dw[1]+Dw[0]) + dxt*Dw[2],
            cux = (2*Dw[0]+Dw[1])*dx + Dw[2]*dy,
            cuy = (2*Dw[1]+Dw[0])*dy + Dw[2]*dx;
            
            // u-u coupling
            K[k*16+l*2]     += jw*(beta*ht*h + dbeta*bvx*bux + nuH*(2*dxt*dx + dyt*0.5*dy) + dNuH*cvx*cux);
            // u-v coupling
            K[k*16+l*2+1]   += jw*(dbeta*bvx*buy + nuH*(0.5*dyt*dx + dxt*dy) + dNuH*cvx*cuy);
            // v-u coupling
            K[k*16+8+l*2]   += jw*(dbeta*bvy*bux + nuH*(0.5*dxt*dy + dyt*dx) + dNuH*cvy*cux);
            // v-v coupling
            K[k*16+8+l*2+1] += jw*(beta*ht*h + dbeta*bvy*buy + nuH*(2*dyt*dy + dxt*0.5*dx) + dNuH*cvy*cuy);


            // // u-u coupling
            // PetscReal dJ[4];
            // dJ[0] = jw*(nuH*(2*dxt*dx + dyt*0.5*dy)+ dNuH*cvx*cux);
            // dJ[1] = jw*(nuH*(0.5*dyt*dx + dxt*dy) + dNuH*cvx*cuy);
            // dJ[2] = jw*(nuH*(0.5*dxt*dy + dyt*dx) + dNuH*cvy*cux);
            // dJ[3] = jw*(nuH*(2*dyt*dy + dxt*0.5*dx) + dNuH*cvy*cuy);
            // K[k*16+l*2]     +=  dJ[0];
            //   // u-v coupling
            // K[k*16+l*2+1]   +=  dJ[1];
            //   // v-u coupling
            // K[k*16+8+l*2]   += dJ[2];
            //   // v-v coupling
            // K[k*16+8+l*2+1] += dJ[3];
            // 
            // 
            // PetscScalar Jx = jw*nuH*(derivx[qk] * jinvDiag[0] * (2*Dw[0]+Dw[1]) + derivy[qk] * jinvDiag[1] * Dw[2]);
            // PetscScalar Jxh = jw*nuHh*(derivx[qk] * jinvDiag[0] * (2*Dwh[0]+Dwh[1]) + derivy[qk] * jinvDiag[1] * Dwh[2]);
            // PetscScalar dJcode = jw*(nuH*(2*dxt*dx + dyt*0.5*dy)+ dNuH*cvx*cux);
            // 
            // // printf("x[0] %g %g x[1] %g %g x[2] %g %g x[3] %g %g\n",x[0].u,x[0].v,x[1].u,x[1].v,x[2].u,x[2].v,x[3].u,x[3].v );
            // if(i<=1 && j<=1)
            // {
            //   // printf("Element (%d,%d) test %d trial %d\nuu: %g fd %g\nuv: %g fd %g\nvu: %g fd %g\nvv: %g fd %g\n",i,j,k,l,
            //   //   dJ[0],dJFD[0],dJ[1],dJFD[1],dJ[2],dJFD[2],dJ[3],dJFD[3]);
            // }
            // // printf("xh %g %g %g %g Dwh: %g Dj: %g code: %g\n",xh[0].u,xh[1].u,xh[2].u,xh[3].u,Dw[2],(Jxh-Jx)/dh,dJcode);
          }
        }
      }
      ierr = MatSetValuesBlockedStencil(J,4,row,4,col,K,ADD_VALUES);CHKERRQ(ierr);
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

  ierr = PetscOptionsHasName(NULL,"-ssa_monitor_jacobian",&flg);CHKERRQ(ierr);
  if (flg) {
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
  if (iascii) {
    ierr = PetscViewerASCIIPrintf(viewer,"FE: Store block size %d\n",
                                  sbs);CHKERRQ(ierr);
  } else {
    SETERRQ1(PETSC_ERR_SUP,"Viewer type %s not supported for SSA FE",
             ((PetscObject)viewer)->type_name);
  }
  ierr = SNESView(snes,viewer);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

//! \brief Computes vertically-averaged ice hardness on the standard grid.
PetscErrorCode SSAFEM::compute_hardav(IceModelVec2S &result) {
  PetscErrorCode ierr;
  const PetscScalar fillval = -0.01;
  PetscScalar *Eij; // columns of enthalpy values

  ierr = enthalpy->begin_access(); CHKERRQ(ierr);
  ierr = thickness->begin_access(); CHKERRQ(ierr);
  ierr = result.begin_access(); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      ierr = enthalpy->getInternalColumn(i,j,&Eij); CHKERRQ(ierr);
      const PetscScalar H = (*thickness)(i,j);
      if (H > 0.0) {
        result(i,j) = ice.averagedHardness_from_enth(H, grid.kBelowHeight(H),
                                                     grid.zlevels, Eij);
      } else { // put negative value below valid range
        result(i,j) = fillval;
      }
    }
  }
  ierr = enthalpy->end_access(); CHKERRQ(ierr);
  ierr = thickness->end_access(); CHKERRQ(ierr);
  ierr = result.end_access(); CHKERRQ(ierr);

  return 0;
}


#undef __FUNCT__
#define __FUNCT__ "SSAFEFunction"
PetscErrorCode SSAFEFunction(DALocalInfo *info,
                             const PISMVector2 **xg, PISMVector2 **yg,
                             FECTX *fe)
{
  return fe->ssa->compute_local_function(info,xg,yg);
}

#undef __FUNCT__
#define __FUNCT__ "SSAFEJacobian"
PetscErrorCode SSAFEJacobian(DALocalInfo *info,
                             const PISMVector2 **xg, Mat J,
                             FECTX *fe)
{
  return fe->ssa->compute_local_jacobian(info,xg,J);
}

