// Copyright (C) 2009--2012 Jed Brown and Ed Bueler and Constantine Khroulev and David Maxwell
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
#include "Mask.hh"
#include "basal_resistance.hh"

SSA *SSAFEMFactory(IceGrid &g, IceBasalResistancePlasticLaw &b,
                   EnthalpyConverter &ec, const NCConfigVariable &c)
{
  return new SSAFEM(g,b,ec,c);
}

//! \brief Allocating SSAFEM-specific objects; called by the constructor.
PetscErrorCode SSAFEM::allocate_fem() {
  PetscErrorCode ierr;

  dirichletScale = 1.0;
  ocean_rho = config.get("sea_water_density");
  earth_grav = config.get("standard_gravity");
  m_beta_ice_free_bedrock = config.get("beta_ice_free_bedrock");

  ierr = DMCreateGlobalVector(SSADA, &r);CHKERRQ(ierr);
  ierr = DMGetMatrix(SSADA, "baij", &J); CHKERRQ(ierr);

  ierr = SNESCreate(grid.com,&snes);CHKERRQ(ierr);

  // Default of maximum 200 iterations; possibly overridded by commandline
  PetscInt snes_max_it = 200;
  ierr = SNESSetTolerances(snes,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT,
                           snes_max_it,PETSC_DEFAULT);
  // ierr = SNESSetOptionsPrefix(snes,((PetscObject)this)->prefix);CHKERRQ(ierr);

  ierr = SNESSetFromOptions(snes);CHKERRQ(ierr);

  // Allocate feStore, which contains coefficient data at the quadrature points of all the elements.
  // There are nElement elements, and FEQuadrature::Nq quadrature points.
  PetscInt nElements = element_index.element_count();
  feStore = new FEStoreNode[FEQuadrature::Nq*nElements];

  // hardav IceModelVec2S is not used (so far).
  const PetscScalar power = 1.0 / flow_law->exponent();
  char unitstr[TEMPORARY_STRING_LENGTH];
  snprintf(unitstr, sizeof(unitstr), "Pa s%f", power);
  ierr = hardav.create(grid, "hardav", true); CHKERRQ(ierr);
  ierr = hardav.set_attrs("internal", "vertically-averaged ice hardness", unitstr, ""); CHKERRQ(ierr);

  return 0;
}

//! Undo the allocations of SSAFEM::allocate_fem; called by the destructor.
PetscErrorCode SSAFEM::deallocate_fem() {
  PetscErrorCode ierr;

  ierr = SNESDestroy(&snes);CHKERRQ(ierr);
  delete feStore;
  ierr = VecDestroy(&r); CHKERRQ(ierr);
  ierr = MatDestroy(&J); CHKERRQ(ierr);

  return 0;
}

// Initialize the solver, called once by the client before use.
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

//! Opportunity to modify behaviour based on command-line options.
/*! Called from SSAFEM::init */
PetscErrorCode SSAFEM::setFromOptions()
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = PetscOptionsHead("SSA FEM options");CHKERRQ(ierr);
  dirichletScale = 1.0e9;
  ierr = PetscOptionsReal("-ssa_fe_dirichlet_scale",
                          "Enforce Dirichlet conditions with this additional scaling",
                          "",
                          dirichletScale,
                          &dirichletScale,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsTail();CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

//! Solve the SSA.
PetscErrorCode SSAFEM::solve()
{
  PetscErrorCode ierr;
  PetscViewer    viewer;
  char           filename[PETSC_MAX_PATH_LEN];
  PetscBool     flg;

  m_epsilon_ssa = config.get("epsilon_ssa");

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
    ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);
  }

  // Set the SNES callbacks to call into our compute_local_function and compute_local_jacobian
  // methods via SSAFEFunction and SSAFEJ
  ierr = DMDASetLocalFunction(SSADA,(DMDALocalFunction1)SSAFEFunction);CHKERRQ(ierr);
  ierr = DMDASetLocalJacobian(SSADA,(DMDALocalFunction1)SSAFEJacobian);CHKERRQ(ierr);
  callback_data.da = SSADA;  callback_data.ssa = this;
  ierr = SNESSetDM(snes, SSADA); CHKERRQ(ierr);
  ierr = SNESSetFunction(snes, r,    SNESDAFormFunction,   &callback_data);CHKERRQ(ierr);
  ierr = SNESSetJacobian(snes, J, J, SNESDAComputeJacobian,&callback_data);CHKERRQ(ierr);

  stdout_ssa.clear();
  if (getVerbosityLevel() >= 2)
    stdout_ssa = "  SSA: ";

  // Set up the system to solve (store coefficient data at the quadrature points):
  ierr = setup(); CHKERRQ(ierr);

  // Solve:
  ierr = SNESSolve(snes,NULL,SSAX);CHKERRQ(ierr);

  // See if it worked.
  SNESConvergedReason reason;
  ierr = SNESGetConvergedReason( snes, &reason); CHKERRQ(ierr);
  if(reason < 0)
  {
    SETERRQ1(grid.com, 1,
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
    ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);
  }
  return 0;
}

//! Initialize stored data from the coefficients in the SSA.  Called by SSAFEM::solve.
/* This method is should be called after SSAFEM::init and whenever
any geometry or temperature related coefficients have changed. The method
stores the values of the coefficients at the quadrature points of each
element so that these interpolated values do not need to be computed
during each outer iteration of the nonlinear solve.*/
PetscErrorCode SSAFEM::setup()
{
  PetscReal      **h,
                 **H,
                 **topg,
                 **tauc_array,
                  *Enth_e[4],
                  *Enth_q[4],
                  **ds_x,
                  **ds_y;
  PetscInt         i,j,k,q,p,
                   Mz = grid.Mz;
  PetscErrorCode   ierr;

  PetscReal ice_rho = config.get("ice_density");

  for(q=0;q<FEQuadrature::Nq;q++)
  {
    Enth_q[q] = new PetscReal[Mz];
  }

  GeometryCalculator gc(sea_level, config);

  ierr = enthalpy->begin_access();CHKERRQ(ierr);
  bool driving_stress_explicit;
  if(surface != NULL) {
    driving_stress_explicit = false;
    ierr = surface->get_array(h);CHKERRQ(ierr);
  } else {
    driving_stress_explicit = true;
    ierr = driving_stress_x->get_array(ds_x);CHKERRQ(ierr);
    ierr = driving_stress_y->get_array(ds_y);CHKERRQ(ierr);
  }

  ierr = thickness->get_array(H);CHKERRQ(ierr);
  ierr = bed->get_array(topg);CHKERRQ(ierr);
  ierr = tauc->get_array(tauc_array);CHKERRQ(ierr);


  PetscInt xs = element_index.xs, xm = element_index.xm,
           ys = element_index.ys, ym = element_index.ym;
  for (i=xs; i<xs+xm; i++) {
    for (j=ys;j<ys+ym; j++) {
      PetscReal hq[FEQuadrature::Nq],hxq[FEQuadrature::Nq],hyq[FEQuadrature::Nq];
      PetscReal ds_xq[FEQuadrature::Nq], ds_yq[FEQuadrature::Nq];
      if(driving_stress_explicit) {
        quadrature.computeTrialFunctionValues(i,j,dofmap,ds_x,ds_xq);
        quadrature.computeTrialFunctionValues(i,j,dofmap,ds_y,ds_yq);
      } else {
        // Extract coefficient values at the quadrature points.
        quadrature.computeTrialFunctionValues(i,j,dofmap,h,hq,hxq,hyq);
      }

      PetscReal Hq[FEQuadrature::Nq], bq[FEQuadrature::Nq], taucq[FEQuadrature::Nq];
      quadrature.computeTrialFunctionValues(i,j,dofmap,H,Hq);
      quadrature.computeTrialFunctionValues(i,j,dofmap,topg,bq);
      quadrature.computeTrialFunctionValues(i,j,dofmap,tauc_array,taucq);

      const PetscInt ij = element_index.flatten(i,j);
      FEStoreNode *feS = &feStore[4*ij];
      for (q=0; q<4; q++) {
        feS[q].H  = Hq[q];
        feS[q].b  = bq[q];
        feS[q].tauc = taucq[q];
        if(driving_stress_explicit) {
          feS[q].driving_stress.u = ds_xq[q];
          feS[q].driving_stress.v = ds_yq[q];
        } else {
          feS[q].driving_stress.u = -ice_rho*earth_grav*Hq[q]*hxq[q];
          feS[q].driving_stress.v = -ice_rho*earth_grav*Hq[q]*hyq[q];
        }
        // feS[q].hx = hxq[q];
        // feS[q].hy = hyq[q];

        feS[q].mask = gc.mask(feS[q].b, feS[q].H);
      }

      // In the following, we obtain the averaged hardness value from enthalpy by
      // interpolating enthalpy in each column over a quadrature point and then
      // taking the average over the column.  A faster approach would be to take
      // the column average over each element nodes and then interpolate to the
      // quadrature points. Does this make a difference?

      // Obtain the values of enthalpy at each vertical level at each of the vertices
      // of the current element.
      ierr = enthalpy->getInternalColumn(i,j,&Enth_e[0]);
      ierr = enthalpy->getInternalColumn(i+1,j,&Enth_e[1]);CHKERRQ(ierr);
      ierr = enthalpy->getInternalColumn(i+1,j+1,&Enth_e[2]);CHKERRQ(ierr);
      ierr = enthalpy->getInternalColumn(i,j+1,&Enth_e[3]);CHKERRQ(ierr);

      // We now want to interpolate to the quadrature points at each of the
      // vertical levels.  It would be nice to use quadrature::computeTestFunctionValues,
      // but the way we have just obtained the values at the element vertices
      // using getInternalColumn doesn't make this straightforward.  So we compute the values
      // by hand.
      const FEFunctionGerm (*test)[FEQuadrature::Nk] = quadrature.testFunctionValues();
      for (k=0; k<Mz; k++) {
        Enth_q[0][k] = Enth_q[1][k] = Enth_q[2][k] = Enth_q[3][k] = 0;
        for (q=0; q<FEQuadrature::Nq; q++) {
          for (p=0; p<FEQuadrature::Nk; p++) {
            Enth_q[q][k] += test[q][p].val * Enth_e[p][k];
          }
        }
      }

      // Now, for each column over a quadrature point, find the averaged_hardness.
      for (q=0; q<FEQuadrature::Nq; q++) {
        // Evaluate column integrals in flow law at every quadrature point's column
        feS[q].B = flow_law->averaged_hardness(feS[q].H, grid.kBelowHeight(feS[q].H),
                                                  &grid.zlevels[0], Enth_q[q]);
      }
    }
  }
  if(surface != NULL) {
    ierr = surface->end_access();CHKERRQ(ierr);
  } else {
    ierr = driving_stress_x->end_access();CHKERRQ(ierr);
    ierr = driving_stress_y->end_access();CHKERRQ(ierr);
  }
  ierr = thickness->end_access();CHKERRQ(ierr);
  ierr = bed->end_access();CHKERRQ(ierr);
  ierr = tauc->end_access();CHKERRQ(ierr);
  ierr = enthalpy->end_access();CHKERRQ(ierr);

  for(q=0;q<4;q++)
  {
    delete [] Enth_q[q];
  }

  return 0;
}

//!\brief Compute the "2 x (effective viscosity) x height' and effective viscous bed strength from
//! the current solution, at a single quadrature point.
/*! The coefficient data at the quadrature point comes from \a feS.
The value of the solution and its symmetric gradient comes \a u and \a Du.
The function returns the values and the derivatives with respect
to the solution in the output variables \a nuH, \a dNuH, \a beta, and \a dbeta.
Use NULL pointers if no derivatives are desired.
*/
inline PetscErrorCode SSAFEM::PointwiseNuHAndBeta(const FEStoreNode *feS,
                                                  const PISMVector2 *u,const PetscReal Du[],
                                                  PetscReal *nuH, PetscReal *dNuH,
                                                  PetscReal *beta, PetscReal *dbeta)
{

  Mask M;

  if (feS->H < strength_extension->get_min_thickness()) {
    *nuH = strength_extension->get_notional_strength();
    if (dNuH) *dNuH = 0;
  } else {
    flow_law->effective_viscosity_with_derivative(feS->B, Du, nuH, dNuH);
    *nuH  *= feS->H;
    *nuH  += m_epsilon_ssa;
    if (dNuH) *dNuH *= feS->H;
  }
  *nuH  *=  2;
  if (dNuH) *dNuH *= 2;

  if( M.grounded_ice(feS->mask) )
  {
    basal.dragWithDerivative(feS->tauc,u->u,u->v,beta,dbeta);
  } else {  
    *beta = 0;
    if( M.ice_free_land(feS->mask) )
    {
      *beta = m_beta_ice_free_bedrock;
    }
    if(dbeta) *dbeta = 0;
  }
  return 0;
}

//! \brief Sets Dirichlet boundary conditions. Called from SSAFEFunction and
//! SSAFEJacobian.
/*! If for some vertex \a local_bc_mask indicates that it
is an explicit Dirichlet node, the values of x for that node is set from the Dirichlet
data BC_vel. The row and column in the \a dofmap are set as invalid.
This last step ensures that the residual and Jacobian entries
corresponding to a Dirichlet unknown are not set in the main loops of
SSAFEM::compute_local_function and SSSAFEM:compute_local_jacobian.
*/
void SSAFEM::FixDirichletValues( PetscReal local_bc_mask[],PISMVector2 **BC_vel,
                                PISMVector2 x[], FEDOFMap &my_dofmap)
{
  for (PetscInt k=0; k<4; k++) {
    if (PismIntMask(local_bc_mask[k]) == 1) { // Dirichlet node
      PetscInt ii, jj;
      my_dofmap.localToGlobal(k,&ii,&jj);
      x[k].u = BC_vel[ii][jj].u;
      x[k].v = BC_vel[ii][jj].v;
      // Mark any kind of Dirichlet node as not to be touched
      my_dofmap.markRowInvalid(k);
      my_dofmap.markColInvalid(k);
    }
  }
}

//! Implements the callback for computing the SNES local function.
/*! Compute the residual \f[r_{ij}= G(x,\psi_{ij}) \f] where \f$G\f$ is the weak form of the SSA, \f$x\f$
is the current approximate solution, and the \f$\psi_{ij}\f$ are test functions. */
PetscErrorCode SSAFEM::compute_local_function(DMDALocalInfo *info, const PISMVector2 **xg, PISMVector2 **yg)
{
  PetscInt         i,j,k,q;
  PetscReal        **bc_mask;
  PISMVector2        **BC_vel;
  PetscErrorCode   ierr;

  (void) info; // Avoid compiler warning.

  // Zero out the portion of the function we are responsible for computing.
  for (i=grid.xs; i<grid.xs+grid.xm; i++) {
    for (j=grid.ys; j<grid.ys+grid.ym; j++) {
      yg[i][j].u = yg[i][j].v = 0;
    }
  }

  // Start access of Dirichlet data, if present.
  if (bc_locations && vel_bc) {
    ierr = bc_locations->get_array(bc_mask);CHKERRQ(ierr);
    ierr = vel_bc->get_array(BC_vel); CHKERRQ(ierr);
  }

  // Jacobian times weights for quadrature.
  PetscScalar JxW[FEQuadrature::Nq];
  quadrature.getWeightedJacobian(JxW);

  // Storage for the current solution at quadrature points.
  PISMVector2 u[FEQuadrature::Nq];
  PetscScalar Du[FEQuadrature::Nq][3];

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
      // Storage for element-local solution and residuals.
      PISMVector2     x[4],y[4];
      // Index into coefficient storage in feStore
      const PetscInt ij = element_index.flatten(i,j);

      // Initialize the map from global to local degrees of freedom for this element.
      dofmap.reset(i,j,grid);

      // Obtain the value of the solution at the nodes adjacent to the element.
      dofmap.extractLocalDOFs(i,j,xg,x);

      // These values now need to be adjusted if some nodes in the element have
      // Dirichlet data.
      if(bc_locations && vel_bc) {
        dofmap.extractLocalDOFs(i,j,bc_mask,local_bc_mask);
        FixDirichletValues(local_bc_mask,BC_vel,x,dofmap);
      }

      // Zero out the element-local residual in prep for updating it.
      for(k=0;k<FEQuadrature::Nk;k++){
        y[k].u = 0; y[k].v = 0;
      }

      // Compute the solution values and symmetric gradient at the quadrature points.
      quadrature.computeTrialFunctionValues(x,u,Du);

      for (q=0; q<FEQuadrature::Nq; q++) {     // loop over quadrature points on this element.

        // Symmetric gradient at the quadrature point.
        PetscScalar *Duq = Du[q];

        // Coefficients and weights for this quadrature point.
        const FEStoreNode *feS = &feStore[ij*FEQuadrature::Nq+q];
        const PetscReal    jw  = JxW[q];
        PetscReal nuH, beta;
        ierr = PointwiseNuHAndBeta(feS,u+q,Duq,&nuH,NULL,&beta,NULL);CHKERRQ(ierr);

        // The next few lines compute the actual residual for the element.
        PISMVector2 f;
        f.u = beta*u[q].u - feS->driving_stress.u;
        f.v = beta*u[q].v - feS->driving_stress.v;

        for(k=0; k<4;k++) {  // loop over the test functions.
          const FEFunctionGerm &testqk = test[q][k];
          y[k].u += jw*(nuH*(testqk.dx*(2*Duq[0]+Duq[1]) + testqk.dy*Duq[2]) + testqk.val*f.u);
          y[k].v += jw*(nuH*(testqk.dy*(2*Duq[1]+Duq[0]) + testqk.dx*Duq[2]) + testqk.val*f.v);
        }
      } // q

      dofmap.addLocalResidualBlock(y,yg);
    } // j-loop
  } // i-loop

  // Until now we have not touched rows in the residual corresponding to Dirichlet data.
  // We fix this now.
  if (bc_locations && vel_bc) {
    // Enforce Dirichlet conditions strongly
    for (i=grid.xs; i<grid.xs+grid.xm; i++) {
      for (j=grid.ys; j<grid.ys+grid.ym; j++) {
        if (bc_locations->as_int(i,j) == 1) {
          // Enforce explicit dirichlet data.
          yg[i][j].u = dirichletScale * (xg[i][j].u - BC_vel[i][j].u);
          yg[i][j].v = dirichletScale * (xg[i][j].v - BC_vel[i][j].v);
        }
      }
    }
    ierr = bc_locations->end_access();CHKERRQ(ierr);
    ierr = vel_bc->end_access(); CHKERRQ(ierr);
  }

  PetscBool monitorFunction;
  ierr = PetscOptionsHasName(NULL,"-ssa_monitor_function",&monitorFunction);CHKERRQ(ierr);
  if (monitorFunction) {
    ierr = PetscPrintf(grid.com,"SSA Solution and Function values (pointwise residuals)\n");CHKERRQ(ierr);
    for (i=grid.xs; i<grid.xs+grid.xm; i++) {
      for (j=grid.ys; j<grid.ys+grid.ym; j++) {
        ierr = PetscSynchronizedPrintf(grid.com,
                                       "[%2d,%2d] u=(%12.10e,%12.10e)  f=(%12.4e,%12.4e)\n",
                                       i,j,xg[i][j].u,xg[i][j].v,yg[i][j].u,yg[i][j].v);CHKERRQ(ierr);
      }
    }
    ierr = PetscSynchronizedFlush(grid.com);CHKERRQ(ierr);
  }

  return 0;
}



//! Implements the callback for computing the SNES local Jacobian.
/*! Compute the Jacobian \f[J_{ij}{kl} \frac{d r_{ij}}{d x_{kl}}= G(x,\psi_{ij}) \f]
where \f$G\f$ is the weak form of the SSA, \f$x\f$ is the current approximate solution, and
the \f$\psi_{ij}\f$ are test functions. */
PetscErrorCode SSAFEM::compute_local_jacobian(DMDALocalInfo *info, const PISMVector2 **xg, Mat Jac )
{

  PetscReal      **bc_mask;
  PISMVector2    **BC_vel;
  PetscInt         i,j;
  PetscErrorCode   ierr;

  // Avoid compiler warning.
  (void) info;

  // Zero out the Jacobian in preparation for updating it.
  ierr = MatZeroEntries(Jac);CHKERRQ(ierr);

  // Start access to Dirichlet data if present.
  if (bc_locations && vel_bc) {
    ierr = bc_locations->get_array(bc_mask);CHKERRQ(ierr);
    ierr = vel_bc->get_array(BC_vel); CHKERRQ(ierr);
  }

  // Jacobian times weights for quadrature.
  PetscScalar JxW[FEQuadrature::Nq];
  quadrature.getWeightedJacobian(JxW);

  // Storage for the current solution at quadrature points.
  PISMVector2 w[FEQuadrature::Nq];
  PetscScalar Dw[FEQuadrature::Nq][3];

  // Values of the finite element test functions at the quadrature points.
  // This is an Nq by Nk array of function germs (Nq=#of quad pts, Nk=#of test functions).
  const FEFunctionGerm (*test)[FEQuadrature::Nk] = quadrature.testFunctionValues();

  // Flags for each vertex in an element that determine if explicit Dirichlet data has
  // been set.
  PetscReal local_bc_mask[FEQuadrature::Nk];

  // Loop through all the elements.
  PetscInt xs = element_index.xs, xm = element_index.xm,
           ys = element_index.ys, ym = element_index.ym;
  for (i=xs; i<xs+xm; i++) {
    for (j=ys; j<ys+ym; j++) {
      // Values of the solution at the nodes of the current element.
      PISMVector2    x[FEQuadrature::Nk];

      // Element-local Jacobian matrix (there are FEQuadrature::Nk vector valued degrees
      // of freedom per elment, for a total of (2*FEQuadrature::Nk)*(2*FEQuadrature::Nk) = 16
      // entries in the local Jacobian.
      PetscReal      K[(2*FEQuadrature::Nk)*(2*FEQuadrature::Nk)];

      // Index into the coefficient storage array.
      const PetscInt ij = element_index.flatten(i,j);

      // Initialize the map from global to local degrees of freedom for this element.
      dofmap.reset(i,j,grid);

      // Obtain the value of the solution at the adjacent nodes to the element.
      dofmap.extractLocalDOFs(i,j,xg,x);

      // These values now need to be adjusted if some nodes in the element have
      // Dirichlet data.
      if(bc_locations && vel_bc) {
        dofmap.extractLocalDOFs(i,j,bc_mask,local_bc_mask);
        FixDirichletValues(local_bc_mask,BC_vel,x,dofmap);
      }

      // Compute the values of the solution at the quadrature points.
      quadrature.computeTrialFunctionValues(x,w,Dw);

      // Build the element-local Jacobian.
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

            // FIXME (DAM 2/28/11) The following computations could be a little better documented.
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

            if(nuH==0)
            {
              verbPrintf(1,grid.com,"nuh=0 i %d j %d q %d k %d\n",i,j,q,k);
            }
            // u-u coupling
            K[k*16+l*2]     += jw*(beta*ht*h + dbeta*bvx*bux + nuH*(2*dxt*dx + dyt*0.5*dy) + dNuH*cvx*cux);
            // u-v coupling
            K[k*16+l*2+1]   += jw*(dbeta*bvx*buy + nuH*(0.5*dyt*dx + dxt*dy) + dNuH*cvx*cuy);
            // v-u coupling
            K[k*16+8+l*2]   += jw*(dbeta*bvy*bux + nuH*(0.5*dxt*dy + dyt*dx) + dNuH*cvy*cux);
            // v-v coupling
            K[k*16+8+l*2+1] += jw*(beta*ht*h + dbeta*bvy*buy + nuH*(2*dyt*dy + dxt*0.5*dx) + dNuH*cvy*cuy);
          } // l
        } // k
      } // q
      ierr = dofmap.addLocalJacobianBlock(K,Jac);
    } // j
  } // i


  // Until now, the rows and columns correspoinding to Dirichlet data have not been set.  We now
  // put an identity block in for these unknowns.  Note that because we have takes steps to not touching these
  // columns previously, the symmetry of the Jacobian matrix is preserved.
  if (bc_locations && vel_bc) {
    for (i=grid.xs; i<grid.xs+grid.xm; i++) {
      for (j=grid.ys; j<grid.ys+grid.ym; j++) {
        if (bc_locations->as_int(i,j) == 1) {
          const PetscReal ident[4] = {dirichletScale,0,0,dirichletScale};
          MatStencil row;
          // FIXME: Transpose shows up here!
          row.j = i; row.i = j;
          ierr = MatSetValuesBlockedStencil(Jac,1,&row,1,&row,ident,ADD_VALUES);CHKERRQ(ierr);
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

  ierr = MatAssemblyBegin(Jac,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(Jac,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

  PetscBool monitor_jacobian;
  ierr = PetscOptionsHasName(NULL,"-ssa_monitor_jacobian",&monitor_jacobian);CHKERRQ(ierr);
  if (monitor_jacobian) {
    PetscViewer    viewer;

    char           file_name[PETSC_MAX_PATH_LEN];
    PetscInt iter;
    ierr = SNESGetIterationNumber(snes,&iter);
    snprintf(file_name,  PETSC_MAX_PATH_LEN, "PISM_SSAFEM_J%d.m",iter);

      ierr = verbPrintf(2, grid.com,
                 "writing Matlab-readable file for SSAFEM system A xsoln = rhs to file `%s' ...\n",
                 file_name); CHKERRQ(ierr);
      ierr = PetscViewerCreate(grid.com, &viewer);CHKERRQ(ierr);
      ierr = PetscViewerSetType(viewer, PETSCVIEWERASCII);CHKERRQ(ierr);
      ierr = PetscViewerSetFormat(viewer, PETSC_VIEWER_ASCII_MATLAB);CHKERRQ(ierr);
      ierr = PetscViewerFileSetName(viewer, file_name);CHKERRQ(ierr);

      ierr = PetscObjectSetName((PetscObject) Jac,"A"); CHKERRQ(ierr);
      ierr = MatView(Jac, viewer);CHKERRQ(ierr);
  }

  ierr = MatSetOption(Jac,MAT_NEW_NONZERO_LOCATION_ERR,PETSC_TRUE);CHKERRQ(ierr);


  PetscFunctionReturn(0);
}

//!
PetscErrorCode SSAFEFunction(DMDALocalInfo *info,
                             const PISMVector2 **xg, PISMVector2 **yg,
                             SSAFEM_SNESCallbackData *fe)
{
  return fe->ssa->compute_local_function(info,xg,yg);
}


PetscErrorCode SSAFEJacobian(DMDALocalInfo *info,
                             const PISMVector2 **xg, Mat J,
                             SSAFEM_SNESCallbackData *fe)
{
  return fe->ssa->compute_local_jacobian(info,xg,J);
}

