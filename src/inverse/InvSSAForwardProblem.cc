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

#include "InvSSAForwardProblem.hh"
#include "Mask.hh"
#include "basal_resistance.hh"
#include "PISMVars.hh"
#include <unistd.h> // for 'access'
#include <errno.h>
#include <sstream>
#include "PIO.hh"

/*! \brief Allocate PETSC structures needed for solving the various linearized
problems associated with the forward problem. */
PetscErrorCode InvSSAForwardProblem::allocate_ksp()
{
  PetscErrorCode ierr;

  // Storage for vector unknowns.
  ierr = DMGetMatrix(SSADA, "baij", &m_MatA); CHKERRQ(ierr);
  ierr = DMCreateLocalVector(SSADA, &m_VecU); CHKERRQ(ierr);
  ierr = DMCreateGlobalVector(SSADA, &m_VecZ2); CHKERRQ(ierr);
  ierr = DMCreateLocalVector(SSADA, &m_VecZ); CHKERRQ(ierr);
  ierr = DMCreateGlobalVector(SSADA, &m_VecRHS2); CHKERRQ(ierr);

  // Storage for scalar unknowns.
  ierr = DMGetMatrix(grid.da2, "baij", &m_MatB); CHKERRQ(ierr);
  ierr = DMCreateGlobalVector(grid.da2, &m_VecRHS); CHKERRQ(ierr);
  ierr = DMCreateGlobalVector(grid.da2, &m_VecV); CHKERRQ(ierr);

  ierr = KSPCreate(grid.com, &m_KSP); CHKERRQ(ierr);
  PetscReal ksp_rtol = 1e-12;
  ierr = KSPSetTolerances(m_KSP,ksp_rtol,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);
  PC pc;
  ierr = KSPGetPC(m_KSP,&pc); CHKERRQ(ierr);
  ierr = PCSetType(pc,PCBJACOBI); CHKERRQ(ierr);
  ierr = KSPSetFromOptions(m_KSP); CHKERRQ(ierr);


  ierr = KSPCreate(grid.com, &m_KSP_B); CHKERRQ(ierr);
  ierr = KSPSetTolerances(m_KSP_B,ksp_rtol,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);
  PC pc_B;
  ierr = KSPGetPC(m_KSP_B,&pc_B); CHKERRQ(ierr);
  ierr = PCSetType(pc_B,PCBJACOBI); CHKERRQ(ierr);
  ierr = KSPSetFromOptions(m_KSP_B); CHKERRQ(ierr);

  return 0;
}

/*! \brief Undoes allocation in allocate_ksp */
PetscErrorCode InvSSAForwardProblem::deallocate_ksp()
{
  PetscErrorCode ierr;

  if (m_KSP != PETSC_NULL) {
    ierr = KSPDestroy(&m_KSP); CHKERRQ(ierr);
  }

  if (m_KSP_B != PETSC_NULL) {
    ierr = KSPDestroy(&m_KSP_B); CHKERRQ(ierr);
  }

  if (m_MatA != PETSC_NULL) {
    ierr = MatDestroy(&m_MatA); CHKERRQ(ierr);
  }

  if (m_MatB != PETSC_NULL) {
    ierr = MatDestroy(&m_MatB); CHKERRQ(ierr);
  }

  if (m_VecU != PETSC_NULL) {
    ierr = VecDestroy(&m_VecU); CHKERRQ(ierr);
  }

  if (m_VecRHS2 != PETSC_NULL) {
    ierr = VecDestroy(&m_VecRHS2); CHKERRQ(ierr);
  }

  if (m_VecRHS != PETSC_NULL) {
    ierr = VecDestroy(&m_VecRHS); CHKERRQ(ierr);
  }

  if (m_VecZ != PETSC_NULL) {
    ierr = VecDestroy(&m_VecZ); CHKERRQ(ierr);
  }

  if (m_VecZ2 != PETSC_NULL) {
    ierr = VecDestroy(&m_VecZ2); CHKERRQ(ierr);
  }

  if (m_VecV != PETSC_NULL) {
    ierr = VecDestroy(&m_VecV); CHKERRQ(ierr);
  }

  return 0;
}

PetscErrorCode InvSSAForwardProblem::allocate_store()
{
  // We own an FEElementMap that knows how many finite element
  // elements we need to access. We use it to determine the
  // size our per-quadrature point storage.
  PetscInt nElements = element_index.element_count();
  m_dtauc_dzeta_store = new PetscReal[FEQuadrature::Nq*nElements];
  return 0;
}

PetscErrorCode InvSSAForwardProblem::deallocate_store(){
  delete [] m_dtauc_dzeta_store;
  return 0;
}


// Initialize the solver, called once by the client before use.
PetscErrorCode InvSSAForwardProblem::init(PISMVars &vars) {
  PetscErrorCode ierr;

  ierr = SSAFEM::init(vars); CHKERRQ(ierr);

  m_misfit_weight = dynamic_cast<IceModelVec2S*>(vars.get("vel_misfit_weight"));
  if (m_misfit_weight == NULL){
    verbPrintf(3,grid.com,"Weight for inverse problem L2 norm not available; using standard L2 norm.\n");
  }

  m_misfit_element_mask = dynamic_cast<IceModelVec2Int*>(vars.get("misfit_element_mask"));
  if (m_misfit_element_mask == NULL){
    verbPrintf(3,grid.com,"Misfit element mask not available; using all elements.\n");
  }

  return 0;
}

PetscErrorCode InvSSAForwardProblem::set_initial_velocity_guess(  IceModelVec2V &v )
{
  PetscErrorCode ierr;
  ierr = v.copy_to(SSAX); CHKERRQ(ierr);
  m_reassemble_T_matrix_needed = true;
  return 0;
}

// FIXME
/*! \brief apparently unused method! */
PetscErrorCode InvSSAForwardProblem::setup_vars()
{
  PetscErrorCode ierr;
  ierr = setup(); CHKERRQ(ierr);
  ierr = assemble_DomainNorm_matrix(); CHKERRQ(ierr);
  ierr = compute_range_l2_area(&m_range_l2_area);
  return 0;
}

PetscErrorCode InvSSAForwardProblem::set_zeta(IceModelVec2S &new_zeta )
{
  PetscErrorCode ierr;
  PetscInt i,j,q;


  PetscReal **zeta_array;
  ierr = new_zeta.get_array(zeta_array);CHKERRQ(ierr);
  PetscInt xs = element_index.xs, xm = element_index.xm,
           ys = element_index.ys, ym = element_index.ym;
  for (i=xs; i<xs+xm; i++) {
    for (j=ys;j<ys+ym; j++) {
      PetscReal zetaq[FEQuadrature::Nq];

      quadrature.computeTrialFunctionValues(i,j,dofmap,zeta_array,zetaq);
      const PetscInt ij = element_index.flatten(i,j);
      FEStoreNode *feS = &feStore[ij*FEQuadrature::Nq];
      PetscReal *dtauc_dzeta = m_dtauc_dzeta_store + ij*FEQuadrature::Nq;
      for (q=0; q<4; q++) {
        m_tauc_param.toTauc(zetaq[q],&(feS[q].tauc),dtauc_dzeta+q);
      }
    }
  }

  ierr = tauc->begin_access(); CHKERRQ(ierr);
  for (i=xs; i<xs+xm; i++) {
    for (j=ys;j<ys+ym; j++) {
      PetscReal tmp;
      m_tauc_param.toTauc(new_zeta(i,j),&tmp,NULL);
      (*tauc)(i,j) = tmp;
    }
  }
  
  ierr = tauc->end_access(); CHKERRQ(ierr);
  ierr = new_zeta.end_access();CHKERRQ(ierr);

  m_reassemble_T_matrix_needed = true;
  m_forward_F_needed = true;

  return 0;
}

int findNextFile(const char *basename)
{
  int n=0;
  while(true)
  {
    std::ostringstream os_nc;
    os_nc << basename << n << ".nc";
    bool does_not_exist = (access(os_nc.str().c_str(),W_OK) == -1) && (errno==ENOENT);
    if( does_not_exist )
    {
      std::ostringstream os_petsc;
      os_petsc << basename << n << ".petsc";
      does_not_exist = (access(os_petsc.str().c_str(),W_OK) == -1) && (errno==ENOENT);
    }

    if(does_not_exist) break;
    n++;
  }
  return n;
}

PetscErrorCode InvSSAForwardProblem::solveF_core()
{
  PetscErrorCode ierr;

  // FIXME (DAM 10/25/11): All this lousy code duplication with SSAFEM::solve.  Why
  // didn't I just call it when I wrote this in the first place???
  m_epsilon_ssa = config.get("epsilon_ssa");
  const PetscScalar DEFAULT_EPSILON_MULTIPLIER_SSA = 4.0;

  if(m_forward_F_needed)
  {
    // Set the SNES callbacks to call into our compute_local_function and compute_local_jacobian
    // methods via SSAFEFunction and SSAFEJ
    ierr = DMDASetLocalFunction(SSADA,(DMDALocalFunction1)SSAFEFunction);CHKERRQ(ierr);
    ierr = DMDASetLocalJacobian(SSADA,(DMDALocalFunction1)SSAFEJacobian);CHKERRQ(ierr);
    callback_data.da = SSADA;  callback_data.ssa = this;
    ierr = SNESSetDM(snes, SSADA); CHKERRQ(ierr);
    ierr = SNESSetFunction(snes, r,    SNESDAFormFunction,   &callback_data);CHKERRQ(ierr);
    ierr = SNESSetJacobian(snes, J, J, SNESDAComputeJacobian,&callback_data);CHKERRQ(ierr);

    SNESConvergedReason reason;
    while(true) {
    // Solve:
    ierr = SNESSolve(snes,NULL,SSAX);CHKERRQ(ierr);

      // See if it worked.
      ierr = SNESGetConvergedReason( snes, &reason); CHKERRQ(ierr);
      if(reason >= 0) {
        break;
      }

      ierr = verbPrintf(1,grid.com,
          "\nPISM WARNING:  SNESSolve() reports 'diverged'; reason = %d = '%s'\n",
          reason,SNESConvergedReasons[reason]); CHKERRQ(ierr);

      const char *savefile = "InvSSA_snesdivergederror";
      int index = findNextFile(savefile);

      std::ostringstream os_petscfile;
      os_petscfile << savefile << index << ".petsc";
      std::string petscfile = os_petscfile.str();
      const char *c_petscfile = petscfile.c_str();

      // char filename[PETSC_MAX_PATH_LEN] = "SSAFEM_snesdivergederror.petsc";
      ierr = verbPrintf(1,grid.com,
          "  writing linear system to PETSc binary file %s ...\n",c_petscfile); CHKERRQ(ierr);
      PetscViewer    viewer;
      ierr = PetscViewerBinaryOpen(grid.com, c_petscfile, FILE_MODE_WRITE, &viewer); CHKERRQ(ierr);
      ierr = MatView(J,viewer); CHKERRQ(ierr);
      ierr = VecView(r,viewer); CHKERRQ(ierr);
      ierr = VecView(SSAX,viewer); CHKERRQ(ierr);
      ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);




      std::ostringstream os_ncfile;
      os_ncfile << savefile << index << ".nc";
      std::string ncfile = os_ncfile.str();

      PIO pio(grid.com, grid.rank, grid.config.get_string("output_format"));
      pio.open(ncfile, PISM_WRITE);
      pio.def_time(grid.config.get_string("time_dimension_name"),
                   grid.config.get_string("calendar"),
                   grid.time->CF_units());
      pio.append_time(grid.config.get_string("time_dimension_name"),0.0);
      pio.close();

      mask->write(ncfile);
      thickness->write(ncfile);
      bed->write(ncfile);
      tauc->write(ncfile);
      enthalpy->write(ncfile);

      if(surface != NULL)
      {
        surface->write(ncfile);
      }
      if(driving_stress_x != NULL)
      {
        driving_stress_x->write(ncfile);
      }
      if(driving_stress_y != NULL)
      {
        driving_stress_y->write(ncfile);
      }


      SETERRQ1(grid.com, 1,
        "InvSSAForwardProblem solve failed to converge (SNES reason %s).\nSo we're giving up.\n\n", SNESConvergedReasons[reason]);


      if(m_epsilon_ssa <= 0.)
      {
        SETERRQ1(grid.com, 1,
          "InvSSAForwardProblem solve failed to converge (SNES reason %s).\nRegularization parameter epsilon_ssa = %f <=0, so we're giving up.\n\n", SNESConvergedReasons[reason]);
      }
      else
      {
        // Bump the regularization and try it again.
        m_epsilon_ssa *=  DEFAULT_EPSILON_MULTIPLIER_SSA;
      }
    }

    verbPrintf(3,grid.com,"SSAFEM solve converged (SNES reason %s)\n\n", SNESConvergedReasons[reason]);

    ierr =  DMGlobalToLocalBegin(SSADA, SSAX, INSERT_VALUES, m_VecU);  CHKERRQ(ierr);
    ierr =    DMGlobalToLocalEnd(SSADA, SSAX, INSERT_VALUES, m_VecU);  CHKERRQ(ierr);
    // ierr = DALocalToLocalBegin(SSADA, m_VecU, INSERT_VALUES, m_VecU);  CHKERRQ(ierr);
    // ierr = DALocalToLocalEnd(SSADA, m_VecU, INSERT_VALUES, m_VecU);  CHKERRQ(ierr);

    m_forward_F_needed = false;
  }
  return 0;
}

/* \brief Solves the nonlinear forward problem taking tauc to SSA velocities. */
PetscErrorCode InvSSAForwardProblem::solveF(IceModelVec2V &result)
{
  PetscErrorCode ierr;
  ierr = solveF_core(); CHKERRQ(ierr);

  ierr = result.copy_from(SSAX); CHKERRQ(ierr);
  ierr = result.beginGhostComm(); CHKERRQ(ierr);
  ierr = result.endGhostComm(); CHKERRQ(ierr);
  return 0;
}

/* \brief Solves the linearized forward problem */
PetscErrorCode InvSSAForwardProblem::solveT( IceModelVec2S &dtauc, IceModelVec2V &result)
{
  KSPConvergedReason  reason;
  PetscErrorCode ierr;
  PetscScalar **dtauc_a;
  PISMVector2 **vel;
  PISMVector2 **rhs;

  // We need SSA velocities to compute the linearization
  ierr = solveF_core(); CHKERRQ(ierr);

  // If tauc has been touched, we need to reassemble the linearzation matrix.
  if(m_reassemble_T_matrix_needed)
  {
    assemble_T_matrix();
  }

  // Assemble the right-hand side for the linearized forward problem.
  dtauc.get_array(dtauc_a);
  ierr = DMDAVecGetArray(SSADA,m_VecU,&vel); CHKERRQ(ierr);
  ierr = DMDAVecGetArray(SSADA,m_VecRHS2,&rhs); CHKERRQ(ierr);
  ierr = assemble_T_rhs(vel,dtauc_a,rhs); CHKERRQ(ierr);
  dtauc.end_access();
  ierr = DMDAVecRestoreArray(SSADA, m_VecU, &vel); CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(SSADA, m_VecRHS2, &rhs); CHKERRQ(ierr);


  // call PETSc to solve linear system by iterative method.
  ierr = KSPSetOperators(m_KSP, m_MatA, m_MatA, SAME_NONZERO_PATTERN); CHKERRQ(ierr);
  ierr = KSPSolve(m_KSP, m_VecRHS2, m_VecZ2); CHKERRQ(ierr); // SOLVE

  ierr = KSPGetConvergedReason(m_KSP, &reason); CHKERRQ(ierr);
  if (reason < 0) {
    SETERRQ1(grid.com, 1,
      "InvSSAForwardProblem::solveT solve failed to converge (KSP reason %s)\n\n", KSPConvergedReasons[reason]);
  }
  else
  {
    verbPrintf(4,grid.com,"InvSSAForwardProblem::solveT converged (KSP reason %s)\n", KSPConvergedReasons[reason] );
  }

  // Extract the solution and communicate.
  ierr = result.copy_from(m_VecZ2); CHKERRQ(ierr);
  ierr = result.beginGhostComm(); CHKERRQ(ierr);
  ierr = result.endGhostComm(); CHKERRQ(ierr);

  return 0;
}

//! \brief Solves the adjoint of the linearized forward problem.
PetscErrorCode InvSSAForwardProblem::solveTStar( IceModelVec2V &residual, IceModelVec2S &result)
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
  ierr = DMDAVecGetArray(SSADA,m_VecRHS2,&RHS2); CHKERRQ(ierr);
  ierr = assemble_TStarA_rhs(R,RHS2); CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(SSADA, m_VecRHS2, &RHS2); CHKERRQ(ierr);
  residual.end_access();

  // call PETSc to solve linear system by iterative method.
  ierr = KSPSetOperators(m_KSP, m_MatA, m_MatA, SAME_NONZERO_PATTERN); CHKERRQ(ierr);
  ierr = KSPSolve(m_KSP, m_VecRHS2, m_VecZ2); CHKERRQ(ierr); // SOLVE
  ierr = KSPGetConvergedReason(m_KSP, &reason); CHKERRQ(ierr);
  if (reason < 0) {
    SETERRQ1(grid.com, 1,
      "InvSSAForwardProblem::solveTStarA solve failed to converge (KSP reason %s)\n\n", KSPConvergedReasons[reason]);
  }
  else
  {
    verbPrintf(4,grid.com,"InvSSAForwardProblem::solveTStarA converged (KSP reason %s)\n", KSPConvergedReasons[reason] );
  }


  ierr = DMGlobalToLocalBegin(SSADA, m_VecZ2, INSERT_VALUES, m_VecZ);  CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(SSADA, m_VecZ2, INSERT_VALUES, m_VecZ);  CHKERRQ(ierr);

  // Assemble the right-hand side for the second step.
  ierr = DMDAVecGetArray(SSADA,m_VecZ,&Z); CHKERRQ(ierr);
  ierr = DMDAVecGetArray(SSADA,m_VecU,&U); CHKERRQ(ierr);
  ierr = DMDAVecGetArray(grid.da2,m_VecRHS,&RHS); CHKERRQ(ierr);
  ierr = assemble_TStarB_rhs(Z,U,RHS); CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(grid.da2,m_VecRHS,&RHS); CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(SSADA,m_VecU,&U); CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(SSADA,m_VecZ,&Z); CHKERRQ(ierr);


  // call PETSc to solve linear system by iterative method.
  ierr = KSPSetOperators(m_KSP_B, m_MatB, m_MatB, SAME_NONZERO_PATTERN); CHKERRQ(ierr);
  ierr = KSPSolve(m_KSP_B, m_VecRHS, m_VecV); CHKERRQ(ierr); // SOLVE
  ierr = KSPGetConvergedReason(m_KSP, &reason); CHKERRQ(ierr);
  if (reason < 0) {
    SETERRQ1(grid.com, 1,
      "InvSSAForwardProblem::solveTStarB solve failed to converge (KSP reason %s)\n\n", KSPConvergedReasons[reason]);
  }
  else  {
    verbPrintf(4,grid.com,"InvSSAForwardProblem::solveTStarB converged (KSP reason %s)\n", KSPConvergedReasons[reason] );
  }

  // Extract the solution and communicate.
  ierr = result.copy_from(m_VecV); CHKERRQ(ierr);
  ierr = result.beginGhostComm(); CHKERRQ(ierr);
  ierr = result.endGhostComm(); CHKERRQ(ierr);

  return 0;
}

//! \brief Builds the matrix for the linearized forward PDE.
PetscErrorCode InvSSAForwardProblem::assemble_T_matrix()
{
  PISMVector2 **vel;
  PetscErrorCode ierr;

  ierr = DMDAVecGetArray(SSADA,m_VecU,&vel); CHKERRQ(ierr);

  DMDALocalInfo *info = NULL;
  ierr = compute_local_jacobian(info,const_cast<const PISMVector2**>(vel),m_MatA);

  ierr = DMDAVecRestoreArray(SSADA, m_VecU, &vel); CHKERRQ(ierr);
  return 0;
}

PetscErrorCode InvSSAForwardProblem::assemble_DomainNorm_matrix()
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

  // Flag at each local node that determines if this value of zeta is unchanging.
  PetscReal **zeta_fixed_mask;
  PetscReal local_zeta_fixed_mask[FEQuadrature::Nk];
  for(int k=0; k<FEQuadrature::Nk;k++) {
    local_zeta_fixed_mask[k] = 0.;
  }
  // Start access of Dirichlet data, if present.
  if (m_zeta_fixed_locations ) {
    ierr = m_zeta_fixed_locations->get_array(zeta_fixed_mask);CHKERRQ(ierr);
  }

  PetscReal cH1 = config.get("inv_ssa_domain_h1_coeff");
  PetscReal cL2 = config.get("inv_ssa_domain_l2_coeff");

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

      // These values now need to be adjusted if some nodes in the element have
      // Dirichlet data.
      if(m_zeta_fixed_locations) {
        dofmap.extractLocalDOFs(i,j,zeta_fixed_mask,local_zeta_fixed_mask);
      }

      // Build the element-local Jacobian.
      ierr = PetscMemzero(K,sizeof(K));CHKERRQ(ierr);
      for (PetscInt q=0; q<FEQuadrature::Nq; q++) {
        for (PetscInt k=0; k<4; k++) {   // Test functions
          for (PetscInt l=0; l<4; l++) { // Trial functions
            if( (local_zeta_fixed_mask[k]==0.) && (local_zeta_fixed_mask[l]==0.) ) {
              const FEFunctionGerm &test_qk=test[q][k];
              const FEFunctionGerm &test_ql=test[q][l];
              K[k][l]     += JxW[q]*(cL2*test_qk.val*test_ql.val
              +  cH1*(test_qk.dx*test_ql.dx + test_qk.dy*test_ql.dy) );
            } // zeta fixed
          } // l
        } // k
      } // q
      ierr = dofmap.addLocalJacobianBlock(&K[0][0],m_MatB);
    } // j
  } // i

  // Until now, the rows and columns corresponding to fixed zeta values have not been set.  We now
  // put an identity block in for these unknowns.  
  if (m_zeta_fixed_locations) {
    for (i=grid.xs; i<grid.xs+grid.xm; i++) {
      for (j=grid.ys; j<grid.ys+grid.ym; j++) {
        if (m_zeta_fixed_locations->as_int(i,j) == 1) {
          const PetscReal ident[4] = {dirichletScale,0,0,dirichletScale};
          MatStencil row;
          // FIXME: Transpose shows up here!
          row.j = i; row.i = j;
          ierr = MatSetValuesBlockedStencil(m_MatB,1,&row,1,&row,ident,ADD_VALUES);CHKERRQ(ierr);
        }
      }
    }
    ierr = m_zeta_fixed_locations->end_access(); CHKERRQ(ierr);
  }

  ierr = MatAssemblyBegin(m_MatB,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(m_MatB,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

  return 0;
}

PetscErrorCode InvSSAForwardProblem::domainIP(Vec a, Vec b, PetscScalar *OUTPUT)
{
  PetscErrorCode ierr;
  PetscScalar **A, **B;
  ierr = DMDAVecGetArray(grid.da2,a,&A); CHKERRQ(ierr);
  ierr = DMDAVecGetArray(grid.da2,b,&B); CHKERRQ(ierr);
  ierr = domainIP_core(A,B,OUTPUT);
  ierr = DMDAVecRestoreArray(grid.da2,a,&A); CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(grid.da2,b,&B); CHKERRQ(ierr);
  return 0;
}
PetscErrorCode InvSSAForwardProblem::domainIP(IceModelVec2S &a, IceModelVec2S &b, PetscScalar *OUTPUT)
{
  PetscErrorCode ierr;
  PetscReal **A, **B;
  ierr = a.get_array(A); CHKERRQ(ierr);
  ierr = b.get_array(B); CHKERRQ(ierr);
  ierr = domainIP_core(A,B,OUTPUT);
  ierr = a.end_access(); CHKERRQ(ierr);
  ierr = b.end_access(); CHKERRQ(ierr);
  return 0;
}

PetscErrorCode InvSSAForwardProblem::rangeIP(Vec a, Vec b, PetscScalar *OUTPUT)
{
  PetscErrorCode ierr;
  PISMVector2 **A, **B;
  ierr = DMDAVecGetArray(SSADA,a,&A); CHKERRQ(ierr);
  ierr = DMDAVecGetArray(SSADA,b,&B); CHKERRQ(ierr);
  ierr = rangeIP_core(A,B,OUTPUT);
  ierr = DMDAVecRestoreArray(SSADA,a,&A); CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(SSADA,b,&B); CHKERRQ(ierr);
  return 0;
}

PetscErrorCode InvSSAForwardProblem::rangeIP(IceModelVec2V &a, IceModelVec2V &b, PetscScalar *OUTPUT)
{
  PetscErrorCode ierr;
  PISMVector2 **A, **B;
  ierr = a.get_array(A); CHKERRQ(ierr);
  ierr = b.get_array(B); CHKERRQ(ierr);
  ierr = rangeIP_core(A,B,OUTPUT);
  ierr = a.end_access(); CHKERRQ(ierr);
  ierr = b.end_access(); CHKERRQ(ierr);
  return 0;
}

PetscErrorCode InvSSAForwardProblem::domainIP_core(PetscReal **A, PetscReal**B, PetscScalar *OUTPUT)
{

  // Just L2 matrix for now.
  PetscInt         i,j;
  PetscErrorCode   ierr;

  // The value of the inner product.
  PetscReal IP = 0;

  PetscReal a[FEQuadrature::Nq], b[FEQuadrature::Nq],
    ax[FEQuadrature::Nq], bx[FEQuadrature::Nq],
    ay[FEQuadrature::Nq], by[FEQuadrature::Nq];

  // Jacobian times weights for quadrature.
  PetscScalar JxW[FEQuadrature::Nq];
  quadrature.getWeightedJacobian(JxW);

  // Flag at each local node that determines if this value of zeta is unchanging.
  PetscReal **zeta_fixed_mask;
  PetscReal local_zeta_fixed_mask[FEQuadrature::Nk];
  for(int k=0; k<FEQuadrature::Nk;k++) {
    local_zeta_fixed_mask[k] = 0.;
  }
  // Start access of Dirichlet data, if present.
  if (m_zeta_fixed_locations ) {
    ierr = m_zeta_fixed_locations->get_array(zeta_fixed_mask);CHKERRQ(ierr);
  }

  PetscReal cH1 = config.get("inv_ssa_domain_h1_coeff");
  PetscReal cL2 = config.get("inv_ssa_domain_l2_coeff");

  // Loop through all LOCAL elements.
  PetscInt xs = element_index.lxs, xm = element_index.lxm,
           ys = element_index.lys, ym = element_index.lym;
  for (i=xs; i<xs+xm; i++) {
    for (j=ys; j<ys+ym; j++) {


      // Storage for element-local solution and residuals.
      PetscReal     a_local[FEQuadrature::Nk], b_local[FEQuadrature::Nk];

      // Obtain the value of the solution at the nodes adjacent to the element.
      dofmap.extractLocalDOFs(i,j,A,a_local);
      dofmap.extractLocalDOFs(i,j,B,b_local);
      if(m_zeta_fixed_locations) {
        dofmap.extractLocalDOFs(i,j,zeta_fixed_mask,local_zeta_fixed_mask);
        for (PetscInt k=0; k<FEQuadrature::Nq; k++) {
          if (PismIntMask(local_zeta_fixed_mask[k]) != 0) { // Dirichlet node
            a_local[k] = 0;
            b_local[k] = 0;
          }
        }
      }
      // Compute the solution values and symmetric gradient at the quadrature points.
      quadrature.computeTrialFunctionValues(a_local,a,ax,ay);
      quadrature.computeTrialFunctionValues(b_local,b,bx,by);

      for (PetscInt q=0; q<FEQuadrature::Nq; q++) {
        IP += JxW[q]*(cL2*a[q]*b[q]+ cH1*(ax[q]*bx[q]+ay[q]*by[q]));
      } // q
    } // j
  } // i

  // End access of Dirichlet data, if present.
  if (m_zeta_fixed_locations ) {
    ierr = m_zeta_fixed_locations->end_access(); CHKERRQ(ierr);
  }

  ierr = PISMGlobalSum(&IP, OUTPUT, grid.com); CHKERRQ(ierr);
  return 0;
}

PetscErrorCode InvSSAForwardProblem::rangeIP_core(PISMVector2 **A, PISMVector2**B, PetscScalar *OUTPUT)
{
  PetscInt         i,j;
  PetscErrorCode   ierr;

  // The value of the inner product on the local part of the domain.
  PetscReal IP = 0;

  PISMVector2 a[FEQuadrature::Nq], b[FEQuadrature::Nq];

  PetscReal **W;
  PetscReal misfit_weight[FEQuadrature::Nq];
  if(m_misfit_weight!=NULL) {
    ierr = m_misfit_weight->get_array(W);CHKERRQ(ierr);
  } else {
    for(int q=0;q<FEQuadrature::Nq;q++) {
      misfit_weight[q]=1;
    }
  }

  // Jacobian times weights for quadrature.
  PetscScalar JxW[FEQuadrature::Nq];
  quadrature.getWeightedJacobian(JxW);

  if(m_misfit_element_mask!=NULL) {
    ierr = m_misfit_element_mask->begin_access();CHKERRQ(ierr);    
  }

  // Loop through all LOCAL elements.
  PetscInt xs = element_index.lxs, xm = element_index.lxm,
           ys = element_index.lys, ym = element_index.lym;
  for (i=xs; i<xs+xm; i++) {
    for (j=ys; j<ys+ym; j++) {
      if(m_misfit_element_mask != NULL) {
        if( m_misfit_element_mask->as_int(i,j) != 1) {
          continue;
        }
      }
      PISMVector2 tmp[FEQuadrature::Nq];
      dofmap.extractLocalDOFs(i,j,A,tmp);
      quadrature.computeTrialFunctionValues(tmp,a);
      dofmap.extractLocalDOFs(i,j,B,tmp);
      quadrature.computeTrialFunctionValues(tmp,b);

      if(m_misfit_weight != NULL) {
        quadrature.computeTrialFunctionValues(i,j,dofmap,W,misfit_weight);
      }
      for (PetscInt q=0; q<FEQuadrature::Nq; q++) {
        IP += JxW[q]*(a[q].u*b[q].u + a[q].v*b[q].v)*misfit_weight[q];
      } // q
    } // j
  } // i

  if(m_misfit_weight!=NULL) {
    ierr = m_misfit_weight->end_access();CHKERRQ(ierr);
  }
  if(m_misfit_element_mask!=NULL) {
    ierr = m_misfit_element_mask->end_access();CHKERRQ(ierr);    
  }


  IP /= m_range_l2_area;
  ierr = PISMGlobalSum(&IP, OUTPUT, grid.com); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode InvSSAForwardProblem::compute_range_l2_area(PetscScalar *OUTPUT)
{
  PetscInt         i,j;
  PetscErrorCode   ierr;

  // The value of the inner product on the local part of the domain.
  PetscReal IP = 0;

  PetscReal **W;
  PetscReal misfit_weight[FEQuadrature::Nq];
  if(m_misfit_weight!=NULL) {
    ierr = m_misfit_weight->get_array(W);CHKERRQ(ierr);
  } else {
    for(int q=0;q<FEQuadrature::Nq;q++) {
      misfit_weight[q]=1.;
    }
  }

  // Jacobian times weights for quadrature.
  PetscScalar JxW[FEQuadrature::Nq];
  quadrature.getWeightedJacobian(JxW);

  // Loop through all LOCAL elements.
  PetscInt xs = element_index.lxs, xm = element_index.lxm,
           ys = element_index.lys, ym = element_index.lym;
  for (i=xs; i<xs+xm; i++) {
    for (j=ys; j<ys+ym; j++) {
      if(m_misfit_weight != NULL) {
        quadrature.computeTrialFunctionValues(i,j,dofmap,W,misfit_weight);
      }
      for (PetscInt q=0; q<FEQuadrature::Nq; q++) {
        IP += JxW[q]*misfit_weight[q];
      } // q
    } // j
  } // i

  if(m_misfit_weight!=NULL) {
    ierr = m_misfit_weight->end_access();CHKERRQ(ierr);
  }

  ierr = PISMGlobalSum(&IP, OUTPUT, grid.com); CHKERRQ(ierr);

  return 0;
}


PetscErrorCode InvSSAForwardProblem::assemble_T_rhs( PISMVector2 **gvel, PetscReal **gdtauc, PISMVector2 **grhs)
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
  Mask M;

  // Iterate over the elements.
  PetscInt xs = element_index.xs, xm = element_index.xm,
           ys = element_index.ys, ym = element_index.ym;
  for (i=xs; i<xs+xm; i++) {
    for (j=ys; j<ys+ym; j++) {
      // Storage for element-local data
      PISMVector2     y[4];

      // Index into coefficient storage in feStore
      const PetscInt ij = element_index.flatten(i,j);

      // Initialize the map from global to local degrees of freedom for this element.
      dofmap.reset(i,j,grid);

      // These values now need to be adjusted if some nodes in the element have
      // Dirichlet data.
      if(bc_locations && vel_bc) {
        dofmap.extractLocalDOFs(i,j,bc_mask,local_bc_mask);
        for (k=0; k<4; k++) {
          if (PismIntMask(local_bc_mask[k]) == 1) { // Dirichlet node
            // Mark any kind of Dirichlet node as not to be touched
            dofmap.markRowInvalid(k);
            dofmap.markColInvalid(k);
          }
        }
      }

      // Obtain the value of the solution at the adjacent nodes to the element.
      quadrature.computeTrialFunctionValues(i,j,dofmap,gvel,u);
      quadrature.computeTrialFunctionValues(i,j,dofmap,gdtauc,dtauc);

      // Zero out the element-local residual in prep for updating it.
      for(k=0;k<FEQuadrature::Nk;k++){
        y[k].u = 0; y[k].v = 0;
      }

      for (q=0; q<FEQuadrature::Nq; q++) {     // loop over quadrature points on this element.

        // Coefficients and weights for this quadrature point.
        const FEStoreNode *feS = &feStore[ij*FEQuadrature::Nq+q];
        const PetscReal    jw  = JxW[q];

        // Determine dbeta/dp at the quadrature point
        PetscReal dbeta_dp = 0;
        if( M.grounded_ice(feS->mask) ) {
          PetscReal dbeta_dtauc = basal.drag(dtauc[q],u[q].u,u[q].v);
          PetscReal dtauc_dzeta = m_dtauc_dzeta_store[ij*FEQuadrature::Nq+q];
          dbeta_dp = dbeta_dtauc*dtauc_dzeta;
        }

        for(k=0; k<FEQuadrature::Nk;k++) {  // loop over the test functions.
          const FEFunctionGerm &testqk = test[q][k];
          //Note the -= (not +=) in the following lines.
          y[k].u -= jw*dbeta_dp*testqk.val*u[q].u;
          y[k].v -= jw*dbeta_dp*testqk.val*u[q].v;
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
        if (bc_locations->as_int(i,j) == 1) {
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

PetscErrorCode InvSSAForwardProblem::assemble_TStarA_rhs( PISMVector2 **R, PISMVector2 **RHS)
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

  if(m_misfit_element_mask!=NULL) {
    ierr = m_misfit_element_mask->begin_access();CHKERRQ(ierr);    
  }

  // Jacobian times weights for quadrature.
  PetscScalar JxW[FEQuadrature::Nq];
  quadrature.getWeightedJacobian(JxW);

  // Flags for each vertex in an element that determine if explicit Dirichlet data has
  // been set.
  PetscReal local_bc_mask[FEQuadrature::Nk];

  // Storage for R at quadrature points.
  PISMVector2 res[FEQuadrature::Nq];

  // An Nq by Nk array of test function values.
  const FEFunctionGerm (*test)[FEQuadrature::Nk] = quadrature.testFunctionValues();

  PetscReal **W;
  PetscReal misfit_weight[FEQuadrature::Nq];
  if(m_misfit_weight!=NULL) {
    ierr = m_misfit_weight->get_array(W);CHKERRQ(ierr);
  } else {
    for(q=0;q<FEQuadrature::Nq;q++) {
      misfit_weight[q]=1;
    }
  }

  // Iterate over the elements.
  PetscInt xs = element_index.xs, xm = element_index.xm,
           ys = element_index.ys, ym = element_index.ym;
  for (i=xs; i<xs+xm; i++) {
    for (j=ys; j<ys+ym; j++) {

      if(m_misfit_element_mask != NULL) {
        if( m_misfit_element_mask->as_int(i,j) != 1) {
          continue;
        }
      }

      // Storage for element-local data
      PISMVector2 y[FEQuadrature::Nk];

      // Initialize the map from global to local degrees of freedom for this element.
      dofmap.reset(i,j,grid);


      // These values now need to be adjusted if some nodes in the element have
      // Dirichlet data.
      if(bc_locations && vel_bc) {
        dofmap.extractLocalDOFs(i,j,bc_mask,local_bc_mask);
        for (k=0; k<4; k++) {
          if (PismIntMask(local_bc_mask[k]) == 1) { // Dirichlet node
            // Mark any kind of Dirichlet node as not to be touched
            dofmap.markRowInvalid(k);
            dofmap.markColInvalid(k);
          }
        }
      }

      // Obtain the value of the solution at the adjacent nodes to the element.
      quadrature.computeTrialFunctionValues(i,j,dofmap,R,res);

      // Zero out the element-local residual in prep for updating it.
      for(k=0;k<FEQuadrature::Nk;k++){
        y[k].u = 0; y[k].v = 0;
      }

      if(m_misfit_weight != NULL) {
        quadrature.computeTrialFunctionValues(i,j,dofmap,W,misfit_weight);
      }

      for (q=0; q<FEQuadrature::Nq; q++) {     // loop over quadrature points on this element.
        // Coefficients and weights for this quadrature point.
        const PetscReal    jw  = JxW[q]/m_range_l2_area;
        for(k=0; k<FEQuadrature::Nk;k++) {  // loop over the test functions.
          const FEFunctionGerm &testqk = test[q][k];
          y[k].u += jw*testqk.val*res[q].u*misfit_weight[q];
          y[k].v += jw*testqk.val*res[q].v*misfit_weight[q];
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
        if (bc_locations->as_int(i,j) == 1) {
          // Enforce explicit homogeneous dirichlet data.
          RHS[i][j].u = 0;
          RHS[i][j].v = 0;
        }
      }
    }
    ierr = bc_locations->end_access();CHKERRQ(ierr);
  }

  if(m_misfit_weight!=NULL) {
    ierr = m_misfit_weight->end_access();CHKERRQ(ierr);
  }
  if(m_misfit_element_mask!=NULL) {
    ierr = m_misfit_element_mask->end_access();CHKERRQ(ierr);    
  }

  return 0;
}


PetscErrorCode InvSSAForwardProblem::assemble_TStarB_rhs( PISMVector2 **Z,
                                                    PISMVector2 **U,
                                                    PetscScalar **RHS)
{
  PetscInt         i,j,k,q;

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

  Mask M;

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
        // Coefficients and weights for this quadrature point.
        const FEStoreNode *feS = &feStore[ij*FEQuadrature::Nq+q];
        const PetscReal    jw  = JxW[q];

        // Determine "beta" at the quadrature point
        PetscReal notquitebeta = 0;
        if( M.grounded_ice(feS->mask) ) {
          notquitebeta = basal.drag(1.,u[q].u,u[q].v);
          PetscReal dtauc_dzeta = m_dtauc_dzeta_store[ij*FEQuadrature::Nq+q];
          notquitebeta*=dtauc_dzeta;
        }

        for(k=0; k<FEQuadrature::Nk;k++) {  // loop over the test functions.
          const FEFunctionGerm &testqk = test[q][k];
          y[k] -= jw*notquitebeta*testqk.val*(u[q].u*z[q].u+u[q].v*z[q].v);
        }
      } // q

      dofmap.addLocalResidualBlock(y,RHS);
    } // j-loop
  } // i-loop

  // In fact we want the right-hand side to be zeros in locations where zeta
  // is fixed. So we adjust the right-hand side now.
  if (m_zeta_fixed_locations) {
    // Enforce Dirichlet conditions strongly
    PetscErrorCode ierr = m_zeta_fixed_locations->begin_access();CHKERRQ(ierr);
    for (i=grid.xs; i<grid.xs+grid.xm; i++) {
      for (j=grid.ys; j<grid.ys+grid.ym; j++) {
        if (m_zeta_fixed_locations->as_int(i,j) == 1) {
          RHS[i][j] = 0;
        }
      }
    }
    ierr = m_zeta_fixed_locations->end_access();CHKERRQ(ierr);
  }

  return 0;

}

