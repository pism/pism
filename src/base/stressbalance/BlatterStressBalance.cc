// Copyright (C) 2010-2011 Jed Brown, Ed Bueler and Constantine Khroulev
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

#include "BlatterStressBalance.hh"
#include "PISMOcean.hh"
#include "IceGrid.hh"

BlatterStressBalance::BlatterStressBalance(IceGrid &g, PISMOceanModel *ocean_model, const NCConfigVariable &conf)
  : PISMStressBalance(g,NULL,NULL,ocean_model,conf)
{
  if (allocate_blatter() != 0) {
    PetscPrintf(grid.com, "FATAL ERROR: BlatterStressBalance allocation failed.\n");
    PISMEnd();
  }
}

BlatterStressBalance::~BlatterStressBalance()
{
  if (deallocate_blatter() != 0) {
    PetscPrintf(grid.com, "FATAL ERROR: BlatterStressBalance deallocation failed.\n");
    PISMEnd();
  }
}

PetscErrorCode BlatterStressBalance::allocate_blatter() {
  PetscErrorCode ierr;
  PetscInt i;

  ierr = THICreate(grid.com,&thi);CHKERRQ(ierr);
  ierr = DMMGCreate(PETSC_COMM_WORLD,thi->nlevels,thi,&dmmg);CHKERRQ(ierr);
  {
    DA da;
    int P = 2;
    ierr = PetscOptionsBegin(grid.com,NULL,"Grid resolution options","");CHKERRQ(ierr);
    {
      if (thi->coarse2d) {
        ierr = PetscOptionsInt("-blatter_zlevels",
                               "Number of elements in z-direction on fine level","",thi->zlevels,&thi->zlevels,NULL); CHKERRQ(ierr);
      } else {
        ierr = PetscOptionsInt("-blatter_P",
                               "Number of elements in z-direction on coarse level","",P,&P,NULL); CHKERRQ(ierr);
      }
    }
    ierr = PetscOptionsEnd();CHKERRQ(ierr);
    if (thi->coarse2d) {
      ierr = DACreate2d(grid.com,DA_XYPERIODIC,DA_STENCIL_BOX,
                        grid.My,grid.Mx,
                        PETSC_DETERMINE,PETSC_DETERMINE,
                        sizeof(Node)/sizeof(PetscScalar),1,0,0,&da);CHKERRQ(ierr);
      da->ops->refinehierarchy  = DARefineHierarchy_THI;
      da->ops->getinterpolation = DAGetInterpolation_THI;
      ierr = PetscObjectCompose((PetscObject)da,"THI",(PetscObject)thi);CHKERRQ(ierr);
    } else {
      ierr = DACreate3d(grid.com,DA_YZPERIODIC,DA_STENCIL_BOX,
                        P,grid.My,grid.Mx,
                        1,PETSC_DETERMINE,PETSC_DETERMINE,
                        sizeof(Node)/sizeof(PetscScalar),1,0,0,0,&da);CHKERRQ(ierr);
    }
    ierr = DASetFieldName(da,0,"x-velocity");CHKERRQ(ierr);
    ierr = DASetFieldName(da,1,"y-velocity");CHKERRQ(ierr);
    ierr = DMMGSetDM(dmmg,(DM)da);CHKERRQ(ierr);
    ierr = DADestroy(da);CHKERRQ(ierr);
  }

  if (thi->tridiagonal) {
    (DMMGGetDA(dmmg))->ops->getmatrix = DAGetMatrix_THI_Tridiagonal;
  }

  {
    /* Use the user-defined matrix type on all but the coarse level */
    ierr = DMMGSetMatType(dmmg,thi->mattype);CHKERRQ(ierr);
    /* PCREDUNDANT only works with AIJ, and so do the third-party direct solvers.  So when running in parallel, we can't
     * use the faster (S)BAIJ formats on the coarse level. */
    ierr = PetscFree(dmmg[0]->mtype);CHKERRQ(ierr);
    ierr = PetscStrallocpy(MATAIJ,&dmmg[0]->mtype);CHKERRQ(ierr);
  }

  /* Spectacularly ugly API, our function evaluation provides ghost values */
  ierr = PetscOptionsSetValue("-dmmg_form_function_ghost","1");CHKERRQ(ierr);

  ierr = DMMGSetSNESLocal(dmmg,THIFunctionLocal,THIJacobianLocal_3D_Full,0,0);CHKERRQ(ierr);

  if (thi->tridiagonal) {
    ierr = DASetLocalJacobian(DMMGGetDA(dmmg),(DALocalFunction1)THIJacobianLocal_3D_Tridiagonal);CHKERRQ(ierr);
  }

  if (thi->coarse2d) {
    for (i=0; i<DMMGGetLevels(dmmg)-1; i++) {
      ierr = DASetLocalJacobian((DA)dmmg[i]->dm,(DALocalFunction1)THIJacobianLocal_2D);CHKERRQ(ierr);
    }
  }

  for (i=0; i<DMMGGetLevels(dmmg); i++) {
    /* This option is only valid for the SBAIJ format.  The matrices we assemble are symmetric, but the SBAIJ assembly
     * functions will complain if we provide lower-triangular entries without setting this option. */
    Mat B = dmmg[i]->B;
    PetscTruth flg1,flg2;
    ierr = PetscTypeCompare((PetscObject)B,MATSEQSBAIJ,&flg1);CHKERRQ(ierr);
    ierr = PetscTypeCompare((PetscObject)B,MATMPISBAIJ,&flg2);CHKERRQ(ierr);
    if (flg1 || flg2) {
      ierr = MatSetOption(B,MAT_IGNORE_LOWER_TRIANGULAR,PETSC_TRUE);CHKERRQ(ierr);
    }
  }

  ierr = MatSetOptionsPrefix(DMMGGetB(dmmg),"thi_");CHKERRQ(ierr);
  ierr = DMMGSetFromOptions(dmmg);CHKERRQ(ierr);
  ierr = THISetDMMG(thi,dmmg);CHKERRQ(ierr);

  ierr = DMMGSetInitialGuess(dmmg,THIInitial);CHKERRQ(ierr);

  return 0;
}

PetscErrorCode BlatterStressBalance::deallocate_blatter() {
  PetscErrorCode ierr;

  ierr = DMMGDestroy(dmmg);CHKERRQ(ierr);
  ierr = THIDestroy(thi);CHKERRQ(ierr);

  return 0;
}


PetscErrorCode BlatterStressBalance::init(PISMVars &vars) {
  //PetscErrorCode ierr;

  variables = &vars;

  // FIXME: nontrivial initialization needed 

  return 0;
}

PetscErrorCode BlatterStressBalance::update(bool fast) {
  PetscErrorCode ierr;

  if (fast) {
    ierr = verbPrintf(1,grid.com,
       "PISM ERROR:  'fast' mode not meaningful for BlatterStressBalance\n"
       "  ENDING ...\n\n"); CHKERRQ(ierr);
    PISMEnd();
  }

  if (ocean) {
    PetscReal sea_level;
    ierr = ocean->sea_level_elevation(sea_level); CHKERRQ(ierr);
    // FIXME: blatter object should know about sea level
  }

  ierr = DMMGSolve(dmmg);CHKERRQ(ierr);
  ierr = THISolveStatistics(thi,dmmg,0,"Full");CHKERRQ(ierr);

  // Transfer solution from the FEM mesh to the regular grid used in the rest
  // of PISM and compute the vertically-averaged velocity.
  ierr = mesh_to_regular_grid(); CHKERRQ(ierr);
  
  ierr = compute_vertical_velocity(&u, &v, basal_melt_rate, w); CHKERRQ(ierr); 

  return 0;
}

PetscErrorCode BlatterStressBalance::mesh_to_regular_grid() {
  PetscErrorCode ierr;

  return 0;
}

PetscErrorCode BlatterStressBalance::regular_grid_to_mesh() {
  PetscErrorCode ierr;

  return 0;
}



PetscErrorCode BlatterStressBalance::get_max_2d_velocity(PetscReal &maxu, PetscReal &maxv) {
  // compute in mesh_to_regular_grid()
  return 0;
}

PetscErrorCode BlatterStressBalance::get_3d_velocity(IceModelVec3* &u_out, IceModelVec3* &v_out, IceModelVec3* &w_out) {
  u_out = &u;
  v_out = &v;
  w_out = &w;
  return 0;
}

PetscErrorCode BlatterStressBalance::get_max_3d_velocity(PetscReal &maxu, PetscReal &maxv, PetscReal &maxw) {
  // FIXME: this is probably not right *NOR* efficient.
  // revise for correctness by looking at FEM nodal values of horizontal velocity
  //   grid values on IceModel grid for vertical velocity
  // revise for efficiency if: finding this max could be done earlier or as part of update() or
  //   with less communication or better memory locality

  PetscErrorCode ierr;
  
  ierr = u.begin_access(); CHKERRQ(ierr);
  ierr = v.begin_access(); CHKERRQ(ierr);
  ierr = w.begin_access(); CHKERRQ(ierr);
  PetscReal my_umax = 0, my_vmax = 0, my_wmax = 0;
  PetscReal *ucol,*vcol,*wcol;
  for (PetscInt   i = grid.xs; i < grid.xs+grid.xm; ++i) {
    for (PetscInt j = grid.ys; j < grid.ys+grid.ym; ++j) {
      ierr = u.getInternalColumn(i, j, &ucol); CHKERRQ(ierr);
      ierr = v.getInternalColumn(i, j, &vcol); CHKERRQ(ierr);
      ierr = w.getInternalColumn(i, j, &wcol); CHKERRQ(ierr);
      for (PetscInt k = 0; k < grid.Mz; ++k) {
        my_umax = PetscMax(my_umax, PetscAbs(ucol[k]));
        my_vmax = PetscMax(my_vmax, PetscAbs(vcol[k]));
        my_wmax = PetscMax(my_wmax, PetscAbs(wcol[k]));
      }
    }
  }
  ierr = w.end_access(); CHKERRQ(ierr);
  ierr = v.end_access(); CHKERRQ(ierr);
  ierr = u.end_access(); CHKERRQ(ierr);  

  ierr = PetscGlobalMax(&my_umax, &maxu, grid.com); CHKERRQ(ierr);
  ierr = PetscGlobalMax(&my_vmax, &maxv, grid.com); CHKERRQ(ierr);
  ierr = PetscGlobalMax(&my_wmax, &maxw, grid.com); CHKERRQ(ierr);

  return 0;
}


PetscErrorCode BlatterStressBalance::extend_the_grid(PetscInt old_Mz) {
  PetscErrorCode ierr;
  ierr = u.extend_vertically(old_Mz, 0.0); CHKERRQ(ierr);
  ierr = v.extend_vertically(old_Mz, 0.0); CHKERRQ(ierr);
  ierr = w.extend_vertically(old_Mz, 0.0); CHKERRQ(ierr);
  return 0;
}

