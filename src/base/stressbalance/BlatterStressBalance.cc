// Copyright (C) 2010-2011 Jed Brown and Ed Bueler
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
#include "THI.c"
#include "PISMOcean.hh"
#include "IceGrid.hh"

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

  ierr = verbPrintf(1,grid.com,
    "\n\nJed's THI running now: ignoring PISM context! ...\n\n"); CHKERRQ(ierr);

  MPI_Comm       comm=grid.com;
  DMMG           *dmmg;
  THI            thi;
  PetscInt       i;
  PetscLogStage  stages[3];
  PetscTruth     repeat_fine_solve = PETSC_FALSE;

  /* We define two stages.  The first includes all setup costs and solves from a naive initial guess.  The second solve
  * is more indicative of what might occur during time-stepping.  The initial guess is interpolated from the next
  * coarser (as in the last step of grid sequencing), and so requires fewer Newton steps. */
  ierr = PetscOptionsGetTruth(NULL,"-blatter_repeat_fine_solve",&repeat_fine_solve,NULL);CHKERRQ(ierr);
  ierr = PetscLogStageRegister("Full solve",&stages[0]);CHKERRQ(ierr);
  if (repeat_fine_solve) {
    ierr = PetscLogStageRegister("Fine-1 solve",&stages[1]);CHKERRQ(ierr);
    ierr = PetscLogStageRegister("Fine-only solve",&stages[2]);CHKERRQ(ierr);
  }

  ierr = PetscLogStagePush(stages[0]);CHKERRQ(ierr);

  ierr = THICreate(comm,&thi);CHKERRQ(ierr);
  ierr = DMMGCreate(PETSC_COMM_WORLD,thi->nlevels,thi,&dmmg);CHKERRQ(ierr);
  {
    DA da;
    PetscInt M = 3,N = 3,P = 2;
    ierr = PetscOptionsBegin(comm,NULL,"Grid resolution options","");CHKERRQ(ierr);
    {
      ierr = PetscOptionsInt("-blatter_M","Number of elements in x-direction on coarse level","",M,&M,NULL);CHKERRQ(ierr);
      N = M;
      ierr = PetscOptionsInt("-blatter_N","Number of elements in y-direction on coarse level (if different from M)","",N,&N,NULL);CHKERRQ(ierr);
      if (thi->coarse2d) {
        ierr = PetscOptionsInt("-blatter_zlevels","Number of elements in z-direction on fine level","",thi->zlevels,&thi->zlevels,NULL);CHKERRQ(ierr);
      } else {
        ierr = PetscOptionsInt("-blatter_P","Number of elements in z-direction on coarse level","",P,&P,NULL);CHKERRQ(ierr);
      }
    }
    ierr = PetscOptionsEnd();CHKERRQ(ierr);
    if (thi->coarse2d) {
      ierr = DACreate2d(comm,DA_XYPERIODIC,DA_STENCIL_BOX,N,M,PETSC_DETERMINE,PETSC_DETERMINE,sizeof(Node)/sizeof(PetscScalar),1,0,0,&da);CHKERRQ(ierr);
      da->ops->refinehierarchy  = DARefineHierarchy_THI;
      da->ops->getinterpolation = DAGetInterpolation_THI;
      ierr = PetscObjectCompose((PetscObject)da,"THI",(PetscObject)thi);CHKERRQ(ierr);
    } else {
      ierr = DACreate3d(comm,DA_YZPERIODIC,DA_STENCIL_BOX,P,N,M,1,PETSC_DETERMINE,PETSC_DETERMINE,sizeof(Node)/sizeof(PetscScalar),1,0,0,0,&da);CHKERRQ(ierr);
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
  ierr = PetscOptionsSetValue("-dmmg_form_function_ghost","1");CHKERRQ(ierr); /* Spectacularly ugly API, our function evaluation provides ghost values */
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
  ierr = DMMGSolve(dmmg);CHKERRQ(ierr);

  ierr = PetscLogStagePop();CHKERRQ(ierr);
  ierr = THISolveStatistics(thi,dmmg,0,"Full");CHKERRQ(ierr);
  /* The first solve is complete */

  if (repeat_fine_solve && DMMGGetLevels(dmmg) > 1) {
    PetscInt nlevels = DMMGGetLevels(dmmg);
    DMMG dmmgc = dmmg[nlevels-2],dmmgf = dmmg[nlevels-1];
    Vec Xc = dmmgc->x,Xf = dmmgf->x;
    ierr = MatRestrict(dmmgf->R,Xf,Xc);CHKERRQ(ierr);
    ierr = VecPointwiseMult(Xc,Xc,dmmgf->Rscale);CHKERRQ(ierr);

    /* Solve on the level with one coarsening, this is a more stringent test of latency */
    ierr = PetscLogStagePush(stages[1]);CHKERRQ(ierr);
    ierr = (*dmmgc->solve)(dmmg,nlevels-2);CHKERRQ(ierr);
    ierr = PetscLogStagePop();CHKERRQ(ierr);
    ierr = THISolveStatistics(thi,dmmg,1,"Fine-1");CHKERRQ(ierr);

    ierr = MatInterpolate(dmmgf->R,Xc,Xf);CHKERRQ(ierr);

    /* Solve again on the finest level, this is representative of what is needed in a time-stepping code */
    ierr = PetscLogStagePush(stages[2]);CHKERRQ(ierr);
    ierr = (*dmmgf->solve)(dmmg,nlevels-1);CHKERRQ(ierr);
    ierr = PetscLogStagePop();CHKERRQ(ierr);
    ierr = THISolveStatistics(thi,dmmg,0,"Fine");CHKERRQ(ierr);
  }

  {
    PetscTruth flg;
    char filename[PETSC_MAX_PATH_LEN] = "";
    ierr = PetscOptionsGetString(PETSC_NULL,"-blatter_o",filename,sizeof(filename),&flg);CHKERRQ(ierr);
    if (flg) {
      ierr = THIDAVecView_VTK_XML(thi,DMMGGetDA(dmmg),DMMGGetx(dmmg),filename);CHKERRQ(ierr);
    }
  }

  ierr = DMMGDestroy(dmmg);CHKERRQ(ierr);
  ierr = THIDestroy(thi);CHKERRQ(ierr);

  // FIXME:  for now we quit with ISMIP-HOM output
  ierr = verbPrintf(1,grid.com,
    "\nJed's THI has now run.  ENDING ...\n\n"); CHKERRQ(ierr);
  PISMEnd();

  // FIXME:  are u,v set?
  
  ierr = compute_vertical_velocity(&u, &v, basal_melt_rate, w); CHKERRQ(ierr); 

  // vertically-average to get uvbar
  return 0;
}


PetscErrorCode BlatterStressBalance::get_max_2d_velocity(PetscReal &maxu, PetscReal &maxv) {
  // FIXME:  evaluate for efficiency 
  PetscErrorCode ierr;
  ierr = uvbar.begin_access(); CHKERRQ(ierr);
  PetscReal my_umax = 0, my_vmax = 0;
  for (PetscInt   i = grid.xs; i < grid.xs+grid.xm; ++i) {
    for (PetscInt j = grid.ys; j < grid.ys+grid.ym; ++j) {
      PISMVector2 bar = uvbar(i,j);
      my_umax = PetscMax(my_umax, PetscAbs(bar.u));
      my_vmax = PetscMax(my_vmax, PetscAbs(bar.v));
    }
  }
  ierr = uvbar.end_access(); CHKERRQ(ierr);
  ierr = PetscGlobalMax(&my_umax, &maxu, grid.com); CHKERRQ(ierr);
  ierr = PetscGlobalMax(&my_vmax, &maxv, grid.com); CHKERRQ(ierr);
  return 0;
}

PetscErrorCode BlatterStressBalance::get_3d_velocity(IceModelVec3* &u_out, IceModelVec3* &v_out, IceModelVec3* &w_out) {
  // FIXME:  this is probably not right; we need to interpolate from FEM nodal 
  //   values to get IceModel grid values
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

