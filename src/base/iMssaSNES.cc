// Copyright (C) 2004-2008 Jed Brown and Ed Bueler
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

#include <cmath>
#include <petscsnes.h>
#include "iceModel.hh"

/* follows
IceModel::computeYieldStressFromBasalShear() in src/base/iMinverseSNES.cc
and
branches/snes-macayeal/src/iMsnesmacayeal.cc?rev=31
*/

PetscErrorCode IceModel::mapUVbarSSAToSSASNESVec(DA ssasnesda, Vec &ssasnesX) {
  PetscErrorCode ierr;
  PetscScalar **u, **v;
  SSASNESNode **x;

  ierr = vubarSSA.get_array(u); CHKERRQ(ierr);
  ierr = vvbarSSA.get_array(v); CHKERRQ(ierr);
  ierr = DAVecGetArray(ssasnesda, ssasnesX, &x); CHKERRQ(ierr);

  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; i++) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; j++) {
      x[i][j].u = u[i][j];
      x[i][j].v = v[i][j];
    }
  }

  ierr = DAVecRestoreArray(ssasnesda, ssasnesX, &x); CHKERRQ(ierr);
  ierr = vubarSSA.end_access(); CHKERRQ(ierr);
  ierr = vvbarSSA.end_access(); CHKERRQ(ierr);

  return 0;
}


PetscErrorCode IceModel::mapSSASNESVecToUVbarSSA(DA ssasnesda, Vec ssasnesX) {
  PetscErrorCode ierr;
  PetscScalar **u, **v;
  SSASNESNode **x;

  ierr = vubarSSA.get_array(u); CHKERRQ(ierr);
  ierr = vvbarSSA.get_array(v); CHKERRQ(ierr);
  ierr = DAVecGetArray(ssasnesda, ssasnesX, &x); CHKERRQ(ierr);

  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; i++) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; j++) {
      u[i][j] = x[i][j].u;
      v[i][j] = x[i][j].v;
    }
  }

  ierr = DAVecRestoreArray(ssasnesda, ssasnesX, &x); CHKERRQ(ierr);
  ierr = vubarSSA.end_access(); CHKERRQ(ierr);
  ierr = vvbarSSA.end_access(); CHKERRQ(ierr);

  return 0;
}


PetscErrorCode IceModel::setbdryvalSSA(DA ssasnesda, Vec &ssasnesBV) {
  PetscErrorCode ierr;
  PetscScalar **uvbar[2];
  SSASNESNode **xBV;

  ierr = vuvbar[0].get_array(uvbar[0]); CHKERRQ(ierr);
  ierr = vuvbar[1].get_array(uvbar[1]); CHKERRQ(ierr);
  ierr = DAVecGetArray(ssasnesda, ssasnesBV, &xBV); CHKERRQ(ierr);

  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; i++) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; j++) {
      // this value will only be used where mask[i][j] == MASK_SHEET
      xBV[i][j].u = 0.5*(uvbar[0][i-1][j] + uvbar[0][i][j]);
      xBV[i][j].v = 0.5*(uvbar[1][i][j-1] + uvbar[1][i][j]);
    }
  }

  ierr = DAVecRestoreArray(ssasnesda, ssasnesBV, &xBV); CHKERRQ(ierr);
  ierr = vuvbar[0].end_access(); CHKERRQ(ierr);
  ierr = vuvbar[1].end_access(); CHKERRQ(ierr);

  return 0;
}


extern PetscScalar nu_eff(PetscScalar schoofReg, PetscScalar barB, 
                          PetscScalar u_x, PetscScalar u_y, PetscScalar v_x, PetscScalar v_y);

extern PetscErrorCode SSASNESFormFunctionLocal(DALocalInfo *info, SSASNESNode **x,
                                               SSASNESNode **f, SSASNESCtx *ctx);


PetscErrorCode IceModel::velocitySSA_SNES(IceModelVec2 vNuH[2], PetscInt *its) {
  PetscErrorCode ierr;

  SNES        snes;                 /* nonlinear solver */
  SSASNESCtx  user;
  Vec         X, R;

  ierr = SNESCreate(grid.com,&snes);CHKERRQ(ierr);

  // DA with dof = 2; otherwise same as grid.da2
  ierr = DACreate2d(grid.com, DA_XYPERIODIC, DA_STENCIL_BOX,
  	            grid.My, grid.Mx, PETSC_DECIDE, PETSC_DECIDE, 2, 1,
  	            PETSC_NULL, PETSC_NULL, &user.ssada); CHKERRQ(ierr);
  user.grid = &grid;
  user.model = this;
  
  // space for solution and residual
  ierr = DACreateGlobalVector(user.ssada, &X); CHKERRQ(ierr);
  ierr = VecDuplicate(X, &R); CHKERRQ(ierr);
  
  // fill fields in application context from current geometry
  user.ctxH = vH;
  user.ctxMask = vMask;
  user.ctxtauc = vtauc;
  user.ctxNuH[0] = vNuH[0];
  user.ctxNuH[1] = vNuH[1];

  user.ctxtaudx = vWork2d[0];
  user.ctxtaudy = vWork2d[1];
  ierr = computeDrivingStress(user.ctxtaudx,user.ctxtaudy); CHKERRQ(ierr);

  ierr = VecDuplicate(X, &user.ctxBV); CHKERRQ(ierr);

  // fill in parameters in app ctx
  user.schoofReg = PetscSqr(regularizingVelocitySchoof/regularizingLengthSchoof);
  user.constantHardness = constantHardnessForSSA;
  
  // flags
  user.leaveNuHAlone = leaveNuHAloneSSA;
  user.useConstantHardness = useConstantHardnessForSSA;
  
  // so that first call can return reasonable first guess
  user.callcount = 0;
  
  Mat J;
  ierr = SNESSetFunction(snes,R,SNESDAFormFunction,&user); CHKERRQ(ierr);
  ierr = MatCreateSNESMF(snes,&J); CHKERRQ(ierr);
  ierr = SNESSetJacobian(snes,J,J,SNESDAComputeJacobian,&user);CHKERRQ(ierr);

  ierr = DASetLocalFunction(user.ssada,(DALocalFunction1)SSASNESFormFunctionLocal);
            CHKERRQ(ierr);

  // no Jacobian matrix at all
  ierr = PetscOptionsSetValue("-snes_mf", PETSC_NULL); CHKERRQ(ierr);
 
  ierr = mapUVbarSSAToSSASNESVec(user.ssada, X); CHKERRQ(ierr);
  
  ierr = setbdryvalSSA(user.ssada, user.ctxBV); CHKERRQ(ierr);
  
  PetscViewer viewer;
  ierr = PetscViewerASCIIGetStdout(PETSC_COMM_WORLD,&viewer); CHKERRQ(ierr);
  ierr = PetscViewerSetFormat(viewer,PETSC_VIEWER_ASCII_INDEX); CHKERRQ(ierr);
  ierr = VecView(user.ctxBV,viewer); CHKERRQ(ierr);
  //PetscViewerDestroy(viewer);

  ierr = SNESSolve(snes,PETSC_NULL,X); CHKERRQ(ierr);

  // some feedback appropriate, but FIXME
  PetscReal           resnorm;
  SNESConvergedReason reason;
  ierr = SNESGetIterationNumber(snes,its);     CHKERRQ(ierr); 
  ierr = SNESGetConvergedReason(snes,&reason); CHKERRQ(ierr);
  ierr = VecNorm(R,NORM_INFINITY,&resnorm);    CHKERRQ(ierr);
  ierr = PetscPrintf(grid.com,
     "    done ... # of Newton iters = %d;  conv reason = %d;  |residual|_infty = %9.3e\n",
     *its, (int)reason, resnorm); CHKERRQ(ierr);

  ierr = mapSSASNESVecToUVbarSSA(user.ssada, X); CHKERRQ(ierr);

  // de-allocate
  ierr = VecDestroy(user.ctxBV);CHKERRQ(ierr);      
  ierr = VecDestroy(X);CHKERRQ(ierr);
  ierr = VecDestroy(R);CHKERRQ(ierr);      
  ierr = SNESDestroy(snes);CHKERRQ(ierr);

  return 0;
}


#define sqr PetscSqr
PetscScalar nu_eff(PetscScalar schoofReg, PetscScalar barB, 
                   PetscScalar u_x, PetscScalar u_y, PetscScalar v_x, PetscScalar v_y) {
  // constant \bar B case; for Test I and EISMINT-Ross and MISMIP
  // FIXME: assumes n=3
  return barB * 0.5 *
             pow(schoofReg + PetscSqr(u_x) + PetscSqr(v_y) + 0.25*PetscSqr(u_y+v_x) + u_x*v_y,
                 -(1.0/3.0));
}


PetscErrorCode SSASNESFormFunctionLocal(DALocalInfo *info, SSASNESNode **x,
                                        SSASNESNode **f, SSASNESCtx *ctx) {
  PetscErrorCode ierr;

//  ierr = PetscPrintf(ctx->grid->com, 
//                     "entering SSASNESFormFunctionLocal() with callcount=%d\n",
//                     ctx->callcount); CHKERRQ(ierr);

  PetscScalar **H, **mask, **tauc, **taudx, **taudy, **nuH[2];
  SSASNESNode **xBV;
  PetscInt Mx, My, xs, ys, xm, ym;

  // use transpose as in IceModel
  Mx = info->my; My = info->mx;
  xs = info->ys; ys = info->xs;
  xm = info->ym; ym = info->xm;

  const PetscScalar dx = ctx->grid->dx,
                    dy = ctx->grid->dy,
                    sc = dx * dy;

  ierr =      ctx->ctxH.get_array(H);      CHKERRQ(ierr);
  ierr =   ctx->ctxMask.get_array(mask);   CHKERRQ(ierr);
  ierr =   ctx->ctxtauc.get_array(tauc);   CHKERRQ(ierr);
  ierr =  ctx->ctxtaudx.get_array(taudx);  CHKERRQ(ierr);
  ierr =  ctx->ctxtaudy.get_array(taudy);  CHKERRQ(ierr);
  ierr = ctx->ctxNuH[0].get_array(nuH[0]); CHKERRQ(ierr);
  ierr = ctx->ctxNuH[1].get_array(nuH[1]); CHKERRQ(ierr);

  ierr = DAVecGetArray(ctx->ssada, ctx->ctxBV, &xBV); CHKERRQ(ierr);

  for (PetscInt i=xs; i<xs+xm; i++) {
    for (PetscInt j=ys; j<ys+ym; j++) {

      int maskval = static_cast<int>(floor(mask[i][j] + 0.5)); // c.f. IceModel::intMask()
      const int MASK_SHEET    = 1;
      const int MASK_DRAGGING = 2;

      if (maskval == MASK_SHEET) {
        // here we assign boundary values
        f[i][j].u = x[i][j].u - xBV[i][j].u;
        f[i][j].v = x[i][j].v - xBV[i][j].v;
      } else {
        // main case: f.d. approx of PDE gives sys of nonlinear eqns
        const PetscScalar 
          u_x_im = (x[i][j].u - x[i-1][j].u) / dx,
          u_x_ip = (x[i+1][j].u - x[i][j].u) / dx,
          u_x_jm = (+ (x[i+1][j-1].u + x[i+1][j].u)
                    - (x[i-1][j-1].u + x[i-1][j].u)) / (4 * dx),
          u_x_jp = (+ (x[i+1][j].u + x[i+1][j+1].u)
                    - (x[i-1][j].u + x[i-1][j+1].u)) / (4 * dx),
  
          u_y_im = (+ (x[i-1][j+1].u + x[i][j+1].u)
                    - (x[i-1][j-1].u + x[i][j-1].u)) / (4 * dy),
          u_y_ip = (+ (x[i][j+1].u + x[i+1][j+1].u)
                    - (x[i][j-1].u + x[i+1][j-1].u)) / (4 * dy),
          u_y_jm = (x[i][j].u - x[i][j-1].u) / dy,
          u_y_jp = (x[i][j+1].u - x[i][j].u) / dy,

          v_x_im = (x[i][j].v - x[i-1][j].v) / dx,
          v_x_ip = (x[i+1][j].v - x[i][j].v) / dx,
          v_x_jm = (+ (x[i+1][j-1].v + x[i+1][j].v)
                    - (x[i-1][j-1].v + x[i-1][j].v)) / (4 * dx),
          v_x_jp = (+ (x[i+1][j].v + x[i+1][j+1].v)
                    - (x[i-1][j].v + x[i-1][j+1].v)) / (4 * dx),
  
          v_y_im = (+ (x[i-1][j+1].v + x[i][j+1].v)
                    - (x[i-1][j-1].v + x[i][j-1].v)) / (4 * dy),
          v_y_ip = (+ (x[i][j+1].v + x[i+1][j+1].v)
                    - (x[i][j-1].v + x[i+1][j-1].v)) / (4 * dy),
          v_y_jm = (x[i][j].v - x[i][j-1].v) / dy,
          v_y_jp = (x[i][j+1].v - x[i][j].v) / dy;
  
        const PetscScalar H_im = 0.5 * (H[i-1][j] + H[i][j]),
                          H_ip = 0.5 * (H[i][j] + H[i+1][j]),
                          H_jm = 0.5 * (H[i][j-1] + H[i][j]),
                          H_jp = 0.5 * (H[i][j] + H[i][j+1]);

        PetscScalar nu_im, nu_ip, nu_jm, nu_jp;
      
        if (ctx->leaveNuHAlone == PETSC_TRUE) {
          // unprotected division o.k. for Test J but probably not otherwise
          nu_im = nuH[0][i-1][j] / H_im;
          nu_ip = nuH[0][i][j]   / H_ip;
          nu_jm = nuH[1][i][j-1] / H_jm;
          nu_jp = nuH[1][i][j]   / H_jp;
        } else {
          if (ctx->useConstantHardness == PETSC_FALSE) {
            SETERRQ(1,"need ctx->useConstantHardness == PETSC_TRUE");
          }
          // FIXME:  assumes constant hardness (constant Bbar); fix for thermocoupled
//          const PetscScalar reg = (ctx->callcount == 0) ? 100.0 * ctx->schoofReg : ctx->schoofReg;
/*
          const PetscScalar reg = ctx->schoofReg;
          nu_im = nu_eff(reg, ctx->constantHardness, u_x_im, u_y_im, v_x_im, v_y_im);
          nu_ip = nu_eff(reg, ctx->constantHardness, u_x_ip, u_y_ip, v_x_ip, v_y_ip);
          nu_jm = nu_eff(reg, ctx->constantHardness, u_x_jm, u_y_jm, v_x_jm, v_y_jm);
          nu_jp = nu_eff(reg, ctx->constantHardness, u_x_jp, u_y_jp, v_x_jp, v_y_jp);
*/
          const PetscScalar constnu = 30.0 * 1e6 * secpera; // = 9.45e14
          nu_im = constnu;
          nu_ip = constnu;
          nu_jm = constnu;
          nu_jp = constnu;
        }
      
        f[i][j].u = sc *( + (2 / dx) * (+ (nu_ip * H_ip * (2 * u_x_ip + v_y_ip))
                                        - (nu_im * H_im * (2 * u_x_im + v_y_im)))
                          + (1 / dy) * (+ (nu_jp * H_jp * (u_y_jp + v_x_jp))
                                        - (nu_jm * H_jm * (u_y_jm + v_x_jm))) )
                     + sc * taudx[i][j];
        f[i][j].v = sc *( + (2 / dy) * (+ (nu_jp * H_jp * (2 * v_y_jp + u_x_jp))
                                        - (nu_jm * H_jm * (2 * v_y_jm + u_x_jm)))
                          + (1 / dx) * (+ (nu_ip * H_ip * (u_y_ip + v_x_ip))
                                        - (nu_im * H_im * (u_y_im + v_x_im))) )
                     + sc * taudy[i][j];
  
        if (maskval == MASK_DRAGGING) {
          // this bypasses IceModel::basalDragx(), IceModel::basalDragy()
//         const PetscScalar smallVel = (ctx->callcount == 0) ? 10.0 / secpera : 0.01 / secpera;
          const PetscScalar smallVel = 0.01 / secpera;
          const PetscScalar denom = sqrt(PetscSqr(smallVel) 
                                         + PetscSqr(x[i][j].u) + PetscSqr(x[i][j].v));
          f[i][j].u -= sc * tauc[i][j] * x[i][j].u / denom;
          f[i][j].v -= sc * tauc[i][j] * x[i][j].v / denom;
          //f[i][j].u -= sc * ctx->model->basal->drag(tauc[i][j], x[i][j].u, x[i][j].v) * x[i][j].u;
          //f[i][j].v -= sc * ctx->model->basal->drag(tauc[i][j], x[i][j].u, x[i][j].v) * x[i][j].v;
        }

      }
#if 0
      if (isnan(f[i][j].u) || isnan(f[i][j].v)
          || (i == 1 && j == 30 && false)) {
        ierr = PetscPrintf(ctx->model->grid.com, "i, j, F.u, F.v = %4d %4d %8e %8e\n",
                           i, j, f[i][j].u, f[i][j].v); CHKERRQ(ierr);
        ierr = PetscPrintf(ctx->model->grid.com, "nu = %8e %8e %8e %8e\n",
                           nu_im, nu_ip, nu_jm, nu_jp); CHKERRQ(ierr);
        ierr = PetscPrintf(ctx->model->grid.com, "u_x = %8e %8e %8e %8e\n",
                           u_x_im, u_x_ip, u_x_jm, u_x_jp); CHKERRQ(ierr);
        ierr = PetscPrintf(ctx->model->grid.com, "u_y = %8e %8e %8e %8e\n",
                           u_y_im, u_y_ip, u_y_jm, u_y_jp); CHKERRQ(ierr);
        ierr = PetscPrintf(ctx->model->grid.com, "v_x = %8e %8e %8e %8e\n",
                           v_x_im, v_x_ip, v_x_jm, v_x_jp); CHKERRQ(ierr);
        ierr = PetscPrintf(ctx->model->grid.com, "v_y = %8e %8e %8e %8e\n",
                           v_y_im, v_y_ip, v_y_jm, v_y_jp); CHKERRQ(ierr);
      }
#endif

    }
  }

  ierr = DAVecRestoreArray(ctx->ssada, ctx->ctxBV, &xBV); CHKERRQ(ierr);
  
  ierr =      ctx->ctxH.end_access(); CHKERRQ(ierr);
  ierr =   ctx->ctxMask.end_access(); CHKERRQ(ierr);
  ierr =   ctx->ctxtauc.end_access(); CHKERRQ(ierr);
  ierr =  ctx->ctxtaudx.end_access(); CHKERRQ(ierr);
  ierr =  ctx->ctxtaudy.end_access(); CHKERRQ(ierr);
  ierr = ctx->ctxNuH[0].end_access(); CHKERRQ(ierr);
  ierr = ctx->ctxNuH[1].end_access(); CHKERRQ(ierr);

  ctx->callcount++;
  return 0;
}

