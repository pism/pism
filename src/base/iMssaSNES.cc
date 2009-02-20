// Copyright (C) 2004-2009 Jed Brown and Ed Bueler
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

  // Communicate so that we have stencil width for evaluation of effective viscosity (and geometry)
  ierr = vubarSSA.beginGhostComm(); CHKERRQ(ierr);
  ierr = vubarSSA.endGhostComm(); CHKERRQ(ierr);
  ierr = vvbarSSA.beginGhostComm(); CHKERRQ(ierr);
  ierr = vvbarSSA.endGhostComm(); CHKERRQ(ierr);
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


PetscErrorCode IceModel::solvefeedback(SNES snes, Vec residual) {
  PetscErrorCode      ierr;
  PetscInt            its;
  PetscReal           resnorm;
  SNESConvergedReason reason;
  ierr = SNESGetIterationNumber(snes,&its);     CHKERRQ(ierr); 
  ierr = SNESGetConvergedReason(snes,&reason); CHKERRQ(ierr);
  ierr = VecNorm(residual,NORM_INFINITY,&resnorm);    CHKERRQ(ierr);
  ierr = PetscPrintf(grid.com,
     "    snes done ... # of Newton iters = %d;  conv reason = %d;  |residual|_infty = %9.3e\n",
     its, (int)reason, resnorm); CHKERRQ(ierr);
  return 0;
}


PetscErrorCode IceModel::getNuFromNuH(IceModelVec2 vNuH[2], SSASNESCtx *user) {
  PetscErrorCode ierr;
  PetscScalar **H, **nuH[2], **nu[2];

  // communicate to get ghosts for H
  ierr = vH.beginGhostComm(); CHKERRQ(ierr);
  ierr = vH.endGhostComm(); CHKERRQ(ierr);

  ierr =             vH.get_array(H);      CHKERRQ(ierr);
  ierr =        vNuH[0].get_array(nuH[0]); CHKERRQ(ierr);
  ierr =        vNuH[1].get_array(nuH[1]); CHKERRQ(ierr);
  ierr = user->ctxNu[0].get_array(nu[0]);  CHKERRQ(ierr);
  ierr = user->ctxNu[1].get_array(nu[1]);  CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; i++) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; j++) {
      const PetscScalar H_ip = 0.5 * (H[i][j] + H[i+1][j]),
                        H_jp = 0.5 * (H[i][j] + H[i][j+1]);      
      // unprotected division o.k. for Test J but probably not otherwise
      nu[0][i][j] = nuH[0][i][j] / H_ip;
      nu[1][i][j] = nuH[1][i][j] / H_jp;
    }
  }
  ierr =             vH.end_access(); CHKERRQ(ierr);
  ierr =        vNuH[0].end_access(); CHKERRQ(ierr);
  ierr =        vNuH[1].end_access(); CHKERRQ(ierr);
  ierr = user->ctxNu[0].end_access(); CHKERRQ(ierr);
  ierr = user->ctxNu[1].end_access(); CHKERRQ(ierr);

  // communicate to get ghosts for ctxNu[2]
  ierr = user->ctxNu[0].beginGhostComm(); CHKERRQ(ierr);
  ierr = user->ctxNu[1].beginGhostComm(); CHKERRQ(ierr);
  ierr = user->ctxNu[0].endGhostComm(); CHKERRQ(ierr);
  ierr = user->ctxNu[1].endGhostComm(); CHKERRQ(ierr);  
  return 0;
}


extern PetscScalar nu_eff(PetscScalar schoofReg, PetscScalar barB, 
                          PetscScalar u_x, PetscScalar u_y, PetscScalar v_x, PetscScalar v_y);
extern PetscErrorCode basalstress(PetscTruth useIceModelBasal, PlasticBasalType *basal,
                                  PetscScalar u, PetscScalar v, PetscScalar tauc,
                                  PetscScalar &taubx, PetscScalar &tauby);
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

  // space for solution and residual
  ierr = DACreateGlobalVector(user.ssada, &X); CHKERRQ(ierr);
  ierr = VecDuplicate(X, &R); CHKERRQ(ierr);
  
  // build application context:
  //   get needed members of IceModel
  user.grid = &grid;
  user.basal = basal;
  //   fill fields in application context from current geometry
  user.ctxH = vH;
  user.ctxMask = vMask;
  user.ctxtauc = vtauc;
  //   get driving stress
  user.ctxtaudx = vWork2d[2];  // note that vNuH[2] might be {vWork2d[0],vWork2d[1]}
  user.ctxtaudy = vWork2d[3];
  ierr = computeDrivingStress(user.ctxtaudx,user.ctxtaudy); CHKERRQ(ierr);
  //   space for effective viscosity under iteration
  ierr = user.ctxNu[0].create(grid, "nu[0]_SSASNES", true); CHKERRQ(ierr);
  ierr = user.ctxNu[1].create(grid, "nu[1]_SSASNES", true); CHKERRQ(ierr);
  //   fill in boundary values
  ierr = VecDuplicate(X, &user.ctxBV); CHKERRQ(ierr);
  ierr = setbdryvalSSA(user.ssada, user.ctxBV); CHKERRQ(ierr);

  //   fill in parameters and flags in app ctx
  user.schoofReg = PetscSqr(regularizingVelocitySchoof/regularizingLengthSchoof);
  user.constantHardness = constantHardnessForSSA;
  user.useConstantHardness = useConstantHardnessForSSA;
  user.useConstantNu = useConstantNuHForSSA;  
  user.usePlasticBasalType = PETSC_TRUE;
  
  // set up SNES function
  ierr = SNESSetFunction(snes,R,SNESDAFormFunction,&user); CHKERRQ(ierr);
  ierr = DASetLocalFunction(user.ssada,(DALocalFunction1)SSASNESFormFunctionLocal);
            CHKERRQ(ierr);

  // last stage: ask for user options to override these settings
  ierr = SNESSetFromOptions(snes);CHKERRQ(ierr);

  ierr = mapUVbarSSAToSSASNESVec(user.ssada, X); CHKERRQ(ierr);

#if 0
  PetscViewer viewer;
  ierr = PetscViewerASCIIGetStdout(PETSC_COMM_WORLD,&viewer); CHKERRQ(ierr);
  ierr = PetscViewerSetFormat(viewer,PETSC_VIEWER_ASCII_INDEX); CHKERRQ(ierr);
  ierr = VecView(user.ctxBV,viewer); CHKERRQ(ierr);
#endif

  // possible sequence:
  //   picard loop:
  //      compute user.ctxNu[2]
  //      set flag so solve holds user.ctxNu[2] fixed
  //      solve
  //      check weak (1e-2) nu convergence criterion
  //   Newton-Krylov:
  //      set flag so solve recomputes viscosity as needed

#if 0
    // linear till
    basal->pseudo_plastic = PETSC_TRUE;
    basal->pseudo_q = 1.0;
    ierr = user.ctxtauc.set(5.703978e+03); CHKERRQ(ierr);
    //ierr = user.ctxtauc.set(0.0); CHKERRQ(ierr);
    //ierr = user.ctxtauc.scale(1.0e-3); CHKERRQ(ierr);
#endif

  // compute Nu from NuH if that is supposed to be the source
  if (leaveNuHAloneSSA == PETSC_TRUE) {
    ierr = updateNuViewers(vNuH, vNuH, true); CHKERRQ(ierr);
    ierr = getNuFromNuH(vNuH, &user); CHKERRQ(ierr);
    user.useStoredNu = PETSC_TRUE;
  } else {
    // idea here is to use three steps of Picard (just to try something)

#if 1
    KSP ksp = SSAKSP;
    Mat A = SSAStiffnessMatrix;
    Vec x = SSAX, rhs = SSARHS; // solve  A x = rhs
    PetscInt    kspits;
    KSPConvergedReason  kspreason;
    ierr = assembleSSARhs((computeSurfGradInwardSSA == PETSC_TRUE), rhs); CHKERRQ(ierr);
    ierr = computeEffectiveViscosity(vNuH,ssaEpsilon); CHKERRQ(ierr);
    ierr = updateNuViewers(vNuH, vNuH, true); CHKERRQ(ierr);
    ierr = assembleSSAMatrix(true, vNuH, A); CHKERRQ(ierr);
    ierr = verbPrintf(3,grid.com, "A:"); CHKERRQ(ierr);
    // call PETSc to solve linear system by iterative method
    ierr = KSPSetOperators(ksp, A, A, DIFFERENT_NONZERO_PATTERN); CHKERRQ(ierr);
    ierr = KSPSetFromOptions(ksp); CHKERRQ(ierr);
    ierr = KSPSolve(ksp, rhs, x); CHKERRQ(ierr); // SOLVE
    ierr = KSPGetIterationNumber(ksp, &kspits); CHKERRQ(ierr);
    ierr = KSPGetConvergedReason(ksp, &kspreason); CHKERRQ(ierr);
    ierr = verbPrintf(3,grid.com, "S:%d,%d: ", kspits, kspreason); CHKERRQ(ierr);
    // finish iteration and report to standard out
    ierr = moveVelocityToDAVectors(x); CHKERRQ(ierr);
#endif

    ierr = computeEffectiveViscosity(vNuH,ssaEpsilon); CHKERRQ(ierr);
    ierr = updateNuViewers(vNuH, vNuH, true); CHKERRQ(ierr);
    ierr = getNuFromNuH(vNuH, &user); CHKERRQ(ierr);

    user.useStoredNu = PETSC_TRUE;
    ierr = SNESSolve(snes,PETSC_NULL,X); CHKERRQ(ierr);
    ierr = solvefeedback(snes, R); CHKERRQ(ierr);

    ierr = mapSSASNESVecToUVbarSSA(user.ssada, X); CHKERRQ(ierr);
    ierr = computeEffectiveViscosity(vNuH,ssaEpsilon); CHKERRQ(ierr);
    ierr = updateNuViewers(vNuH, vNuH, true); CHKERRQ(ierr);
    ierr = getNuFromNuH(vNuH, &user); CHKERRQ(ierr);

    user.useStoredNu = PETSC_TRUE;
    ierr = SNESSolve(snes,PETSC_NULL,X); CHKERRQ(ierr);
    ierr = solvefeedback(snes, R); CHKERRQ(ierr);

    ierr = mapSSASNESVecToUVbarSSA(user.ssada, X); CHKERRQ(ierr);
    ierr = computeEffectiveViscosity(vNuH,ssaEpsilon); CHKERRQ(ierr);
    ierr = updateNuViewers(vNuH, vNuH, true); CHKERRQ(ierr);
    ierr = getNuFromNuH(vNuH, &user); CHKERRQ(ierr);

    user.useStoredNu = PETSC_TRUE;
    ierr = SNESSolve(snes,PETSC_NULL,X); CHKERRQ(ierr);
    ierr = solvefeedback(snes, R); CHKERRQ(ierr);

    user.useStoredNu = PETSC_FALSE;
  }
  
  ierr = SNESSolve(snes,PETSC_NULL,X); CHKERRQ(ierr);
  ierr = solvefeedback(snes, R); CHKERRQ(ierr);
  ierr = SNESGetIterationNumber(snes,its); CHKERRQ(ierr); 

  ierr = mapSSASNESVecToUVbarSSA(user.ssada, X); CHKERRQ(ierr);

  // de-allocate
  ierr = VecDestroy(user.ctxBV);CHKERRQ(ierr);      
  ierr = VecDestroy(X);CHKERRQ(ierr);
  ierr = VecDestroy(R);CHKERRQ(ierr);      
  ierr = SNESDestroy(snes);CHKERRQ(ierr);
  ierr = user.ctxNu[0].destroy(); CHKERRQ(ierr);
  ierr = user.ctxNu[1].destroy(); CHKERRQ(ierr);

  return 0;
}


#define sqr PetscSqr
static PetscScalar nu_eff(PetscTruth useConstantNu, PetscScalar schoofReg, PetscScalar barB, 
                   PetscScalar u_x, PetscScalar u_y, PetscScalar v_x, PetscScalar v_y) {
  //PetscPrintf(PETSC_COMM_WORLD, "nu_eff() turned off\n"); PetscEnd();
  if (useConstantNu == PETSC_TRUE) {
    //return 30.0 * 1e6 * secpera; // 30 MPa a^-1 = 9.45e14 Pa s^-1
    const PetscScalar strainrate = (100.0 / secpera) / 100.0e3;  // = 1e-3 s^-1; typical strain rate
    return barB * 0.5 * pow(PetscSqr(strainrate), -(1.0/3.0));
  } else {
    // constant \bar B case; for Test I and EISMINT-Ross and MISMIP
    // FIXME: assumes n=3
    return barB * 0.5 *
             pow(schoofReg + PetscSqr(u_x) + PetscSqr(v_y) + 0.25*PetscSqr(u_y+v_x) + u_x*v_y,
                 -(1.0/3.0));
  }
}


PetscErrorCode basalstress(PetscTruth usePlasticBasalType, PlasticBasalType *basal,
                           PetscScalar u, PetscScalar v, PetscScalar tauc,
                           PetscScalar &taubx, PetscScalar &tauby) {
//  PetscPrintf(PETSC_COMM_WORLD, "basalstress() turned off\n"); PetscEnd();
  if (usePlasticBasalType == PETSC_TRUE) {
    taubx = - basal->drag(tauc, u, v) * u;
    tauby = - basal->drag(tauc, u, v) * v;
  } else {
    // value stated in Hulbe&MacAyeal1999 for ice stream E
    const PetscScalar beta = 1.8e9;  // Pa s m^{-1}
    taubx = - beta * u;
    tauby = - beta * v;
  }
  return 0;
}


PetscErrorCode SSASNESFormFunctionLocal(DALocalInfo *info, SSASNESNode **x,
                                        SSASNESNode **f, SSASNESCtx *ctx) {
  PetscErrorCode ierr;

  PetscScalar **H, **mask, **tauc, **taudx, **taudy, **nu[2];
  SSASNESNode **xBV;
  PetscInt xs, ys, xm, ym;

  // use transpose as in IceModel
  xs = info->ys; ys = info->xs;
  xm = info->ym; ym = info->xm;

  const PetscScalar dx = ctx->grid->dx,
                    dy = ctx->grid->dy,
                    sc = dx * dy;

  ierr =     ctx->ctxH.get_array(H);     CHKERRQ(ierr);
  ierr =  ctx->ctxMask.get_array(mask);  CHKERRQ(ierr);
  ierr =  ctx->ctxtauc.get_array(tauc);  CHKERRQ(ierr);
  ierr = ctx->ctxtaudx.get_array(taudx); CHKERRQ(ierr);
  ierr = ctx->ctxtaudy.get_array(taudy); CHKERRQ(ierr);
  ierr = ctx->ctxNu[0].get_array(nu[0]); CHKERRQ(ierr);
  ierr = ctx->ctxNu[1].get_array(nu[1]); CHKERRQ(ierr);

  ierr = DAVecGetArray(ctx->ssada, ctx->ctxBV, &xBV); CHKERRQ(ierr);

  for (PetscInt i=xs; i<xs+xm; i++) {
    for (PetscInt j=ys; j<ys+ym; j++) {

      int maskval = PismIntMask(mask[i][j]);

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
      
        if (ctx->useStoredNu == PETSC_TRUE) {
          nu_im = nu[0][i-1][j];
          nu_ip = nu[0][i][j];
          nu_jm = nu[1][i][j-1];
          nu_jp = nu[1][i][j];
        } else {
          if (ctx->useConstantHardness == PETSC_FALSE) {
            SETERRQ(1,"need ctx->useConstantHardness == PETSC_TRUE");
          }
          nu_im = nu_eff(ctx->useConstantNu, ctx->schoofReg, ctx->constantHardness,
                         u_x_im, u_y_im, v_x_im, v_y_im);
          nu_ip = nu_eff(ctx->useConstantNu, ctx->schoofReg, ctx->constantHardness,
                         u_x_ip, u_y_ip, v_x_ip, v_y_ip);
          nu_jm = nu_eff(ctx->useConstantNu, ctx->schoofReg, ctx->constantHardness,
                         u_x_jm, u_y_jm, v_x_jm, v_y_jm);
          nu_jp = nu_eff(ctx->useConstantNu, ctx->schoofReg, ctx->constantHardness,
                         u_x_jp, u_y_jp, v_x_jp, v_y_jp);          
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
          PetscScalar mytaubx, mytauby;
          ierr = basalstress(ctx->usePlasticBasalType, ctx->basal,
                      x[i][j].u,  x[i][j].v, tauc[i][j], mytaubx, mytauby); CHKERRQ(ierr);
          f[i][j].u += sc * mytaubx;
          f[i][j].v += sc * mytauby;
//          f[i][j].u -= sc * mytaubx;
//          f[i][j].v -= sc * mytauby;
        }

      }

    }
  }

  ierr = DAVecRestoreArray(ctx->ssada, ctx->ctxBV, &xBV); CHKERRQ(ierr);
  
  ierr =     ctx->ctxH.end_access(); CHKERRQ(ierr);
  ierr =  ctx->ctxMask.end_access(); CHKERRQ(ierr);
  ierr =  ctx->ctxtauc.end_access(); CHKERRQ(ierr);
  ierr = ctx->ctxtaudx.end_access(); CHKERRQ(ierr);
  ierr = ctx->ctxtaudy.end_access(); CHKERRQ(ierr);
  ierr = ctx->ctxNu[0].end_access(); CHKERRQ(ierr);
  ierr = ctx->ctxNu[1].end_access(); CHKERRQ(ierr);

  return 0;
}

