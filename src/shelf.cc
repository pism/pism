// Copyright (C) 2004-2007 Jed Brown and Ed Bueler
//
// This file is part of Pism.
//
// Pism is free software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the Free Software
// Foundation; either version 2 of the License, or (at your option) any later
// version.
//
// Pism is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
// details.
//
// You should have received a copy of the GNU General Public License
// along with Pism; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

static char help[] =
  "Driver for testing ice stream (dragging ice shelf) model.  Implements verification tests.";

#include <cmath>
#include "iceModel.hh"
#include "iceCompModel.hh"

#define sec(x) (1.0 / cos((x)))

#define EMBAY_VARS                              \
  PetscScalar x = xx / grid.p->Lx;              \
  PetscScalar x = yy / grid.p->Ly;              \
  PetscScalar B = 135720960;


class ShelfModel : public IceModel {
public:
  ShelfModel(IceGrid &g, GlenIce &i);
  virtual PetscErrorCode initFromOptions();
  virtual PetscErrorCode run();
protected:
  virtual PetscErrorCode createVecs();
  virtual PetscErrorCode destroyVecs();
  virtual PetscErrorCode afterInitHook();
  virtual PetscScalar basalDragx(PetscScalar **u, PetscScalar **v,
                                 PetscInt i, PetscInt j) const;
  virtual PetscScalar basalDragy(PetscScalar **u, PetscScalar **v,
                                 PetscInt i, PetscInt j) const;
  virtual PetscErrorCode computeEffectiveViscosity(Vec vNu[2], PetscReal epsilon);
  PetscErrorCode initFloatingSlab();
  PetscErrorCode initEmbayment();
  PetscErrorCode initBoring();
  Vec vbetax, vbetay;
  Vec vNu[2];
  Vec vu_exact, vv_exact;
  PetscScalar **basalBetax, **basalBetay; // This is sloppy
};


ShelfModel::ShelfModel(IceGrid &g, GlenIce &i)
  : IceModel(g, i) {
  setThermalBedrock(PETSC_FALSE);
  setUseMacayealVelocity(PETSC_TRUE);
  setDoGrainSize(PETSC_FALSE);
}


PetscErrorCode ShelfModel::createVecs() {
  PetscErrorCode ierr;

  if (createVecs_done) {
    ierr = destroyVecs(); CHKERRQ(ierr);
  }

  ierr = IceModel::createVecs(); CHKERRQ(ierr);
  ierr = VecDuplicate(vh, &vbetax); CHKERRQ(ierr);
  ierr = VecDuplicate(vh, &vbetay); CHKERRQ(ierr);
  ierr = VecDuplicate(vh, &vNu[0]); CHKERRQ(ierr);
  ierr = VecDuplicate(vh, &vNu[1]); CHKERRQ(ierr);
  ierr = VecDuplicate(vh, &vu_exact); CHKERRQ(ierr);
  ierr = VecDuplicate(vh, &vv_exact); CHKERRQ(ierr);
  
  return 0;
}


PetscErrorCode ShelfModel::destroyVecs() {
  PetscErrorCode ierr;

  if (basalBetax) {
    ierr = DAVecRestoreArray(grid.da2, vbetax, &basalBetax); CHKERRQ(ierr);
  }
  if (basalBetay) {
    ierr = DAVecRestoreArray(grid.da2, vbetay, &basalBetay); CHKERRQ(ierr);
  }
  
  ierr = VecDestroy(vbetax); CHKERRQ(ierr);
  ierr = VecDestroy(vbetay); CHKERRQ(ierr);
  ierr = VecDestroy(vNu[0]); CHKERRQ(ierr);
  ierr = VecDestroy(vNu[1]); CHKERRQ(ierr);
  ierr = VecDestroy(vu_exact); CHKERRQ(ierr);
  ierr = VecDestroy(vv_exact); CHKERRQ(ierr);

  return 0;
}

PetscScalar ShelfModel::basalDragx(PetscScalar **u, PetscScalar **v,
                                 PetscInt i, PetscInt j) const {
  return basalBetax[i][j];
  return 0;
  return 1e0 * secpera;
  return IceModel::basalDragx(u, v, i, j);
}

PetscScalar ShelfModel::basalDragy(PetscScalar **u, PetscScalar **v,
                                 PetscInt i, PetscInt j) const {
  return basalBetay[i][j];
  return 0;
  return 1e0 * secpera;
  return IceModel::basalDragy(u, v, i, j);
}

PetscErrorCode ShelfModel::computeEffectiveViscosity(Vec im_vNu[2], PetscReal epsilon) {
  PetscErrorCode ierr;

#define CORRECT_VELOCITY 0

  
#if (CORRECT_VELOCITY)
  ierr = VecCopy(vubar, vWork2d[6]); CHKERRQ(ierr);
  ierr = VecCopy(vvbar, vWork2d[7]); CHKERRQ(ierr);
  // Start with correct values
  ierr = VecCopy(vu_exact, vubar); CHKERRQ(ierr);
  ierr = VecCopy(vv_exact, vvbar); CHKERRQ(ierr);
#endif

  // VecPointwiseMult(w, x, y) : w = x * y
#if (0)
  ierr = VecPointwiseMult(im_vNu[0], vNu[0], vH); CHKERRQ(ierr);
  ierr = VecPointwiseMult(im_vNu[1], vNu[1], vH); CHKERRQ(ierr);
#else
  // Start with nothing; it should not be used anyway
  ierr = VecSet(im_vNu[0], 0); CHKERRQ(ierr);
  ierr = VecSet(im_vNu[1], 0); CHKERRQ(ierr);

  ierr = IceModel::computeEffectiveViscosity(im_vNu, epsilon); CHKERRQ(ierr);

  PetscScalar **mask, **H, **nu_exact[2], **nu[2];
  ierr = DAVecGetArray(grid.da2, vMask, &mask); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vH, &H); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vNu[0], &nu_exact[0]); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vNu[1], &nu_exact[1]); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, im_vNu[0], &nu[0]); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, im_vNu[1], &nu[1]); CHKERRQ(ierr);
  for (PetscInt o=0; o < 2; o++) {
    for (PetscInt i=grid.xs; i < grid.xs+grid.xm; i++) {
      for (PetscInt j=grid.ys; j < grid.ys+grid.ym; j++) {
        const PetscInt      oi = 1-o, oj=o;
        if (intMask(mask[i][j]) == MASK_SHEET && intMask(mask[i+oi][j+oj]) == MASK_SHEET) {
          nu[o][i + oi][j + oj] = nu_exact[o][i + oi][j + oj]
            * 0.5 * (H[i][j] + H[i + oi][j + oj]);
        }
      }
    }
  }
  ierr = DAVecRestoreArray(grid.da2, vMask, &mask); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vH, &H); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vNu[0], &nu_exact[0]); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vNu[1], &nu_exact[1]); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, im_vNu[0], &nu[0]); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, im_vNu[1], &nu[1]); CHKERRQ(ierr);

  ierr = IceModel::computeEffectiveViscosity(im_vNu, epsilon); CHKERRQ(ierr);
#endif
  

#if (CORRECT_VELOCITY)
  ierr = VecCopy(vWork2d[6], vubar); CHKERRQ(ierr);
  ierr = VecCopy(vWork2d[7], vvbar); CHKERRQ(ierr);
#endif
  
  return 0;
}

PetscErrorCode ShelfModel::initFromOptions() {
  PetscErrorCode ierr;

  //ierr = initBoring(); CHKERRQ(ierr);
  ierr = initEmbayment(); CHKERRQ(ierr);
  //ierr = initFloatingSlab(); CHKERRQ(ierr);

  ierr = IceModel::initFromOptions(); CHKERRQ(ierr);
  
  return 0;
}

PetscErrorCode ShelfModel::initFloatingSlab() {
  PetscErrorCode ierr;
  const PetscScalar G_geothermal = 0.0;             // J/m^2 s; geo. heat flux
  const PetscScalar L = 500e3;
  const PetscScalar slab_thickness   = 800;
  const PetscScalar null_thickness   = DEFAULT_MINH_MACAYEAL;
  PetscScalar **h, **H, **mask;

  ierr = initIceParam(grid.com, &grid.p, &grid.bag); CHKERRQ(ierr);
  grid.p->Mbz = 0;
  ierr = grid.createDA(); CHKERRQ(ierr);
  ierr = createVecs(); CHKERRQ(ierr);
  ierr = grid.rescale(L, L, 1000); CHKERRQ(ierr);
  ierr = VecSet(vbed, -1e3); CHKERRQ(ierr);
  ierr = VecSet(vGhf, G_geothermal); CHKERRQ(ierr);
  ierr = VecSet(vAccum, 0); CHKERRQ(ierr);
  ierr = VecSet(vTs, 260); CHKERRQ(ierr);
  ierr = VecSet(vT, 260); CHKERRQ(ierr);
  ierr = VecSet(vtau, DEFAULT_INITIAL_AGE_YEARS); CHKERRQ(ierr);
  setConstantGrainSize(DEFAULT_GRAIN_SIZE);
  setInitialAgeYears(DEFAULT_INITIAL_AGE_YEARS);

  ierr = DAVecGetArray(grid.da2, vH, &H); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vh, &h); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vMask, &mask); CHKERRQ(ierr);
  for (PetscInt i = grid.xs; i < grid.xs + grid.xm; i++) {
    for (PetscInt j = grid.ys; j < grid.ys + grid.ym; j++) {
      const PetscScalar x = -grid.p->Lx + i * grid.p->dx;
      const PetscScalar y = -grid.p->Ly + j * grid.p->dy;
      const PetscScalar r = sqrt(PetscSqr(x) + PetscSqr(y));
      if (r < grid.p->Lx / 1.5)
        H[i][j] = slab_thickness;
      else
        H[i][j] = null_thickness;
      // H[i][j] = slab_thickness * PetscSqr(1 - r / (4 * L));
      if (false) {
        if (r < grid.p->dx * grid.p->Mx / 8)
          mask[i][j] = MASK_SHEET;
        else
          mask[i][j] = MASK_FLOATING;
      }
        
      h[i][j] = H[i][j] * (1 - ice.rho / ocean.rho);
    }
  }
  ierr = DAVecRestoreArray(grid.da2, vH, &H); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vh, &h); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vMask, &mask); CHKERRQ(ierr);  
  ierr = DALocalToLocalBegin(grid.da2, vH, INSERT_VALUES, vH); CHKERRQ(ierr);
  ierr = DALocalToLocalEnd(grid.da2, vH, INSERT_VALUES, vH); CHKERRQ(ierr);
  ierr = DALocalToLocalBegin(grid.da2, vh, INSERT_VALUES, vh); CHKERRQ(ierr);
  ierr = DALocalToLocalEnd(grid.da2, vh, INSERT_VALUES, vh); CHKERRQ(ierr);
  ierr = DALocalToLocalBegin(grid.da2, vMask, INSERT_VALUES, vMask); CHKERRQ(ierr);
  ierr = DALocalToLocalEnd(grid.da2, vMask, INSERT_VALUES, vMask); CHKERRQ(ierr);

  initialized_p = PETSC_TRUE;
  return 0;
}


PetscErrorCode ShelfModel::initEmbayment() {
  const PetscScalar G_geothermal = 0.0;
  const PetscScalar B = 135720960.0;
  // 135720880.827
  const PetscScalar L = 500e3;
  const PetscScalar scale_factor = 1.0;
  PetscErrorCode ierr;
  PetscScalar **modelH, **modelh, **modelb, **mask, **myNu[2];
  PetscScalar **uvbar[2];
  PetscScalar **u_exact, **v_exact;
  
  ierr = initIceParam(grid.com, &grid.p, &grid.bag); CHKERRQ(ierr);
  grid.p->Mbz = 0;
  ierr = grid.createDA(); CHKERRQ(ierr);
  ierr = createVecs(); CHKERRQ(ierr);
  ierr = grid.rescale(0.5 * scale_factor * L, scale_factor * L, 1000); CHKERRQ(ierr);
  ierr = VecSet(vbed, 0); CHKERRQ(ierr);
  ierr = VecSet(vMask, MASK_FLOATING); CHKERRQ(ierr);
  ierr = VecSet(vGhf, G_geothermal); CHKERRQ(ierr);
  ierr = VecSet(vAccum, 0); CHKERRQ(ierr);
  ierr = VecSet(vTs, 260); CHKERRQ(ierr);
  ierr = VecSet(vT, 260); CHKERRQ(ierr);
  ierr = VecSet(vtau, DEFAULT_INITIAL_AGE_YEARS); CHKERRQ(ierr);
  setConstantGrainSize(DEFAULT_GRAIN_SIZE);
  setInitialAgeYears(DEFAULT_INITIAL_AGE_YEARS);

  ierr = DAVecGetArray(grid.da2, vH, &modelH); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vh, &modelh); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vbed, &modelb); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vMask, &mask); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vbetax, &basalBetax); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vbetay, &basalBetay); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vNu[0], &myNu[0]); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vNu[1], &myNu[1]); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vuvbar[0], &uvbar[0]); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vuvbar[1], &uvbar[1]); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vu_exact, &u_exact); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vv_exact, &v_exact); CHKERRQ(ierr);

  for (PetscInt i = grid.xs; i < grid.xs + grid.xm; i++) {
    for (PetscInt j = grid.ys; j < grid.ys + grid.ym; j++) {
      const PetscScalar xx = i * grid.p->dx;
      const PetscScalar yy = (j - (grid.p->My - 1) / 2) * grid.p->dy;
      const PetscScalar x0 = xx / (2 * grid.p->Lx);
      const PetscScalar y0 = yy / grid.p->Ly;
      {
        const PetscScalar x = x0;
        const PetscScalar y = y0;
#include "embayment/test.h"
        const PetscScalar nu = (B / 2) / cbrt(alpha / (L * L));
        const PetscScalar nu_x = -nu * alpha_x / (3 * alpha);
        const PetscScalar nu_y = -nu * alpha_y / (3 * alpha);
        const PetscScalar rhsx = ((2 * nu_x * H * (2 * u_x + v_y)
                                   + 2 * nu * H_x * (2 * u_x + v_y)
                                   + 2 * nu * H * (2 * u_xx + v_xy)) / (L * L)
                                  + (nu_y * H * (u_y + v_x)
                                     + nu * H_y * (u_y + v_x)
                                     + nu * H * (u_yy + v_xy)) / (L * L)
                                  - IceType::rho * MaterialType::grav * H * h_x / L);
        const PetscScalar rhsy = ((2 * nu_y * H * (u_x + 2 * v_y)
                                   + 2 * nu * H_y * (u_x + 2 * v_y)
                                   + 2 * nu * H * (u_xy + 2 * v_yy)) / (L * L)
                                  + (nu_x * H * (u_y + v_x)
                                     + nu * H_x * (u_y + v_x)
                                     + nu * H * (u_xy + v_xx)) / (L * L)
                                  - IceType::rho * MaterialType::grav * H * h_y / L);
        modelH[i][j] = H;
        modelh[i][j] = h;
        modelb[i][j] = b;
#if 1
        basalBetax[i][j] = (u == 0) ? 0 : rhsx / u;
        basalBetay[i][j] = (v == 0) ? 0 : rhsy / v;
#else
        basalBetax[i][j] = rhsx;
        basalBetay[i][j] = rhsy;
#endif
        u_exact[i][j] = u;
        v_exact[i][j] = v;
      }

      for (PetscInt o = 0; o < 2; o++) {
        for (PetscInt op = 0; op < 2; op++) {
          const PetscInt oi = (o == 0 && op == 0) ? -1 : 0;
          const PetscInt oj = (o == 1 && op == 0) ? -1 : 0;
          //xx = (i - 0.5 + oi) * grid.p->dx;
          //yy = (j - 0.5 + oj - (grid.p->My - 1) / 2) * grid.p->dy;
          //x = xx / (2 * grid.p->Lx);
          //y = yy / grid.p->Ly;
          const PetscScalar dix = (o == 0) ? oi + 0.5 : oi;
          const PetscScalar diy = (o == 1) ? oj + 0.5 : oj;
          const PetscScalar x = x0 + dix / (grid.p->Mx - 1);
          const PetscScalar y = y0 + diy * 2 / (grid.p->My - 1);
          //#include "embayment/uv.h"
          //#include "embayment/nu.h"
#include "embayment/test.h"
          const PetscScalar nu = (B / 2) / cbrt(alpha / (L * L));
          myNu[o][i + oi][j + oj] = nu;
          //printf("[%4.3f, %4.3f] ", dix, diy);
          if (false) printf("(%4.3f+%4.3f, %4.3f+%4.3f; %2d, %2d) %e %e\n",
                            x0, x, y0, y, oi, oj, u * secpera, v * secpera);
          uvbar[o][i + oi][j + oj] = (o == 0) ? u : v;
        }
      }

      //if (0 < x0 && x0 < 1 && PetscAbs(y0) < 1 && true) {
      //if (0.1 < x0 && x0 < 0.8 && PetscAbs(y0) < 0.5 && true) {
      if (0.1 < x0 && x0 < 0.9 && PetscAbs(y0) < 0.9 && true) {
        mask[i][j] = MASK_DRAGGING;
      } else {
        mask[i][j] = MASK_SHEET;
      }
    }
  }

  ierr = DAVecRestoreArray(grid.da2, vH, &modelH); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vh, &modelh); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vbed, &modelb); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vMask, &mask); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vuvbar[0], &uvbar[0]); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vuvbar[1], &uvbar[1]); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vNu[0], &myNu[0]); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vNu[1], &myNu[1]); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vu_exact, &u_exact); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vv_exact, &v_exact); CHKERRQ(ierr);
  ierr = VecCopy(vu_exact, vubar); CHKERRQ(ierr);
  ierr = VecCopy(vv_exact, vvbar); CHKERRQ(ierr);
  // We will leave these open so we can use them in the model.
  // ierr = DAVecRestoreArray(grid.da2, vbetax, &basalBetax); CHKERRQ(ierr);
  // ierr = DAVecRestoreArray(grid.da2, vbetay, &basalBetay); CHKERRQ(ierr);

//   ierr = DAVecRestoreArray(grid.da2, vbetax, &basalBetax); CHKERRQ(ierr);
//   ierr = DAVecRestoreArray(grid.da2, vbetay, &basalBetay); CHKERRQ(ierr);
//   ierr = DALocalToGlobal(grid.da2, vbetax, INSERT_VALUES, g2); CHKERRQ(ierr);
//   ierr = VecView(g2, PETSC_VIEWER_DRAW_WORLD); CHKERRQ(ierr); 

  initialized_p = PETSC_TRUE;
  return 0;
}

PetscErrorCode ShelfModel::initBoring() {
  const PetscScalar G_geothermal = 0.0;
  const PetscScalar B = 135720960.0;
  // 135720880.827
  const PetscScalar L = 500e3;
  const PetscScalar scale_factor = 1.0;
  const PetscScalar thickness = 500.0;
  const PetscScalar velocity_grad = 200 / secpera;
  PetscErrorCode ierr;
  PetscScalar **mask, **myNu;
  PetscScalar **uvbar[2];
  PetscScalar **u_exact, **v_exact;
  
  ierr = initIceParam(grid.com, &grid.p, &grid.bag); CHKERRQ(ierr);
  grid.p->Mbz = 0;
  ierr = grid.createDA(); CHKERRQ(ierr);
  ierr = createVecs(); CHKERRQ(ierr);
  ierr = grid.rescale(scale_factor * L, scale_factor * L, 1000); CHKERRQ(ierr);
  ierr = VecSet(vbed, 0); CHKERRQ(ierr);
  ierr = VecSet(vMask, MASK_FLOATING); CHKERRQ(ierr);
  ierr = VecSet(vGhf, G_geothermal); CHKERRQ(ierr);
  ierr = VecSet(vAccum, 0); CHKERRQ(ierr);
  ierr = VecSet(vh, thickness); CHKERRQ(ierr);
  ierr = VecSet(vH, thickness); CHKERRQ(ierr);
  ierr = VecSet(vbed, 0); CHKERRQ(ierr);
  ierr = VecSet(vTs, 260); CHKERRQ(ierr);
  ierr = VecSet(vT, 260); CHKERRQ(ierr);
  ierr = VecSet(vtau, DEFAULT_INITIAL_AGE_YEARS); CHKERRQ(ierr);
  setConstantGrainSize(DEFAULT_GRAIN_SIZE);
  setInitialAgeYears(DEFAULT_INITIAL_AGE_YEARS);

  ierr = DAVecGetArray(grid.da2, vMask, &mask); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vbetax, &basalBetax); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vbetay, &basalBetay); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vNu[0], &myNu[0]); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vNu[1], &myNu[1]); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vuvbar[0], &uvbar[0]); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vuvbar[1], &uvbar[1]); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vu_exact, &u_exact); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vv_exact, &v_exact); CHKERRQ(ierr);

  for (PetscInt i = grid.xs; i < grid.xs + grid.xm; i++) {
    for (PetscInt j = grid.ys; j < grid.ys + grid.ym; j++) {
      const PetscScalar xx = (i - (grid.p->Mx - 1) / 2) * grid.p->dx;
      const PetscScalar yy = (j - (grid.p->My - 1) / 2) * grid.p->dy;
      const PetscScalar x0 = xx / grid.p->Lx;
      const PetscScalar y0 = yy / grid.p->Ly;
      {
        const PetscScalar x = x0;
        const PetscScalar y = y0;
        basalBetax[i][j] = 0;
        basalBetay[i][j] = 0;
        const PetscScalar u_x = velocity_grad / grid.p->Lx;
        const PetscScalar u_y = 0.0;
        const PetscScalar v_x = 0.0;
        const PetscScalar v_y = velocity_grad / grid.p->Ly;
        myNu[i][j] = thickness * 0.5 * B * pow(0.5 * PetscSqr(u_x) + 0.5 * PetscSqr(v_y)
                                               + 0.5 * PetscSqr(u_x + v_y)
                                               + 0.25 * PetscSqr(u_y + v_x), -1 / 3);
        u_exact[i][j] = x * velocity_grad;
        v_exact[i][j] = y * velocity_grad;
      }

      for (PetscInt o = 0; o < 2; o++) {
        for (PetscInt op = 0; op < 2; op++) {
          const PetscInt oi = (o == 0 && op == 0) ? -1 : 0;
          const PetscInt oj = (o == 1 && op == 0) ? -1 : 0;
          //xx = (i - 0.5 + oi) * grid.p->dx;
          //yy = (j - 0.5 + oj - (grid.p->My - 1) / 2) * grid.p->dy;
          //x = xx / (2 * grid.p->Lx);
          //y = yy / grid.p->Ly;
          const PetscScalar dix = (o == 0) ? oi + 0.5 : oi;
          const PetscScalar diy = (o == 1) ? oj + 0.5 : oj;
          const PetscScalar x = x0 + dix * 2 / (grid.p->Mx - 1);
          const PetscScalar y = y0 + diy * 2 / (grid.p->My - 1);
          const PetscScalar u = velocity_grad * x;
          const PetscScalar v = velocity_grad * y;
          uvbar[o][i + oi][j + oj] = (o == 0) ? u : v;
        }
      }

      if (PetscAbs(x0) < 0.8 && PetscAbs(y0) < 0.8 && true) {
        mask[i][j] = MASK_DRAGGING;
      } else {
        mask[i][j] = MASK_SHEET;
      }
    }
  }

  ierr = DAVecRestoreArray(grid.da2, vMask, &mask); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vuvbar[0], &uvbar[0]); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vuvbar[1], &uvbar[1]); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vNu[0], &myNu[0]); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vNu[1], &myNu[1]); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vu_exact, &u_exact); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vv_exact, &v_exact); CHKERRQ(ierr);
  
  if (false) {
    ierr = DAVecRestoreArray(grid.da2, vbetax, &basalBetax); CHKERRQ(ierr);
    ierr = DAVecRestoreArray(grid.da2, vbetay, &basalBetay); CHKERRQ(ierr);
  }
  
  initialized_p = PETSC_TRUE;
  return 0;
}


PetscErrorCode ShelfModel::afterInitHook() {
  PetscErrorCode ierr;
  ierr = IceModel::afterInitHook(); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode ShelfModel::run() {
  PetscErrorCode ierr;

  //ierr = IceModel::run(); CHKERRQ(ierr);
  ierr = VecSet(vubar, 0.0); CHKERRQ(ierr);
  ierr = VecSet(vvbar, 0.0); CHKERRQ(ierr);

  ierr = setupForMacayeal(DEFAULT_MINH_MACAYEAL,PETSC_FALSE); CHKERRQ(ierr);
//  ierr = mapStaggeredVelocityToStandard(); CHKERRQ(ierr);
  ierr = vertAveragedVelocityToRegular(); CHKERRQ(ierr);
  ierr = velocityMacayeal(); CHKERRQ(ierr);
  ierr = cleanupAfterMacayeal(DEFAULT_MINH_MACAYEAL); CHKERRQ(ierr);
  ierr = broadcastMacayealVelocity(); CHKERRQ(ierr);
  ierr = computeMaxVelocities(); CHKERRQ(ierr);
  ierr = massBalExplicitStep(); CHKERRQ(ierr);

  ierr = updateViewers(); CHKERRQ(ierr);
  
  {
    PetscViewer u_error_view, v_error_view;
    PetscReal norm_inf, norm_2, norm_1;
    ierr = PetscViewerDrawOpen(grid.com, PETSC_NULL, "u error",
                               PETSC_DECIDE, PETSC_DECIDE, PETSC_DECIDE, PETSC_DECIDE,
                               &u_error_view); CHKERRQ(ierr);
    ierr = PetscViewerDrawOpen(grid.com, PETSC_NULL, "v error",
                               PETSC_DECIDE, PETSC_DECIDE, PETSC_DECIDE, PETSC_DECIDE,
                               &v_error_view); CHKERRQ(ierr);
    
    ierr = VecWAXPY(vWork2d[0], -1, vu_exact, vubar); CHKERRQ(ierr);
    ierr = DALocalToGlobal(grid.da2, vWork2d[0], INSERT_VALUES, g2); CHKERRQ(ierr);
    ierr = VecScale(g2, secpera);
    ierr = VecView(g2, u_error_view); CHKERRQ(ierr);
    ierr = VecNorm(g2, NORM_INFINITY, &norm_inf); CHKERRQ(ierr);
    ierr = VecNorm(g2, NORM_2, &norm_2); CHKERRQ(ierr);
    ierr = VecNorm(g2, NORM_1, &norm_1); CHKERRQ(ierr);
    ierr = PetscPrintf(grid.com, "u : norm_infty |Error| = %e\n"
                       "u : norm_2 |Error| = %e\n"
                       "u : norm_1 |Error| = %e\n",
                       norm_inf,
                       norm_2 / (grid.p->Mx * grid.p->My),
                       norm_1 / (grid.p->Mx * grid.p->My)); CHKERRQ(ierr);

    ierr = VecWAXPY(vWork2d[0], -1, vv_exact, vvbar); CHKERRQ(ierr);
    ierr = DALocalToGlobal(grid.da2, vWork2d[0], INSERT_VALUES, g2); CHKERRQ(ierr);
    ierr = VecScale(g2, secpera);
    ierr = VecView(g2, v_error_view); CHKERRQ(ierr);
    ierr = VecNorm(g2, NORM_INFINITY, &norm_inf); CHKERRQ(ierr);
    ierr = VecNorm(g2, NORM_2, &norm_2); CHKERRQ(ierr);
    ierr = VecNorm(g2, NORM_1, &norm_1); CHKERRQ(ierr);
    ierr = PetscPrintf(grid.com, "v : norm_infty |Error| = %e\n"
                       "v : norm_2 |Error| = %e\n"
                       "v : norm_1 |Error| = %e\n",
                       norm_inf,
                       norm_2 / (grid.p->Mx * grid.p->My),
                       norm_1 / (grid.p->Mx * grid.p->My)); CHKERRQ(ierr);
  }
  
  return 0;
}

int main(int argc, char *argv[]) {
  PetscErrorCode ierr;
  MPI_Comm com;
  PetscMPIInt rank, size;
  
  ierr = PetscInitialize(&argc, &argv, PETSC_NULL, help); CHKERRQ(ierr);

  com = PETSC_COMM_WORLD;
  ierr = MPI_Comm_rank(com, &rank); CHKERRQ(ierr);
  ierr = MPI_Comm_size(com, &size); CHKERRQ(ierr);

  if (false) {
    const PetscScalar B = 135720960.;
    const PetscScalar L = 500e3;
    const PetscScalar x = 0.5, y = 0; //0.1e-5 * L;
    //#include "embayment/H.h"
    //#include "embayment/uv.h"
#include "embayment/test.h"
#include "embayment/betax.h"
    //#include "embayment/betay.h"
    ierr = PetscPrintf(com, "B = %e\n", B); CHKERRQ(ierr);
    ierr = PetscPrintf(com, "H(%f, %f) = %e m\n", x, y, H); CHKERRQ(ierr);
    ierr = PetscPrintf(com, "u(%f, %f) = %e %e %e m/a\n", x, y, u * secpera, u_x * secpera / L, u_xx * secpera / (L * L)); CHKERRQ(ierr);
    ierr = PetscPrintf(com, "v(%f, %f) = %e m/a\n", x, y, v * secpera); CHKERRQ(ierr);
    ierr = PetscPrintf(com, "alpha(%f, %f) = %e\n",
                       x, y, alpha / (L * L)); CHKERRQ(ierr);
    ierr = PetscPrintf(com, "alpha_x(%f, %f) = %e\n",
                       x, y, alpha_x / (L * L * L)); CHKERRQ(ierr);
    ierr = PetscPrintf(com, "TEST(%f, %f) = %e\n",
                       x, y, /*pow(L, 0.6667) */ (B / 2) / cbrt(alpha)); CHKERRQ(ierr);
    const PetscScalar nu = (B / 2) / cbrt(alpha / (L * L));
    const PetscScalar nu_x = -nu * alpha_x / (3 * alpha);
    const PetscScalar nu_y = -nu * alpha_y / (3 * alpha);
    const PetscScalar rhsx = ((2 * nu_x * H * (2 * u_x + v_y)
                               + 2 * nu * H_x * (2 * u_x + v_y)
                               + 2 * nu * H * (2 * u_xx + v_xy)) / (L * L)
                              + (nu_y * H * (u_y + v_x)
                                 + nu * H_y * (u_y + v_x)
                                 + nu * H * (u_yy + v_xy)) / (L * L)
                              - IceType::rho * MaterialType::grav * H * h_x / L);
    ierr = PetscPrintf(com, "nu: (%f, %f) %e %e %e\n", x, y, nu, nu_x/L, nu_y/L); CHKERRQ(ierr);
    ierr = PetscPrintf(com, "rhsx(%f, %f) = %e\n", x, y, rhsx); CHKERRQ(ierr);
    ierr = PetscPrintf(com, "my betax(%f, %f) = %e\n", x, y, rhsx / u); CHKERRQ(ierr);

    ierr = PetscPrintf(com, "betax(%f, %f) = %e\n", x, y, betax); CHKERRQ(ierr);
    //ierr = PetscPrintf(com, "betay(%f, %f) = %e\n", x, y, betay); CHKERRQ(ierr);
  }

  /* This explicit scoping forces destructors to be called before PetscFinalize() */
  {
    IceGrid g(com, rank, size);

    GlenIce  ice;
    //ThermoGlenArrIce ice;
    
    ShelfModel m(g, ice);

    ierr = m.setFromOptions(); CHKERRQ(ierr);
    ierr = m.initFromOptions(); CHKERRQ(ierr);
    ierr = m.setSoundingFromOptions(); CHKERRQ(ierr);

    ierr = PetscPrintf(com, "initialization done\n"); CHKERRQ(ierr);

    ierr = m.run(); CHKERRQ(ierr);

    {
      PetscInt pause_time;
      PetscTruth pause_p;
      ierr = PetscOptionsGetInt(PETSC_NULL, "-pause", &pause_time, &pause_p); CHKERRQ(ierr);
      if (pause_p == PETSC_TRUE)
        ierr = PetscSleep(pause_time); CHKERRQ(ierr);
    }

    ierr = PetscPrintf(com, "done with run ... "); CHKERRQ(ierr);
    // see comments in run_ice.cc re default output naming convention
    //ierr = m.writeFiles("shelf"); CHKERRQ(ierr);
    ierr = PetscPrintf(com, " ... done.\n"); CHKERRQ(ierr);
  }

  ierr = PetscFinalize(); CHKERRQ(ierr);
  return 0;
}
